library('coda')
library('data.table')
library('mvtnorm')
library('R2jags')
library('rjags')
library('rstan')

SetDefaultHyperParameters <- function(K = 5, logit.mu.rho = -0.22) {
  V = matrix(logit.mu.rho, nrow = K, ncol = K)
  diag(V) = 1
  V = 1 * V
  list(lambda = 1.4,
       alpha = rep(1, K),
       mu_logit = rep(-0.5, K), V = V,
       ma = 2, mb = 2, 
       aa = 1, bb = 9,
       cc = 4, dd = 2,
       ee = 1, ff = 1,
       tau_theta = 1/10,
       pind_a = 2, pind_b = 8,
       pind_dirich_par = c(5, 1, 5),
       theta2_cat = c(0, 0.5, -1),
       theta2_mu = -3)
}

PrepareData <- function(sim.obj, hyper.pars.list) {
  K = length(sim.obj$pars.baseline$Mu)
  Smax = length(sim.obj$pars.baseline$Pi) - 1
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  design.mat$GSmat = rbind(rep(0, K), design.mat$Lmat)
  design.mat$Umat_ncol = ncol(design.mat$Umat)
  sim.dat = c(list(
    # data
    N_ctrl = nrow(sim.obj$MBS.ctrl),
    K = ncol(sim.obj$MBS.ctrl),
    mbs_ctrl = sim.obj$MBS.ctrl,
    N_case = nrow(sim.obj$MBS.case),
    mbs_case = sim.obj$MBS.case,
    mss_case = sim.obj$MSS.case),
    # hyper parameters
    hyper.pars.list,
    # design matrix for multivariate binary data
    design.mat
  )
  sim.dat
}

FitJagsNoReg <- function(sim.obj, parameters.to.save,
                         model.file, bs.tpr.prefix = NULL,
                         hyper.pars.list = SetDefaultHyperParameters(),
                         n.iter = 6000, n.burnin = 4000,
                         n.thin = 2, n.chains = 2) {
# Fit JAGS non-regression models
# Example:
#   model.file = "./jags/SparseCorr1_BSandSS_NoReg.jags"
#   par.default = SetDefaultSimulationParameter(5)
#   sim.obj = do.call(SimulatePerchData, par.default)
#   jags.result = FitJagsNoReg(sim.obj, c("cell_prob", "bs_tpr", "ss_tpr"),
#     model.file, NULL, SetDefaultHyperParameters(K = 3), 1000, 500, 1, 1)
  sim.dat = PrepareData(sim.obj, hyper.pars.list)
  if (n.chains > 1) {
    mc.fit <- jags.parallel(
      data = sim.dat,
      #inits = bayes.mod.inits,
      parameters.to.save = parameters.to.save,
      n.chains = n.chains,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      model.file = model.file)
  } else {
    mc.fit <- jags(
      data = sim.dat,
      #inits = bayes.mod.inits,
      parameters.to.save = parameters.to.save,
      n.chains = n.chains,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      model.file = model.file)    
  }
  print(mc.fit)
  return(mc.fit)
}

ExtractMu <- function(coda.chains, sim.obj, print.summary = FALSE) {
# This function uses the posterior samples of theta1 and theta2 in __non-reg__ 
# version to compute the marginal mean of the latent variable (etiology pi).
  K = length(sim.obj$pars.baseline$Mu)
  Smax = length(sim.obj$pars.baseline$Pi) - 1
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  varnames = colnames(coda.chains)
  cell_prob = coda.chains[, paste0("cell_prob[", 1:(design.mat$J1 + 1), "]")]
  A = apply(cell_prob, 1, sum)
  cell_prob = cell_prob/A
  Mu.postsample = cell_prob[, -1] %*% t(design.mat$MuMat)
  if (print.summary) {
    print(summary(Mu.postsample))
  }
  Mu.postsample
}

ListEtiology <- function(coda.chains, sim.obj, top5.names,
                         num.keep = NULL, threshold = 0) {
# This function uses the posterior samples of theta1 and theta2 in __non-reg__
# version to compute the probability of each combination of pathogens.
  K = length(sim.obj$pars.baseline$Mu)
  Smax = length(sim.obj$pars.baseline$Pi) - 1
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  cell_prob = coda.chains[, paste0("cell_prob[", 1:(design.mat$J1 + 1), "]")]
  A = apply(cell_prob, 1, sum)
  cell_prob = cell_prob/A
  cell.prob.mean = apply(cell_prob, 2, mean)
  LMAT = rbind(rep(0, 5), design.mat$Lmat)
  EtioMat = LMAT[order(cell.prob.mean, decreasing = TRUE), ]
  EtioComb = apply(EtioMat, 1, function(x) {
    paste(top5.names[x > 0], collapse = "-")
  })
  EtioComb[EtioComb == ""] = "None_Above"
  if (length(num.keep) < 1) {
    num.keep = length(EtioComb)
  }
  EtioCombProb = data.frame(
      EtioComb = EtioComb,
      Probability = round(sort(cell.prob.mean,
                               decreasing = TRUE), 4))[1:num.keep, ]
  rownames(EtioCombProb) = NULL
  subset(EtioCombProb, Probability >= threshold)
}


ListEtiologyPriorSC1 <- function(K, Smax, hyper.pars.list,
                             patho.names,
                             n.sample = 5000,
                             num.keep = NULL, threshold = 0) {
# This function list the Etiology prior of each pathogen combination in the
# Sparse Corr 1 model.
# Example:
#   hyper.pars.list$theta2_mu = -1
#   hyper.pars.list$pind_a = 1
#   ListEtiologyPriorSC1(5, 5, hyper.pars.list, patho.names, 20000)
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  n.theta1 = ncol(design.mat$Lmat)
  n.theta2 = ncol(design.mat$Umat)
  # generate theta1 and theta2
  theta1 = rmvnorm(n = n.sample, mean = hyper.pars.list$mu_logit,
                   sigma = diag(rep(1/hyper.pars.list$tau_theta, n.theta1)))
  pind = rbeta(n.sample, hyper.pars.list$pind_a, hyper.pars.list$pind_b)
  theta2.value = rnorm(n.sample * n.theta2,
                       hyper.pars.list$theta2_mu,
                       1/hyper.pars.list$tau_theta)
  theta2.value = matrix(theta2.value, nrow = n.sample, ncol = n.theta2)
  theta2 = matrix(NA, nrow = n.sample, ncol = n.theta2)
  for (i in 1:n.sample) {
    ind = rbinom(n.theta2, 1, pind[i])
    theta2[i, ] = theta2.value[i, ] * ind
  }
  # calculate cell probabilities
  exp.theta = exp(cbind(0, cbind(theta1, theta2) %*%
                           t(cbind(design.mat$Lmat, design.mat$Umat))))
  
  A = apply(exp.theta, 1, sum)
  cell.prob = exp.theta/A
  cell.prob.mean = apply(cell.prob, 2, mean)
  LMAT = rbind(rep(0, n.theta1), design.mat$Lmat)
  EtioMat = LMAT[order(cell.prob.mean, decreasing = TRUE), ]
  EtioComb = apply(EtioMat, 1, function(x) {
    paste(patho.names[x > 0], collapse = "-")
  })
  EtioComb[EtioComb == ""] = "None_Above"
  if (length(num.keep) < 1) {
    num.keep = length(EtioComb)
  }
  EtioCombProb = data.frame(
    EtioComb = EtioComb,
    Probability = round(sort(cell.prob.mean,
                             decreasing = TRUE), 4))[1:num.keep, ]
  rownames(EtioCombProb) = NULL
  subset(EtioCombProb, Probability >= threshold)  
}
