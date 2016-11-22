library('coda')
library('data.table')
library('mvtnorm')
library('R2jags')
library('rjags')


SetDefaultHyperParameters <- function(K = 5, logit.mu.rho = -0.22) {
  V = matrix(logit.mu.rho, nrow = K, ncol = K)
  diag(V) = 1
  V = 1 * V
  list(lambda = 1.4,
       alpha = rep(1, K),
       mu_logit = rep(-0.5, K), V = V,
       ma = 2, mb = 2, 
       aa = 2, bb = 18,
       cc = 5, dd = 3,
       ee = 1, ff = 1,
       tau_theta = 0.2,
       pind_a = 3, pind_b = 2,
       pind_dirich_par = c(5, 1, 5),
       theta2_cat = c(0, 0.5, -1),
       theta2_mu = -1.5)
}

PrepareData <- function(sim.obj, hyper.pars.list) {
  K = ncol(sim.obj$L)
  Smax = sim.obj$Smax
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  design.mat$GSmat = rbind(rep(0, K), design.mat$Lmat)
  design.mat$Umat_ncol = ncol(design.mat$Umat)
  # new var: I_strata[], N_strata, Beta[,], X_unique[,], D
  X_unique = uniquecombs(sim.obj$X)
  sim.dat = c(list(
    # measurements
    N_ctrl = nrow(sim.obj$MBS.ctrl),
    K = ncol(sim.obj$MBS.ctrl),
    mbs_ctrl = sim.obj$MBS.ctrl,
    N_case = nrow(sim.obj$MBS.case),
    mbs_case = sim.obj$MBS.case,
    mss_case = sim.obj$MSS.case,
    # covariates
    X_unique = X_unique,
    I_strata = attr(X_unique, "index"),
    N_strata = nrow(X_unique),
    D = ncol(X_unique)),
    # hyper parameters
    hyper.pars.list,
    # design matrix for multivariate binary data
    design.mat
  )
  sim.dat
}

FitJags <- function(sim.obj, parameters.to.save,
                    model.file, bs.tpr.prefix = NULL,
                    hyper.pars.list = SetDefaultHyperParameters(),
                    n.iter = 6000, n.burnin = 4000,
                    n.thin = 2, n.chains = 2) {
  # Fit JAGS non-regression models
  # Example:
  #   model.file = "./jags/SparseCorr1_BSandSS_NoReg.jags"
  #   par.default = SetDefaultSimulationParameter(5)
  #   sim.obj = do.call(SimulatePerchData, par.default)
  #   jags.result = FitJags(sim.obj, c("cell_prob", "bs_tpr", "ss_tpr"),
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


