library('coda')
library('data.table')
library('mvtnorm')
library('R2jags')
library('rjags')

TrueParVal <- function(sim.study.obj, par.to.save) {
  mu = sim.study.obj$sim.obj0$pars.baseline$Mu
  ss_tpr = sim.study.obj$sim.pars$ss.tpr
  bs_tpr = sim.study.obj$sim.pars$bs.tpr
  bs_fpr = sim.study.obj$sim.pars$bs.fpr
  K = sim.study.obj$sim.par$K
  par.name = as.vector(sapply(par.to.save, function(x) {
    paste0(x, "[", 1:K, "]")}))
  par.est = as.vector(sapply(par.to.save, function(x) eval(parse(text = x))))
  true.par = data.frame(par.name = par.name,
                        par.est = par.est)
  true.par
}


WideToLong <- function(fit, method.name) {
  fit = as.data.frame(fit)
  var.names = colnames(fit)
  par.names = rep(var.names, each = nrow(fit))
  fit.value = c()
  for (i in 1:ncol(fit)) {
    fit.value = c(fit.value, fit[, i])
  }
  result = data.frame(Estimate = fit.value,
                      Method = rep(method.name, prod(dim(fit))),
                      Parameter = par.names)
  result
}


ListEtiologyPriorSC1 <- function(K, Smax, hyper.pars.list,
                                 patho.names,
                                 n.sample = 5000, round.digit = 4,
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
                       sqrt(1/hyper.pars.list$tau_theta))
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
  cell.prob.upper = apply(cell.prob, 2, quantile, probs = 0.975)
  cell.prob.lower = apply(cell.prob, 2, quantile, probs = 0.025)
  LMAT = rbind(rep(0, n.theta1), design.mat$Lmat)
  s.seq = rowSums(LMAT)
  aggr.upper = aggregate(cell.prob.upper, by = list(s.seq), mean)$x
  aggr.lower = aggregate(cell.prob.lower, by = list(s.seq), mean)$x
  for (i in 1:length(s.seq)) {
    cell.prob.upper[i] = aggr.upper[s.seq[i] + 1]
    cell.prob.lower[i] = aggr.lower[s.seq[i] + 1]
  }
  
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
                             decreasing = TRUE), round.digit),
    Prob.lower = round(cell.prob.lower[order(cell.prob.mean,
                                             decreasing = TRUE)], 4),
    Prob.upper = round(cell.prob.upper[order(cell.prob.mean,
                                             decreasing = TRUE)], 4)    
    )[1:num.keep, ]
  rownames(EtioCombProb) = NULL
  
  Pr.num = cbind(cell.prob[, 1] ,
                 cell.prob[, -1] %*% t(design.mat$PiMat))
  Pr.num.mean = colMeans(Pr.num)
  Pr.num.pathogen = cbind(Pr.num.mean, t(apply(Pr.num, 2, quantile,
                                               c(0.025, 0.975))))
  rownames(Pr.num.pathogen) = paste0("S=", 0:Smax)
  colnames(Pr.num.pathogen) = c("mean", "lower", "upper")
  
  list(Cell.prob = subset(EtioCombProb, Probability >= threshold),
       Pr.num.pathogen = round(Pr.num.pathogen, round.digit),
       EtioMat = EtioMat)
}

VarLogitPriorSC1 <- function(K, Smax, hyper.pars.list,
                             lambda = 0.1,
                             n.sample = 20000) {
# parameters include: pi0, mu1-K, ss_tpr, bs_tpr, bs_fpr 
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  n.theta1 = ncol(design.mat$Lmat)
  n.theta2 = ncol(design.mat$Umat)
  # generate theta1 and theta2
  theta1 = rmvnorm(n = n.sample, mean = hyper.pars.list$mu_logit,
                   sigma = diag(rep(1/hyper.pars.list$tau_theta, n.theta1)))
  pind = rbeta(n.sample, hyper.pars.list$pind_a, hyper.pars.list$pind_b)
  theta2.value = rnorm(n.sample * n.theta2,
                       hyper.pars.list$theta2_mu,
                       sqrt(1/hyper.pars.list$tau_theta))
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
  
  Pi0 = Logit(cell.prob[, 1])
  Mu = Logit(cell.prob[, -1] %*% t(design.mat$MuMat))
  ss_tpr = Logit(rbeta(n.sample, hyper.pars.list$aa, hyper.pars.list$bb))
  bs_tpr = Logit(rbeta(n.sample, hyper.pars.list$cc, hyper.pars.list$dd))
  bs_fpr = Logit(rbeta(n.sample, hyper.pars.list$ee, hyper.pars.list$ff))
  
  v.prior = matrix(0, nrow = 4 * K + 1, ncol = 4 * K + 1)
  v1 = covRob(data = cbind(Pi0, Mu))$cov
  v2 = soft.thresholding(v1, lambda)
  diag(v2) = diag(v1)
  v.prior[1:(K + 1), 1:(K + 1)] = v2
  diag(v.prior)[-(1:(K + 1))] = c(rep(var(bs_tpr), K),
                                  rep(var(ss_tpr), K),
                                  rep(var(bs_fpr), K))
  v.prior
}

ExtractMu <- function(coda.chains, sim.obj, print.summary = FALSE) {
  # This function uses the posterior samples of theta1 and theta2 in __non-reg__ 
  # version to compute the marginal mean of the latent variable (etiology pi).
  K = ncol(sim.obj$L)
  Smax = sim.obj$Smax
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  varnames = colnames(coda.chains)
  D = ncol(sim.obj$X)
  if (D == 1) {
    cell_prob = coda.chains[, paste0("cell_prob[", 1:(design.mat$J1 + 1), "]")]
    A = apply(cell_prob, 1, sum)
    cell_prob = cell_prob/A
    Mu.postsample = cell_prob[, -1] %*% t(design.mat$MuMat)
    if (print.summary) {
      print(summary(Mu.postsample))
    }
    return(Mu.postsample)
  } else {
    X.unique = uniquecombs(sim.obj$X)
    N.strata = nrow(X.unique)
    Mu.postsample.list = list()
    for (s in 1:N.strata) {
      cell_prob = coda.chains[, paste0("cell_prob[",
                                       s, ",",
                                       1:(design.mat$J1 + 1), "]")]
      A = apply(cell_prob, 1, sum)
      cell_prob = cell_prob/A
      Mu.postsample = cell_prob[, -1] %*% t(design.mat$MuMat)
      Mu.postsample.list[[s]] = Mu.postsample
    }
    return(Mu.postsample.list)
  }
}

ExtractPrNumPath <- function(coda.chains, sim.obj, print.summary = FALSE) {
  # This function uses the posterior samples of theta1 and theta2 in __non-reg__ 
  # version to compute the marginal mean of the latent variable (etiology pi).
  K = ncol(sim.obj$L)
  Smax = sim.obj$Smax
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  varnames = colnames(coda.chains)
  D = ncol(sim.obj$X)
  if (D == 1) {
    cell_prob = coda.chains[, paste0("cell_prob[", 1:(design.mat$J1 + 1), "]")]
    A = apply(cell_prob, 1, sum)
    cell_prob = cell_prob/A
    Pi.postsample = cell_prob[, -1] %*% t(design.mat$PiMat)
    if (print.summary) {
      print(summary(Pi.postsample))
    }
    return(Pi.postsample)
  } else {
    X.unique = uniquecombs(sim.obj$X)
    N.strata = nrow(X.unique)
    Pi.postsample.list = list()
    for (s in 1:N.strata) {
      cell_prob = coda.chains[, paste0("cell_prob[",
                                       s, ",",
                                       1:(design.mat$J1 + 1), "]")]
      A = apply(cell_prob, 1, sum)
      cell_prob = cell_prob/A
      Pi.postsample = cell_prob[, -1] %*% t(design.mat$PiMat)
      Pi.postsample.list[[s]] = Pi.postsample
    }
    return(Pi.postsample.list)
  }
}

ListEtiology <- function(coda.chains, sim.obj, etio.names, reorder = TRUE,
                         num.keep = NULL) {
  # This function uses the posterior samples of theta1 and theta2 in __non-reg__
  # version to compute the probability of each combination of pathogens.
  K = ncol(sim.obj$L)
  Smax = sim.obj$Smax
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  EtioList = list()
  n.strata = nrow(uniquecombs(sim.obj$X))
  if (length(n.strata) < 1) {
    n.strata = 1
  }
  for (i in 1:n.strata) {
    if (n.strata == 1) {
      cell_prob = coda.chains[, paste0("cell_prob[",
                                       1:(design.mat$J1 + 1), "]")]
    } else {
      cell_prob = coda.chains[, paste0("cell_prob[", i, ",",
                                       1:(design.mat$J1 + 1), "]")]
    }
    A = apply(cell_prob, 1, sum)
    cell_prob = cell_prob/A
    cell.prob.mean = apply(cell_prob, 2, mean)
    cell.prob.upper = apply(cell_prob, 2, quantile, probs = 0.975)
    cell.prob.lower = apply(cell_prob, 2, quantile, probs = 0.025)
    
    prob.therest = NULL
    if (length(num.keep) > 0 & num.keep < ncol(cell_prob)) {
      prob.therest = rowSums(cell_prob[, -(1:num.keep)])
      prob.therest.mean = mean(prob.therest)
      prob.therest.upper = quantile(prob.therest, probs = 0.975)
      prob.therest.lower = quantile(prob.therest, probs = 0.025)
    }
  
    LMAT = rbind(rep(0, K), design.mat$Lmat)
    if (reorder) {
      EtioMat = LMAT[order(cell.prob.mean, decreasing = TRUE), ]
    } else {
      EtioMat = LMAT  
    }
    EtioComb = apply(EtioMat, 1, function(x) {
      paste(etio.names[x > 0], collapse = "-")
    })
    EtioComb[EtioComb == ""] = "None_Above"
    if (length(num.keep) < 1) {
      num.keep = length(EtioComb)
    }
    if (reorder) {
      EtioCombProb = data.frame(
        EtioComb = EtioComb,
        Probability = round(sort(cell.prob.mean,
                                 decreasing = TRUE), 4),
        Prob.lower = round(cell.prob.lower[order(cell.prob.mean,
                                                 decreasing = TRUE)], 4),
        Prob.upper = round(cell.prob.upper[order(cell.prob.mean,
                                                 decreasing = TRUE)], 4)
      )[1:num.keep, ]      
    } else {
      EtioCombProb = data.frame(
        EtioComb = EtioComb,
        Probability = round(cell.prob.mean, 4),
        Prob.lower = round(cell.prob.lower, 4),
        Prob.upper = round(cell.prob.upper, 4)
      )[1:num.keep, ]      
    }
    if (length(prob.therest) > 0) {
      EtioCombProb = rbind(EtioCombProb,
                           data.frame(EtioComb = "The_Rest",
                                      Probability = prob.therest.mean,
                                      Prob.lower = prob.therest.lower,
                                      Prob.upper = prob.therest.upper))
    }
    rownames(EtioCombProb) = NULL
    EtioList[[i]] = EtioCombProb
  }
  EtioList
}

SweepPosterior <- function(coda.chains, sim.obj, sweep.on = "ss_tpr",
                           model.type = "regression") {
# return residuals given variable sweep.on
  if (model.type == "regression") {
    par.mat = cbind(as.data.frame(
                      coda.chains[, grepl("tpr", colnames(coda.chains))]),
                    as.data.frame(
                      coda.chains[, grepl("fpr", colnames(coda.chains))]))
    par.mat.tr = cbind(as.data.frame(
                         coda.chains[, grepl("Beta", colnames(coda.chains))]),
                       Logit(par.mat))
  } else {
    Mu.fit = ExtractMu(coda.chains, sim.obj)
    colnames(Mu.fit) = paste0("Mu[", 1:ncol(Mu.fit), "]")
    Pr.NumPath.fit = ExtractPrNumPath(coda.chains, sim.obj)
    Pr.NoneAbove = apply(Pr.NumPath.fit, 1, function(x) 1 - sum(x))
  
    par.mat = cbind(Pr.NoneAbove, Mu.fit,
                    as.data.frame(
                      coda.chains[, grepl("tpr", colnames(coda.chains))]),
                    as.data.frame(
                      coda.chains[, grepl("fpr", colnames(coda.chains))]))
    par.mat.tr = Logit(par.mat)
  }
  if ("theta2_value" %in% colnames(coda.chains)) {
    par.mat.tr = cbind(par.mat.tr, coda.chains[, "theta2_value"])
    colnames(par.mat.tr)[ncol(par.mat.tr)] = "theta2.value"
  }
  if (length(sweep.on) > 0) {
    sweep.index = which(grepl(sweep.on, colnames(par.mat.tr)))
    par.a = par.mat.tr[, -sweep.index]
    par.b = par.mat.tr[, sweep.index]
    for (j in 1:ncol(par.a)) {
      y = par.a[, j]
      fit = lm(y ~ ., as.data.frame(cbind(y, par.b)))
      par.a[, j] = residuals(fit)
    }
    return(par.a)
  } else {
    return(par.mat.tr)
  }
}

ExtractMuCellProb <- function(sim.fit) {
  # arg:
  #   sim.fit: a data.frome generated by simulation study
  Mu.fit = list() 
  cell.prob.fit = NULL
  n.strata = 4
  for (i in 1:n.strata) {
    temp = sim.fit[, paste0(paste0("Mu.", i), ".", 1:5, ".")]
    colnames(temp) = paste0("Mu[", 1:5, "]")
    Mu.fit[[i]] = temp
    # sim.fit[, paste0(paste0("PrNumPath.", i), ".", 0:5, ".")]
    temp = sim.fit[, paste0(paste0("cell_prob.", i), ".", 1:32, ".")]
    colnames(temp) = paste0(paste0("cell_prob[", i), ",", 1:32, "]")
    if (length(cell.prob.fit) > 0) {
      cell.prob.fit = cbind(cell.prob.fit, temp)    
    } else {
      cell.prob.fit = temp
    }
  }
  list(Mu.fit = Mu.fit, cell.prob.fit = cell.prob.fit)
}