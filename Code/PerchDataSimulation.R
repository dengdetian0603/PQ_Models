library(BB)
library(doMC)
library(foreach)
library(limSolve)
library(mgcv)
library(nleqslv)



DesignMatrixAppxQuadExp <- function(K, Smax = 3) {
# Generate design matrix based on approximated quadratic exponential family
#
# Example:
#   DesignMatrixAppxQuadExp(5, 3)
#
  index = cumsum(sapply(1:Smax, function(x) choose(K, x)))
  Lmat = matrix(NA, ncol = K, nrow = index[Smax])
  row1 = 1
  for (s in 1:Smax){
    Lmat[row1:index[s], ] =
        t(combn(1:K, s, function(x) {
          a = rep(0, K)
          a[x] = 1
          a
        }))
    row1 = index[s]+1
  }
  Umat = t(apply(Lmat, 1, function(l) apply(combn(l, 2), 2, prod)))
  MuMat = matrix(NA, nrow = K, ncol = index[Smax])
  for (i in 1:K) {
    MuMat[i,] = as.numeric(Lmat[, i] == 1)
  }
  PiMat = matrix(NA, nrow = Smax, ncol = index[Smax])
  S = apply(Lmat, 1, sum)
  for (i in 1:Smax) {
    PiMat[i, ] = as.numeric(S == i)
  }
  return(list(Lmat = Lmat, Umat = Umat, MuMat = MuMat,
              J1 = index[Smax], PiMat = PiMat))
}

GenDefaultCovariates <- function(nsample, num.covariates, interact = FALSE) {
# Simulate covariates, assuming all X are categorical.
# Example:
#   GenDefaultCovariates(10, 1)
  x.cat = list()
  for (i in 1:num.covariates) {
    x.cat[[i]] = c(0,1)
  }
  x.comb = do.call(expand.grid, x.cat)
  x.interact = apply(x.comb, 1, prod)
  if (num.covariates == 2 && interact) {
    x.unique = cbind(1, x.comb, x.interact)
    colnames(x.unique) = c("intercept", "AGE", "SEVERITY", "AGE:SEVERITY")
  } else {
    x.unique = cbind(1, x.comb)
  }
  nx.unique = nrow(x.unique)
  x.index = sample(1:nx.unique, nsample - nx.unique, replace = TRUE)
  X = as.matrix(rbind(x.unique, x.unique[x.index, ]))
  X
}

Logit <- function(x)
{
  log(x/(1 - x))
}

InvLogit <- function(x)
{
      1/(1 + exp(-x))
}

BetaToMu <- function(beta.mat, x.mat, K, nsample) {
# Example:
#   beta.mat = c(1, 0, 0.3)
#   x.mat = rep(1, 10)
#   BetaToMu(beta.mat, x.mat, 3, 10)
# 
  beta.mat = matrix(beta.mat, ncol = K)
  x.mat = matrix(x.mat, nrow = nsample)
  InvLogit(x.mat %*% beta.mat)
}

PiToBetaTheta <- function(K, Smax, Pi.seed, phi.iter = 500) {
# Example:
#   Pi = c(0.05, 0.5, 0.35, 0.1-1e-2-1e-3, 1e-2, 1e-3)
#   PiToBetaTheta(5, 5, Pi)
  Pi = Pi.seed[1:(1 + Smax)]
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  Phis = suppressWarnings(
      xsample(E = design.mat$PiMat, F = Pi[-1]/Pi[1],
              G = diag(rep(1, design.mat$J1)),
              H = rep(0, design.mat$J1),
              iter = phi.iter, burninlength = 5,
              type = "mirror", test = FALSE)$X)
  phis = Phis[phi.iter - 10, ]
  Thetas = solve(qr(cbind(design.mat$Lmat, design.mat$Umat),
                    LAPACK = TRUE), log(phis))
  phis.ls = exp(cbind(design.mat$Lmat, design.mat$Umat) %*% Thetas)
  cell.probs = phis.ls/(1 + sum(phis.ls))
  Pi0 = 1/(1 + sum(phis.ls))
  Mu = as.vector(design.mat$MuMat %*% cell.probs)
  Pi1toSmax = design.mat$PiMat %*% cell.probs
  list(Beta0 = Logit(Mu), Mu = Mu, Pi = c(Pi0, Pi1toSmax), Thetas = Thetas)
}


MuToTheta1 <- function(Mu, Theta2, LUmat, MuMat, K, initvalue=NULL) {
# 
# Example:
#   Pi = c(0.1, 0.5, 0.3, 0.1, 0, 0)
#   dmat = DesignMatrixAppxQuadExp(5, 3)
#   pars = PiToBetaTheta(5, 3, Pi)
#   MuToTheta1(pars$Mu, pars$Thetas[-(1:5)],
#              cbind(dmat$Lmat, dmat$Umat), dmat$MuMat, 5)
#
  Eqt <- function(theta1, mu = Mu, theta2 = Theta2) {
    theta = matrix(c(theta1, theta2), ncol = 1)
    # potentials denotes the un-normalized probability
    potentials = exp(LUmat %*% theta)
    # A is the normalizing constant, i.e. sum of all potentials
    A = sum(potentials) + 1
    # k values of mu
    values = as.vector(MuMat %*% potentials) - mu * A
    values
  }
  if (length(initvalue) < K) initvalue = rep(0, K)
  result = tryCatch(
      nleqslv(initvalue, fn = Eqt, method = "Broyden",
              global = "dbldog", xscalm = "auto", 
              control = list(maxit = 3500, cndtol = 1e-15, ftol = 1e-10))$x, 
      error = function(c) "error solving for theta1 by nleqslv.")
  if (is.character(result)) {
    print(result)
    result = tryCatch(
        BBsolve(par = initvalue, fn = Eqt, quiet = TRUE,
                control = list(NM = FALSE))$par, 
        error = function(c) "error solving for theta1 by BBsolve." )
    if (is.character(result)) {
      stop(result)
    }
  }
  #print(result$message)
  return(result)
}

 
XBetaToCellProb <- function(K, X.unique, Beta, Theta2, LUmat, MuMat,
                            theta1.init = NULL, normalize = TRUE) {
# 
# Example:
#   Pi = c(0.1, 0.5, 0.3, 0.1, 0, 0)
#   dmat = DesignMatrixAppxQuadExp(5, 3)
#   pars = PiToBetaTheta(5, 3, Pi)
#   X = unique(GenDefaultCovariates(10, 1))
#   beta = matrix(c(pars$Beta0, 0, 0, 0.1, -0.1, 0.2), byrow = TRUE, nrow = 2)
#   XBetaToCellProb(5, X, beta, pars$Thetas[-(1:5)],
#                   cbind(dmat$Lmat, dmat$Umat), dmat$MuMat)
  X.index = attr(X.unique, "index")
  Beta = matrix(Beta, ncol = K)
  Mu.unique = InvLogit(X.unique %*% Beta)
  theta1.unique = t(apply(Mu.unique, 1,
                    function(x) {
                      MuToTheta1(Mu = x, Theta2 = Theta2, LUmat, MuMat, K = K,
                                 initvalue = theta1.init)}))
  if (normalize == TRUE) {
  # cell probabilities
    cellprobs.unique =
        t(apply(theta1.unique, 1,
                function(x) {
                  c(1, exp(LUmat %*% c(x, Theta2)))/
                       (1 + sum(exp(LUmat %*% c(x, Theta2))))}))
  } else if (normalize == FALSE) {
  # unnormalized cell probabilities
    cellprobs.unique =
        t(apply(theta1.unique, 1,
                function(x) {c(1, exp(LUmat %*% c(x, Theta2)))}))
  } else {
  # log unnormailized cell probabilities
    cellprobs.unique =
        t(apply(theta1.unique, 1, function(x) {c(0, LUmat %*% c(x, Theta2))}))
  }
  return(cellprobs.unique)
}

XBetaToCellProb2 <- function(K, X.unique, Beta, theta2,
                             lumat, normalize = TRUE){
  X.index = attr(X.unique, "index")
  Beta = matrix(Beta, ncol = K)
  theta1.unique = X.unique %*% Beta
  
  if (normalize == TRUE) {
    cellprobs.unique = 
        t(apply(theta1.unique, 1,
                function(x) {
                  c(1, exp(lumat %*% c(x, theta2))) /
                       (1 + sum(exp(lumat %*% c(x, theta2))))
                }))
  } else if (normalize == FALSE) {
    cellprobs.unique =
        t(apply(theta1.unique, 1,
                function(x) {
                  c(1, exp(lumat %*% c(x, theta2)))
                }))
  } else {
    cellprobs.unique = t(apply(theta1.unique, 1,
                               function(x) {c(0, lumat %*% c(x, theta2))}))
  }
  return(cellprobs.unique)
}

SimulateGSdata <- function(ncase, K, Smax, Pi.seed,
                           num.covariates, Betas, phi.iter = 500) {
#
# Args:
#   Betas: matrix, nrow = num.covariates, ncol = K
#
# Example:
#   Pi = c(0.1, 0.5, 0.3, 0.1, 0, 0)
#   SimulateGSdata(10, 5, 3, Pi, 1, c(0, 0, 0.1, -0.1, 0.2))
  X = GenDefaultCovariates(ncase, num.covariates)
  dmat = DesignMatrixAppxQuadExp(K, Smax)
  pars = PiToBetaTheta(K, Smax, Pi.seed, phi.iter)
  if (length(Betas) > 0) Betas = matrix(Betas, nrow = num.covariates, ncol = K)
  Betas = rbind(pars$Beta0, Betas)
  X.unique = uniquecombs(X)
  X.index = attr(X.unique, "index")
  cell.prob.unique = XBetaToCellProb(K, X.unique, Betas, pars$Thetas[-(1:K)],
                                     cbind(dmat$Lmat, dmat$Umat), dmat$MuMat)
  Lmat.withZero = rbind(rep(0, K), dmat$Lmat)
  dat.GS = t(sapply(X.index,
                    function(x) {
                      t(rmultinom(1, 1, cell.prob.unique[x, ])) %*%
                      Lmat.withZero}))
  list(dat.GS = dat.GS, pars.baseline = pars,
       cell.prob.unique = cell.prob.unique, X = X)
}

SimulateGSdata2 <- function(ncase, K, Smax, Pi.seed,
                            num.covariates, Betas, theta2.value, theta2.pind,
                            phi.iter = 500, has.interact = FALSE) {
  # Simulate GS data using conditional parameterization, with Sparse Corr 1
  #
  # Args:
  #   Betas: matrix, nrow = num.covariates, ncol = K
  #
  # Example:
  #   Pi = c(0.1, 0.5, 0.3, 0.1 - 0.0003, 0.0002, 0.0001)
  #   SimulateGSdata2(10, 5, 3, Pi, 1, c(0.3, 0, 0.1, -0.1, 0.2))
  X = GenDefaultCovariates(ncase, num.covariates, has.interact)
  dmat = DesignMatrixAppxQuadExp(K, Smax)
  pars = PiToBetaTheta(K, Smax, Pi.seed, phi.iter)
  
  theta2 = rep(theta2.value, choose(K, 2)) *
           rbinom(choose(K, 2), 1, theta2.pind)
  theta1 = MuToTheta1(pars$Mu, theta2,
                      cbind(dmat$Lmat, dmat$Umat), dmat$MuMat, K)
  
  if (length(Betas) > 0) {
    Betas = matrix(Betas,
                   nrow = num.covariates + has.interact, 
                   ncol = K)
  }
  Betas = rbind(theta1, Betas)
  X.unique = uniquecombs(X)
  X.index = attr(X.unique, "index")
  cell.prob.unique = XBetaToCellProb2(K, X.unique, Betas, theta2,
                                      cbind(dmat$Lmat, dmat$Umat))
  Lmat.withZero = rbind(rep(0, K), dmat$Lmat)
  dat.GS = t(sapply(X.index,
                    function(x) {
                      t(rmultinom(1, 1, cell.prob.unique[x, ])) %*%
                        Lmat.withZero}))
  Mu.unique = cell.prob.unique[, -1] %*% t(dmat$MuMat)
  Pi.unique = cbind(cell.prob.unique[, 1],
                    cell.prob.unique[, -1] %*% t(dmat$PiMat))
  list(dat.GS = dat.GS,
       pars.baseline = list(Mu = Mu.unique,
                            Pi = Pi.unique,
                            Betas = Betas,
                            theta2 = theta2),
       cell.prob.unique = cell.prob.unique,
       X = X)
}


LtoM <- function(L, TPR, FPR){
  k = ncol(L)
  notL = 1 - L
  registerDoMC(detectCores() - 1) 
  M = foreach(i = 1:nrow(L), .combine = rbind) %dopar% {
    rbinom(k, 1, L[i, ] * TPR + notL[i, ] * FPR)
  }
  M
}


SimulatePerchData <- function(ncase, nctrl, K, Smax, Pi.seed, phi.iter = 500,
                              num.covariates, Betas, theta2.value, theta2.pind,
                              ss.tpr, bs.tpr, bs.fpr, type = "sc1",
                              has.interact = FALSE) {
#
# Example:
#   par.default = SetDefaultSimulationParameter(1)
#   do.call(SimulatePerchData, par.default)
#
  if (type == "sc1") {
    GS.obj = SimulateGSdata2(ncase, K, Smax, Pi.seed,
                             num.covariates, Betas, theta2.value, theta2.pind,
                             phi.iter, has.interact)
  } else {
    GS.obj = SimulateGSdata(ncase, K, Smax, Pi.seed,
                            num.covariates, Betas, phi.iter)
  }

  L = GS.obj$dat.GS
  MSS.case = LtoM(L, ss.tpr, 0)
  MBS.case = LtoM(L, bs.tpr, bs.fpr)
  MBS.ctrl = t(matrix(rbinom(nctrl * K, 1, bs.fpr), nrow = K))
  
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  Mu.unique = GS.obj$cell.prob.unique[, -1] %*% t(design.mat$MuMat)
   
  list(L = L, MSS.case = MSS.case, MBS.case = MBS.case, MBS.ctrl = MBS.ctrl,
       X = GS.obj$X, pars.baseline = GS.obj$pars.baseline, Smax = Smax, 
       cell.prob.unique = GS.obj$cell.prob.unique, Mu.unique = Mu.unique)
}

ReSimulateData <- function(ncase, nctrl, num.covariates, has.interact,
                           cell.prob.unique, K, Smax, 
                           ss.tpr, bs.tpr, bs.fpr, ...) {
#
  X = GenDefaultCovariates(ncase, num.covariates, has.interact)
  dmat = DesignMatrixAppxQuadExp(K, Smax)
  X.unique = uniquecombs(X)
  X.index = attr(X.unique, "index")
  Lmat.withZero = rbind(rep(0, K), dmat$Lmat)
  dat.GS = t(sapply(X.index,
                    function(x) {
                      t(rmultinom(1, 1, cell.prob.unique[x, ])) %*%
                        Lmat.withZero}))
  L = dat.GS
  MSS.case = LtoM(L, ss.tpr, 0)
  MBS.case = LtoM(L, bs.tpr, bs.fpr)
  MBS.ctrl = t(matrix(rbinom(nctrl * K, 1, bs.fpr), nrow = K))
  list(L = L, MSS.case = MSS.case, MBS.case = MBS.case, MBS.ctrl = MBS.ctrl,
       X = X, K = K, Smax = Smax, cell.prob.unique = cell.prob.unique)
}

SimulateNoRegData <- function(ncase, nctrl, theta1, theta2, 
                              ss.tpr, bs.tpr, bs.fpr) {
  # simulate data from quadratic exponential data explicitly
  K = length(theta1)
  theta.vec = c(theta1, theta2) # of length: K + choose(K, 2)
  d.mat = DesignMatrixAppxQuadExp(K, K)
  LU.mat = cbind(d.mat$Lmat, d.mat$Umat)
  potentials = c(1, exp((LU.mat %*% cbind(theta.vec))[, 1]))
  cell.probs = potentials/sum(potentials)  
  
  Lmat.withZero = rbind(rep(0, K), d.mat$Lmat)
  dat.GS = t(rmultinom(ncase, 1, cell.probs)) %*% Lmat.withZero
  L = dat.GS
  MSS.case = LtoM(L, ss.tpr, 0)
  MBS.case = LtoM(L, bs.tpr, bs.fpr)
  MBS.ctrl = t(matrix(rbinom(nctrl * K, 1, bs.fpr), nrow = K))
  list(L = L, MSS.case = MSS.case, MBS.case = MBS.case, MBS.ctrl = MBS.ctrl,
       X = matrix(1, nrow = ncase, ncol = 1), K = K,
       cell.probs = cell.probs,
       Mu = round(d.mat$MuMat %*% cbind(cell.probs[-1]), 4),
       Pr.NumPathogen = round(rbind(cell.probs[1],
                                    d.mat$PiMat %*% cbind(cell.probs[-1])), 4))
}



SetDefaultSimulationParameter <- function(option = 1) {
  set.seed(124)
  if (option == 1) {
    par.config = list(ncase = 2500, nctrl = 1500, K = 5, Smax = 3,
                      Pi.seed = c(0.05, 0.55, 0.3, 0.1, 0, 0),
                      phi.iter = 500,
                      num.covariates = 0, Betas = NULL,
                      ss.tpr = c(0.11, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 2) {
    par.config = list(ncase = 2500, nctrl = 1500, K = 5, Smax = 5,
                      Pi.seed = c(0.05, 0.5, 0.35, 0.1 - 1e-2 - 1e-3,
                                  1e-2, 1e-3),
                      phi.iter = 500,
                      num.covariates = 1,
                      Betas = c(-0.5, -0.5, 0.5, 0.5, -0.5),
                      ss.tpr = c(0.11, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 3) {
    par.config = list(ncase = 2500, nctrl = 1500, K = 10, Smax = 4,
                      Pi.seed = c(0.1, 0.4, 0.3, 0.15, 0.05, rep(0, 6)),
                      phi.iter = 500,
                      num.covariates = 0,
                      Betas = NULL,
                      ss.tpr = sample(c(0.11, 0.12, 0.08, 0.15, 0.10),
                                      10, TRUE),
                      bs.tpr = sample(c(0.8, 0.6 ,0.7 ,0.7 ,0.5), 10, TRUE),
                      bs.fpr = sample(c(0.5, 0.55, 0.40, 0.35, 0.45),
                                      10, TRUE))
  } else if (option == 4) {
    par.config = list(ncase = 2500, nctrl = 1500, K = 5, Smax = 5,
                      Pi.seed = c(0.05, 0.5, 0.35, 0.1 - 1e-2 - 1e-3,
                                  1e-2, 1e-3),
                      phi.iter = 500,
                      num.covariates = 0,
                      Betas = NULL,
                      ss.tpr = c(0.11, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 5) {
    par.config = list(ncase = 2500, nctrl = 1500, K = 3, Smax = 3,
                      Pi.seed = c(0.05, 0.7, 0.2, 0.05),
                      phi.iter = 500,
                      num.covariates = 0,
                      Betas = NULL,
                      ss.tpr = c(0.11, 0.12, 0.08),
                      bs.tpr = c(0.8, 0.6 ,0.7),
                      bs.fpr = c(0.5, 0.55, 0.40))
  } else if (option == 6) {
    par.config = list(ncase = 3000, nctrl = 500, K = 5, Smax = 5,
                      Pi.seed = c(0.05, 0.5, 0.35, 0.1 - 1e-2 - 1e-3,
                                  1e-2, 1e-3),
                      phi.iter = 500,
                      num.covariates = 2,
                      Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2,
                                0.3, 0.5, -0.2, -0.4, 0.2),
                      ss.tpr = c(0.11, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 7) {
    par.config = list(ncase = 3000, nctrl = 500, K = 5, Smax = 5,
                      Pi.seed = c(0.05, 0.5, 0.35, 0.1 - 1e-2 - 1e-3,
                                  1e-2, 1e-3),
                      phi.iter = 500,
                      num.covariates = 2,
                      has.interact = TRUE,
                      Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2,
                                0.3, 0.5, -0.2, -0.4, 0.2,
                                -0.1, 0.2, 0.1, 0.3, -0.2),
                      ss.tpr = c(0.11, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 8) {
    par.config = list(ncase = 500, nctrl = 1000, K = 5, Smax = 1,
                      Pi.seed = c(0.1, 0.9),
                      phi.iter = 500,
                      num.covariates = 2,
                      has.interact = TRUE,
                      Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2,
                                0.3, 0.5, -0.2, -0.4, 0.2,
                                -0.1, 0.2, 0.1, 0.3, -0.2),
                      ss.tpr = c(0.11, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 9) {
    par.config = list(ncase = 500, nctrl = 1000, K = 5, Smax = 5,
                      Pi.seed = c(0.05, 0.5, 0.35, 0.1 - 1e-2 - 1e-3,
                                  1e-2, 1e-3),
                      phi.iter = 500,
                      num.covariates = 2,
                      has.interact = TRUE,
                      Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2,
                                0.3, 0.5, -0.2, -0.4, 0.2,
                                -0.1, 0.2, 0.1, 0.3, -0.2),
                      ss.tpr = c(0.81, 0.82, 0.88, 0.85, 0.80),
                      bs.tpr = c(0.98, 0.96 ,0.97 ,0.97 ,0.95),
                      bs.fpr = c(0.05, 0.055, 0.040, 0.035, 0.045))
  }

  c(par.config, list(type = "sc1", theta2.value = -0.4, theta2.pind = 0.7))
}


ReorgToNplcmData <- function(sim.obj,
                             etio.names = c("A", "B", "C", "D", "E")) {
  NPPCR = rbind(sim.obj$MBS.case, sim.obj$MBS.ctrl)
  colnames(NPPCR) = etio.names
  BCX = rbind(sim.obj$MSS.case,
              matrix(NA, nrow = nrow(sim.obj$MBS.ctrl),
                     ncol = length(etio.names)))
  colnames(BCX) = etio.names
  Mobs = list(MBS = list(NPPCR = NPPCR), MSS = list(BCX = BCX))
  Y = c(rep(1, nrow(sim.obj$MSS.case)), rep(0, nrow(sim.obj$MBS.ctrl)))
  X = rbind(sim.obj$X, matrix(NA, nrow = nrow(sim.obj$MBS.ctrl),
                              ncol = ncol(sim.obj$X)))
  list(Mobs = Mobs, X = X, Y = Y)
}












