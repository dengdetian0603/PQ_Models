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

GenDefaultCovariates <- function(nsample, num.covariates) {
# Simulate covariates, assuming all X are categorical.
# Example:
#   GenDefaultCovariates(10, 1)
  x0 = rep(1, nsample)
  if (num.covariates > 0) {
    X = rbinom(nsample * num.covariates, 1, 0.5)
    X = matrix(X, nrow = nsample)
    return(cbind(x0, X))
  } else {
    return(matrix(x0, ncol = 1))
  }
}

Logit <- function(x)
{
  log(x/(1-x))
}

InvLogit <- function(x)
{
      1/(1+exp(-x))
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

PiToBetaTheta <- function(K, Smax, Pi.seed) {
# Example:
#   Pi = c(0.1, 0.5, 0.3, 0.1, 0, 0)
#   PiToBetaTheta(5, 3, Pi)
  Pi = Pi.seed[1:(1 + Smax)]
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  Phis = suppressWarnings(
      xsample(E = design.mat$PiMat, F = Pi[-1]/Pi[1],
              G = diag(rep(1, design.mat$J1)),
              H = rep(0, design.mat$J1),
              iter = 500, burnin =50,
              type = "mirror", test=FALSE)$X)
  phis = Phis[300, ]
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
    if(is.character(result)) {
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
  Beta = matrix(Beta, ncol=K)
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


SimulateGSdata <- function(ncase, K, Smax, Pi.seed, num.covariates, Betas) {
#
# Args:
#   Betas: matrix, nrow = num.covariates, ncol = K
#
# Example:
#   Pi = c(0.1, 0.5, 0.3, 0.1, 0, 0)
#   SimulateGSdata(10, 5, 3, Pi, 1, c(0, 0, 0.1, -0.1, 0.2))
  X = GenDefaultCovariates(ncase, num.covariates)
  dmat = DesignMatrixAppxQuadExp(K, Smax)
  pars = PiToBetaTheta(K, Smax, Pi.seed)
  if (length(Betas) > 0) Betas = matrix(Betas, nrow = num.covariates, ncol = K)
  Betas = rbind(pars$Beta0, Betas)
  X.unique = uniquecombs(X)
  X.index = attr(X.unique, "index")
  cell.prob.unique = XBetaToCellProb(K, X, Betas, pars$Thetas[-(1:K)],
                                     cbind(dmat$Lmat, dmat$Umat), dmat$MuMat)
  Lmat.withZero = rbind(rep(0, K), dmat$Lmat)
  dat.GS = t(sapply(X.index,
                    function(x) {
                      t(rmultinom(1, 1, cell.prob.unique[x, ])) %*%
                      Lmat.withZero}))
  list(dat.GS = dat.GS, pars.baseline = pars,
       cell.prob.unique = cell.prob.unique, X = X)
}


LtoM <- function(L, TPR, FPR){
  k = ncol(L)
  notL = 1-L
  registerDoMC(detectCores() - 1) 
  M= foreach(i=1:nrow(L), .combine=rbind) %dopar% {
    rbinom(k,1,L[i,]*TPR+notL[i,]*FPR)
  }
  M
}


SimulatePerchData <- function(ncase, nctrl, K, Smax, Pi.seed,
                              num.covariates, Betas,
                              ss.tpr, bs.tpr, bs.fpr) {
#
# Example:
#   par.default = SetDefaultSimulationParameter(1)
#   do.call(SimulatePerchData, par.default)
#
  GS.obj = SimulateGSdata(ncase, K, Smax, Pi.seed, num.covariates, Betas)
  L = GS.obj$dat.GS
  MSS.case = LtoM(L, ss.tpr, 0)
  MBS.case = LtoM(L, bs.tpr, bs.fpr)
  MBS.ctrl = t(matrix(rbinom(nctrl * K, 1, bs.fpr), nrow = K))
  list(L = L, MSS.case = MSS.case, MBS.case = MBS.case, MBS.ctrl = MBS.ctrl,
       X = GS.obj$X, pars.baseline = GS.obj$pars.baseline,
       cell.prob.unique = GS.obj$cell.prob.unique)
}


SetDefaultSimulationParameter <- function(option = 1) {
  if (option == 1) {
    par.config = list(ncase = 1000, nctrl = 1000, K = 5, Smax = 3,
                      Pi.seed = c(0.1, 0.5, 0.3, 0.1, 0, 0),
                      num.covariates = 0, Betas = NULL,
                      ss.tpr = c(0.05, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 2) {
    par.config = list(ncase = 1000, nctrl = 1000, K = 5, Smax = 3,
                      Pi.seed = c(0.1, 0.5, 0.3, 0.1, 0, 0),
                      num.covariates = 1,
                      Betas = c(0, 0, 0.1, -0.1, 0.2),
                      ss.tpr = c(0.05, 0.12, 0.08, 0.15, 0.10),
                      bs.tpr = c(0.8, 0.6 ,0.7 ,0.7 ,0.5),
                      bs.fpr = c(0.5, 0.55, 0.40, 0.35, 0.45))
  } else if (option == 3) {
    par.config = list(ncase = 1000, nctrl = 1000, K = 10, Smax = 4,
                      Pi.seed = c(0.1, 0.4, 0.3, 0.15, 0.05, rep(0, 6)),
                      num.covariates = 0,
                      Betas = NULL,
                      ss.tpr = sample(c(0.05, 0.12, 0.08, 0.15, 0.10),
                                      10, TRUE),
                      bs.tpr = sample(c(0.8, 0.6 ,0.7 ,0.7 ,0.5), 10, TRUE),
                      bs.fpr = sample(c(0.5, 0.55, 0.40, 0.35, 0.45),
                                      10, TRUE))
  }
  par.config
}















