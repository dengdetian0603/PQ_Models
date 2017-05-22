library(Rcpp)
sourceCpp("./VBEM_Update.cpp")
# -----------------------------------------------------------------------------
SetVBHyperParameters <- function(K = 5) {
  list(aa = rep(2, K), bb = rep(18, K), # ss.tpr
       cc = 5, dd = 3,  # bs.tpr
       ee = 1, ff = 1,  # bs.fpr
       theta_mu = rbind(rep(-0.5, K)), # theta(1)
       rho_mu = -1.5, # theta(2)
       theta_tau = 0.2,
       rho_tau = 0.2,
       pind_a = 6, pind_b = 2 # pi_d
       )
}

SetVBInitialValues <- function(K = 5, input.obj, hyper.par.list) {
  # take input arguments
  ss.avail = input.obj$ss.available
  bs.avail = input.obj$bs.available
  n.case = nrow(input.obj$MBS.case)
  # initialize wth prior
  init.vals = list(A_st = rep(NA, K), B_st = rep(NA, K),
                   A_bt = rep(NA, K), B_bt = rep(NA, K),
                   A_bf = rep(NA, K), B_bf = rep(NA, K),
                   qL = matrix(1/K, nrow = n.case, ncol = K),
                   mu_theta = hyper.par.list$theta_mu,
                   tau_theta = matrix(hyper.par.list$theta_tau,
                                      nrow = ncol(input.obj$X), ncol = K),
                   mu_rho = hyper.par.list$rho_mu,
                   tau_rho = hyper.par.list$theta_tau,
                   qD = matrix(0.5, nrow = K, ncol = K)
                  )
  init.vals$A_st[ss.avail] = hyper.par.list$aa[ss.avail]
  init.vals$B_st[ss.avail] = hyper.par.list$bb[ss.avail]
  init.vals$A_bt[bs.avail] = rep(hyper.par.list$cc, length(bs.avail))
  init.vals$B_bt[bs.avail] = rep(hyper.par.list$dd, length(bs.avail))
  init.vals$A_bf[bs.avail] = colSums(input.obj$MBS.ctrl[, bs.avail],
                                     na.rm = TRUE)
  init.vals$B_bf[bs.avail] = colSums(1 - input.obj$MBS.ctrl[, bs.avail],
                                     na.rm = TRUE)
  diag(init.vals$qD) = 0
  init.vals
}

# -----------------------------------------------------------------------------
LoopyUpdateQL = function(qL, qD, H, mu.rho, nloop = 4) {
  qLD = qL %*% qD
  for (loop in 1:nloop) {
    for (k in 1:ncol(qL)) {
      qL.k = qL[, k]
      qL[, k] = 1/(1 + exp(-(H[, k] + qLD[, k] * mu.rho)))
      d.qL.k =  qL[, k] - qL.k
      qLD = qLD + d.qL.k %o% qD[k, ]
    }
  }
  qL
}


UpdateTheta = function(theta.init, qL, qD, qLD, mu.rho, tau.rho,
                       theta.mu, theta.tau) {
  # UpdateTheta(rep(0, 5), qL, qD, qLD, -2, 0.2, rep(-2, 5), 0.2)
  #
  qL2D2 = (qL ^ 2) %*% (qD ^ 2)
  var.rhoR = (mu.rho ^ 2 + 1/tau.rho) * (qLD - qL2D2) + 1/tau.rho * qLD ^ 2
  rho.qLD = mu.rho * qLD
  K = ncol(qL)
  result = list(mu.theta = rep(NA, K), tau.theta = rep(NA ,K))
  for (k in 1:K) {
    QThetaK = function(x) {
      e.theta = exp(x + rho.qLD[, k])
      fmat = x * qL[, k] - log(1 + e.theta) -
             (e.theta * var.rhoR[, k])/(2 * (1 + e.theta) ^ 2)
      fval = sum(fmat, na.rm = TRUE) - 0.5 * theta.tau * (x - theta.mu[k]) ^ 2
      fval
    }
    thetak = optim(par = theta.init[k], fn = QThetaK, method = "BFGS",
                   control = list(fnscale = -1), hessian = TRUE)
    result$mu.theta[k] = thetak$par
    result$tau.theta[k] = abs(thetak$hessian)
    if (thetak$convergence > 0) {
      print(thetak$message)
    }
  }
  result
}


UpdateRho = function(rho.init, qL, qLD, mu.theta, tau.theta,
                     rho.mu, rho.tau) {
  # UpdateRho(-2, qL, qLD, rep(2, 5), rep(0.5, 5), -2, 0.2)
  #
  K = ncol(qL)
  fmat1 = qL * qLD
  QRho = function(x) {
    fmat2 = matrix(0, nrow = K, ncol = nrow(qL))
    for (j in 0:(K - 1)) {
      e.rho = exp(mu.theta + j * x)
      f.rho = log(1 + e.rho) + e.rho/(2 * tau.theta * (1 + e.rho) ^ 2)
      weight.qLD = dpois(j, qLD)
      fmat2 = fmat2 + t(weight.qLD) * f.rho
    }
    fval = sum(x * fmat1 - t(fmat2), na.rm = TRUE) - 0.5 * rho.tau * (x - rho.mu) ^ 2
    # fval = sum(x * fmat1 - t(fmat2)) + rho.tau * log(-x) - rho.mu * x
    fval
  }
  rho.obj = optim(par = rho.init, fn = QRho, method = "BFGS",
                  control = list(fnscale = -1), hessian = TRUE)
  # rho.init = -abs(rho.init) - 1e-3
  # print(rho.init)
  # rho.obj = constrOptim(theta = rho.init, f = QRho, ui = matrix(-1), ci = 1e-3,
  #                       control = list(fnscale = -1),
  #                       method = NULL, hessian = TRUE)
  if (rho.obj$convergence > 0) {
    print(rho.obj$message)
  }
  list(mu.rho = rho.obj$par, tau.rho = abs(rho.obj$hessian)[1, 1])
}


LoopyUpdateQD = function(par.list, qLD, pind_a, pind_b) {
  # loopyUpdateQD(qL, qD, qLD, qLL, qL2D2, rep(-1, 5), rep(0.4, 5), -2, 0.4, 1.2, 3)
  # 
  psi.A.B = digamma(pind_a) - digamma(pind_b)
  qLL = crossprod(par.list$qL, par.list$qL)
  qL2D2 = (par.list$qL ^ 2) %*% (par.list$qD ^ 2)
  qD = loopyUpdateQD(par.list$qL, par.list$qD, qLD, qLL, qL2D2,
                     par.list$mu_theta, par.list$tau_theta,
                     par.list$mu_rho, par.list$tau_rho, psi.A.B, 
                     nLoop = 3)
  qD
}

# -----------------------------------------------------------------------------
FitVBEMnoReg = function(input.obj, hyper.par.list, init.val,
                        max.iter = 500, tol = 1e-6) {
  # 
  # prepare data
  ss.avail = input.obj$ss.available
  bs.avail = input.obj$bs.available
  MSS.case = as.matrix(input.obj$MSS.case)
  MBS.case = as.matrix(input.obj$MBS.case)
  MBS.ctrl = as.matrix(input.obj$MBS.ctrl)
  # 
  par.list$mu_theta = as.vector(par.list$mu_theta)
  par.list$tau_theta  = as.vector(par.list$tau_theta)
  hyper.pars.list$theta_mu = as.vector(hyper.pars.list$theta_mu)
  #
  n.case = nrow(MBS.case)
  n.ctrl = nrow(MBS.ctrl)
  par.list = init.val
  par.diff = 1
  par.curr = unlist(par.list[-7])
  iter = 0
  while (abs(par.diff) > tol & iter < max.iter) {
    par.prev = par.curr
    # update qL
    H = matrix(par.list$mu_theta, nrow = n.case,
               ncol = ncol(MBS.case), byrow = TRUE)
    psi.A.A = (digamma(par.list$A_bt) - digamma(par.list$A_bf))[bs.avail]
    psi.B.B = (digamma(par.list$B_bt) - digamma(par.list$B_bf))[bs.avail]
    psi.A.B = (digamma(par.list$A_bf + par.list$B_bf) -
               digamma(par.list$A_bt + par.list$B_bt))[bs.avail]
    H.bs = t(t(MBS.case[, bs.avail]) * psi.A.A +
             t(1 - MBS.case[, bs.avail]) * psi.B.B + psi.A.B)
    H[, bs.avail] = H[, bs.avail] + H.bs
    if (length(ss.avail) > 0) {
      H[, ss.avail] = t(t(H[, ss.avail]) +
                        digamma(par.list$B_st[ss.avail]) -
                        digamma(par.list$A_st[ss.avail] +
                                par.list$B_st[ss.avail])) +
                      MSS.case[, ss.avail] * 1e+16
    }
    qL = LoopyUpdateQL(par.list$qL, par.list$qD, H,
                       par.list$mu_rho, nloop = 4)
    par.list$qL = qL

    # update A, B
    q.sum = colSums(par.list$qL, TRUE)
    ## silver
    qss.sum = q.sum[ss.avail]
    if (length(ss.avail) > 0) {
      qMs.sum = colSums(par.list$qL[, ss.avail] * MSS.case[, ss.avail], TRUE)
      par.list$A_st[ss.avail] = qMs.sum + hyper.par.list$aa[ss.avail]
      par.list$B_st[ss.avail] = qss.sum - qMs.sum + hyper.par.list$bb[ss.avail]
    }
    ## bronze
    qbs.sum = q.sum[bs.avail]
    qMb.sum = colSums(par.list$qL[, bs.avail] * MBS.case[, bs.avail], TRUE)
    par.list$A_bt[bs.avail] = qMb.sum + hyper.par.list$cc
    par.list$B_bt[bs.avail] = qbs.sum - qMb.sum + hyper.par.list$dd
    Msum.case = colSums(MBS.case[, bs.avail], TRUE)
    Msum.ctrl = colSums(MBS.ctrl[, bs.avail], TRUE)
    par.list$A_bf[bs.avail] = Msum.case + Msum.ctrl -
                                qMb.sum + hyper.par.list$ee
    par.list$B_bf[bs.avail] = n.case + n.ctrl + qMb.sum - qbs.sum -
                                Msum.case - Msum.ctrl + hyper.par.list$ff

    # update mu_theta, tau_theta  
    qLD = par.list$qL %*% par.list$qD
    theta.update = UpdateTheta(par.list$mu_theta, par.list$qL, par.list$qD, qLD,
                               par.list$mu_rho, par.list$tau_rho,
                               hyper.par.list$theta_mu,
                               hyper.par.list$theta_tau)
    par.list$mu_theta = theta.update$mu.theta
    par.list$tau_theta = theta.update$tau.theta

    # update mu_rho, tau_rho
    rho.update = UpdateRho(par.list$mu_rho, par.list$qL, qLD,
                           par.list$mu_theta, par.list$tau_theta,
                           hyper.par.list$rho_mu, hyper.par.list$rho_tau)
    par.list$mu_rho = rho.update$mu.rho
    par.list$tau_rho = rho.update$tau.rho
    
    # update qD
    qD = LoopyUpdateQD(par.list, qLD,
                       hyper.par.list$pind_a, hyper.par.list$pind_b)
    par.list$qD = qD
    
    # stop signal
    par.curr = unlist(par.list[-7])
    par.diff = sqrt(sum((par.curr - par.prev) ^ 2, na.rm = TRUE))
    iter = iter + 1
    if (iter %% 50 == 0) {
      print(paste(iter, ":", par.diff))
    }
  }
  return(par.list)
}

# -----------------------------------------------------------------------------
VBEtioProbsNoReg = function(par.list, Smax = NULL) {
  K = ncol(par.list$qL)
  if (length(Smax) < 1) {
    Smax = K
  }
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  Lall = rbind(rep(0, K), design.mat$Lmat)
  potentials = exp(Lall %*% par.list$mu_theta +
                   rowSums(Lall %*% (par.list$mu_rho * par.list$qD) * Lall))
  Aconst = sum(potentials)
  etio.probs = potentials/Aconst
  etio.mean = design.mat$MuMat %*% cbind(etio.probs[-1])
  etio.numPi = c(etio.probs[1],
                 (design.mat$PiMat %*% cbind(etio.probs[-1]))[, 1]) 
  
  Apseudo = apply(1 + exp(par.list$mu_theta +
                          t(Lall %*% (par.list$mu_rho *
                                      par.list$qD))), 2, prod)
  etio.probs.pseudo = potentials/Apseudo
  etio.probs.pseudo = etio.probs.pseudo/sum(etio.probs.pseudo)
  etio.mean.pseudo = design.mat$MuMat %*% cbind(etio.probs.pseudo[-1])
  etio.numPi.pseudo = c(etio.probs.pseudo[1],
                        (design.mat$PiMat %*%
                         cbind(etio.probs.pseudo[-1]))[, 1]) 
  # round(cbind(etio.probs, etio.probs.pseudo, sim.obj$cell.prob.unique[1, ]), 3)
  list(etio.probs = etio.probs, etio.mean = etio.mean, etio.numPi = etio.numPi,
       etio.probs.pL = etio.probs.pseudo, etio.mean.pL = etio.mean.pseudo,
       etio.numPi.pL = etio.numPi.pseudo,
       ss.tpr = par.list$A_st/(par.list$A_st + par.list$B_st),
       bs.tpr = par.list$A_bt/(par.list$A_bt + par.list$B_bt),
       bs.fpr = par.list$A_bf/(par.list$A_bf + par.list$B_bf))
}

VBEtioProbsReg = function(par.list, Smax = NULL) {
  K = ncol(par.list$qL)
  if (length(Smax) < 1) {
    Smax = K
  }
  design.mat = DesignMatrixAppxQuadExp(K, Smax)
  Lall = rbind(rep(0, K), design.mat$Lmat)
  n.strata = nrow(par.list$mu_theta)
  result = list()
  for (i in 1:n.strata) {
    potentials = exp(Lall %*% par.list$mu_theta[i, ] +
                       rowSums(Lall %*% (par.list$mu_rho * par.list$qD) * Lall))
    Apseudo = apply(1 + exp(par.list$mu_theta[i, ] +
                              t(Lall %*% (par.list$mu_rho *
                                            par.list$qD))), 2, prod)
    etio.probs.pseudo = potentials/Apseudo
    etio.probs.pseudo = etio.probs.pseudo/sum(etio.probs.pseudo)
    etio.mean.pseudo = design.mat$MuMat %*% cbind(etio.probs.pseudo[-1])
    etio.numPi.pseudo = c(etio.probs.pseudo[1],
                          (design.mat$PiMat %*%
                             cbind(etio.probs.pseudo[-1]))[, 1]) 
    result[[i]] = list(etio.probs.pL = etio.probs.pseudo,
                       etio.mean.pL = etio.mean.pseudo,
                       etio.numPi.pL = etio.numPi.pseudo)
  }
  result = c(result,
             list(ss.tpr = par.list$A_st/(par.list$A_st + par.list$B_st),
                  bs.tpr = par.list$A_bt/(par.list$A_bt + par.list$B_bt),
                  bs.fpr = par.list$A_bf/(par.list$A_bf + par.list$B_bf)))
  result  
}


ResampleCaseCtrl = function(input.obj) {
  nCase = nrow(input.obj$MBS.case)
  nCtrl = nrow(input.obj$MBS.ctrl)
  boot.case = sample(1:nCase, nCase, replace = TRUE)
  boot.ctrl = sample(1:nCtrl, nCtrl, replace = TRUE)
  input.obj$MSS.case = input.obj$MSS.case[boot.case, ]
  input.obj$MBS.case = input.obj$MBS.case[boot.case, ]
  input.obj$X = input.obj$X[boot.case, ]
  input.obj$MBS.ctrl = input.obj$MBS.ctrl[boot.ctrl, ]
  input.obj
}
