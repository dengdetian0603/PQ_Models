library('R2jags')
library('rjags')

# -----------------------------------------------------------------------------
FitMLEgivenGSNoReg <- function(K, sim.obj) {
  bs.tpr.mle = rep(NA, K)
  bs.tpr.hat = rep(NA, K)
  for (j in 1:K) {
    bs.tpr.hat[j] =
        sum(sim.obj$MBS.case[, j] * sim.obj$MSS.case[, j]) /
        sum(sim.obj$MSS.case[, j])
    bs.tpr.mle[j] =
        sum(sim.obj$MBS.case[, j] * sim.obj$L[, j]) /
        sum(sim.obj$L[, j])
  }  
  mu.mle = rep(NA, K)
  mu.mle = apply(sim.obj$L, 2, mean)
  data.frame(mu.mle, bs.tpr.mle, bs.tpr.hat)
}

SetDefaultHyperParameters <- function() {
  list(ma = 3, mb = 7, 
       aa = 1, bb = 9,
       cc = 4, dd = 2,
       ee = 1, ff = 1)
}

FitJagsNoReg <- function(sim.obj, parameters.to.save,
                         model.file, bs.tpr.prefix,
                         hyper.pars.list = SetDefaultHyperParameters(),
                         n.iter = 6000, n.burnin = 4000,
                         n.thin = 2, n.chains = 2) {
  sim.dat <- c(list(
      # data
      N_ctrl = nrow(sim.obj$MBS.ctrl),
      K = ncol(sim.obj$MBS.ctrl),
      mbs_ctrl = sim.obj$MBS.ctrl,
      N_case = nrow(sim.obj$MBS.case),
      mbs_case = sim.obj$MBS.case,
      mss_case = sim.obj$MSS.case,
      # hyper parameters
      bs_tpr_prefix = bs.tpr.prefix), hyper.pars.list)
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
