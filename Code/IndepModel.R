library('coda')
library('R2jags')
library('rjags')

# -----------------------------------------------------------------------------
FitMLEgivenGSNoReg <- function(K, sim.obj) {
# Get bs.tpr and mu estimtates assuming we know GS data
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
# Fit JAGS non-regression models
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


SimStudyNoReg <- function(sim.pars = SetDefaultSimulationParameter(1),
                          par.to.save = c("mu"),
                          model.file,
                          bs.tpr.option = 0,
                          n.iter, n.burnin, n.thin,
                          n.rep, max.core = 40) {
# Simulation study for JAGS non-regression models
# Example:
#   SimStudyNoReg(model.file = "./jags/Indep_BSandSSpos_NoReg.txt",
#                 n.iter = 6000, n.burnin = 3000, n.thin = 3, n.rep = 3)
#
  registerDoMC(min(detectCores() - 1, max.core))
  mc.fit.all = foreach(i = 1:n.rep, .combine = rbind) %dopar% {
    print(paste0(i, "th rep..."))
    set.seed(i*123)
    sim.obj = do.call(SimulatePerchData, sim.pars)
    direct.fit = FitMLEgivenGSNoReg(sim.pars$K, sim.obj)
    if (bs.tpr.option == 0) {
      bs.tpr.prefix = sim.pars$bs.tpr
    } else {
      bs.tpr.prefix = direct.fit$bs.tpr.hat
    }
    mc.fit = FitJagsNoReg(sim.obj, par.to.save,
                          model.file,
                          bs.tpr.prefix,
                          n.iter = n.iter, n.burnin = n.burnin,
                          n.thin = n.thin, n.chains = 1)
    coda.fit = as.mcmc(mc.fit)
    post.mean = summary(coda.fit)[[1]][,1]
    c(post.mean, direct.fit$bs.tpr.hat)
  }
  mc.fit.all
}


