library('R2jags')
library('rjags')

source('./PerchDataSimulation.R')

par.default = SetDefaultSimulationParameter(1)
par.default$ncase = 4000
par.default$nctrl = 4000
sim.obj = do.call(SimulatePerchData, par.default)
# -----------------------------------------------------------------------------

sim.dat <- list(# data
	            N_ctrl = nrow(sim.obj$MBS.ctrl),
	            K = ncol(sim.obj$MBS.ctrl),
                mbs_ctrl = sim.obj$MBS.ctrl,
                N_case = nrow(sim.obj$MBS.case),
                mbs_case = sim.obj$MBS.case,
                mss_case = sim.obj$MSS.case,
                # hyper parameters
                ma = 3, mb = 7, 
                aa = 1, bb = 9,
                cc = 4, dd = 2,
                ee = 1, ff = 1)

# Fit jags model
# use R2jags lib
bayes.mod.params <- c("mu", "bs_tpr")
bayes.mod.fit <- jags.parallel(
	                 data = sim.dat,
                     #inits = bayes.mod.inits,
                     parameters.to.save = bayes.mod.params,
                     n.chains = 3,
                     n.iter = 8000,
                     n.burnin = 4000,
                     n.thin = 2,
                     model.file = './jags/Indep_BSandSSpos_NoReg.txt')
print(bayes.mod.fit)

#bayes.mod.fit.upd <- autojags(bayes.mod.fit)
#print(bayes.mod.fit.upd)

sim.obj$pars.baseline$Mu
#par.default$ss.tpr
par.default$bs.tpr