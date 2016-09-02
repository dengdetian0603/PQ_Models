library('R2jags')
library('rjags')

source('./PerchDataSimulation.R')

par.default = SetDefaultSimulationParameter(1)
par.default$ncase = 1500

set.seed(123)
sim.obj = do.call(SimulatePerchData, par.default)
# -----------------------------------------------------------------------------
sim.dat <- list(# data
	            N_ctrl = nrow(sim.obj$MBS.ctrl),
	            K = ncol(sim.obj$MBS.ctrl),
                mbs_ctrl = sim.obj$MBS.ctrl,
                N_case = nrow(sim.obj$MBS.case),
                mbs_case = sim.obj$MBS.case,
                # hyper parameters
                bs_tpr = par.default$bs.tpr,
                ma = 3, mb = 7, 
                ee = 1, ff = 1)

for (j in 1:K) {
	print(
	    sum(sim.obj$MBS.case[, j] * sim.obj$MSS.case[, j]) /
	    sum(sim.obj$MSS.case[, j])
	)
}

bayes.mod.params <- c("mu", "bs_fpr")
bayes.mod.fit <- jags.parallel(
	                 data = sim.dat,
                     #inits = bayes.mod.inits,
                     parameters.to.save = bayes.mod.params,
                     n.chains = 3,
                     n.iter = 2000,
                     n.burnin = 400,
                     n.thin = 2,
                     model.file = './jags/Indep_BSonly_FixPR_NoReg.txt')
print(bayes.mod.fit)


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
                     n.chains = 2,
                     n.iter = 6000,
                     n.burnin = 3000,
                     n.thin = 2,
                     model.file = './jags/Indep_BSandSSpos_NoReg.txt')
print(bayes.mod.fit)

#bayes.mod.fit.upd <- autojags(bayes.mod.fit)
#print(bayes.mod.fit.upd)

sim.obj$pars.baseline$Mu
#par.default$ss.tpr
par.default$bs.tpr

# -----------------------------------------------------------------------------
par.default = SetDefaultSimulationParameter(1)
par.default$ncase = 1500
par.default$Smax = 2
set.seed(123)
sim.obj1 = do.call(SimulatePerchData, par.default)
sim.obj2 = do.call(SimulatePerchData, par.default)