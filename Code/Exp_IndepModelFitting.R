source("./PerchDataSimulation.R")
source("./IndepModel.R")

par.default = SetDefaultSimulationParameter(1)
par.default$ncase = 1500

set.seed(123)
sim.obj = do.call(SimulatePerchData, par.default)
direct.fit = FitMLEgivenGSNoReg(par.default$K, sim.obj)

mc.fit1.0 = FitJagsNoReg(sim.obj, c("mu"),
                         "./jags/Indep_BSonly_FixPR_NoReg.txt",
                         par.default$bs.tpr, n.chains = 1)

mc.fit1.1 = FitJagsNoReg(sim.obj, c("mu"),
                      "./jags/Indep_BSonly_FixPR_NoReg.txt",
                      direct.fit$bs.tpr.mle, n.chains = 1)

mc.fit1.2 = FitJagsNoReg(sim.obj, c("mu"),
                         "./jags/Indep_BSonly_FixPR_NoReg.txt",
                         direct.fit$bs.tpr.hat, n.chains = 1)

mc.fit2 = FitJagsNoReg(sim.obj, c("mu", "bs_tpr"),
                       "./jags/Indep_BSandSSpos_NoReg.txt",
                       direct.fit$bs.tpr.mle, n.chains = 1)

mc.fit3 = FitJagsNoReg(sim.obj, c("mu", "bs_tpr"),
                       "./jags/Indep_BSandSS_NoReg.txt",
                       direct.fit$bs.tpr.mle, n.chains = 1)

