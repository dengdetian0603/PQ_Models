setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
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

sim.study.obj = SimStudyNoRegJAGS(par.to.save = c("mu", "bs_tpr"),
                                  model.file =  "./jags/Indep_BSandSSpos_NoReg.txt",
                                  bs.tpr.option = 1,
                                  n.iter = 60, n.burnin = 30,
                                  n.thin = 1, n.rep = 2)

#-----------------------------------------------------------------------------
model.file = "./stan/Indep_BSandSSpos_NoReg.stan"
model.file = "./stan/Indep_BSandSSpos_NoReg_v2.stan"
model.file = "./stan/Indep_LogitPrior_BSandSSpos_NoReg_v2.stan"
par.default = SetDefaultSimulationParameter(4)
stan.model1 = stan_model(model.file)
stan.model2 = stan_model(model.file)
stan.model3 = stan_model(model.file)

par.default$ncase = 5000
par.default$nctrl = 2000
set.seed(999)
sim.obj = do.call(SimulatePerchData, par.default)
sim.obj$pars.baseline$Pi
direct.fit = FitMLEgivenGSNoReg(par.default$K, sim.obj)
stan.result = FitStanNoReg(sim.obj, "mu", model.file, method = c("HMC"))
matrix(stan.result$MAP$par[-(1:5)], ncol = 4)
round(cbind(sim.obj$pars.baseline$Mu,
            par.default$ss.tpr,
            par.default$bs.tpr,
            par.default$bs.fpr), 3)

sim.study.2 = SimStudyNoRegSTAN(
    sim.obj0 = sim.obj,
    sim.pars = SetDefaultSimulationParameter(4),
    model.file = model.file,
    par.to.save = c("mu", "bs_tpr", "ss_tpr"), method = "VB",
    n.rep = 50, seed0 = 21209)

true.par = TrueParVal(sim.study.2, par.to.save)
PlotSimStudy(sim.study.2$stan.fit.all, true.par)