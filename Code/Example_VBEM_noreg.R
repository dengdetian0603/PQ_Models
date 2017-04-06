setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingVB.R")
source("../../BAKER/utils.R")

# high-quality data
par.default = SetDefaultSimulationParameter(9)
par.default$ncase = 1000
par.default$nctrl = 1000
par.default$num.covariates = 1
par.default$has.interact = FALSE
par.default$Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2)
par.default$theta2.value = -1.5
par.default$theta2.pind = 0.8
par.default$Pi.seed = c(0.1, 0.7, 0.15, 0.05 - 1e-2 - 1e-3,
                        1e-2, 1e-3)

hyper.pars.list = SetVBHyperParameters(K = 5)
hyper.pars.list$aa = rep(31, 5)
hyper.pars.list$bb = rep(8.26, 5)
hyper.pars.list$cc = 38
hyper.pars.list$dd = 2.4
hyper.pars.list$rho_mu = -2
hyper.pars.list$theta_mu = rep(0, 5)

# beta_parms_from_quantiles(c(0.35, 0.65), p = c(0.025, 0.975), plot = TRUE)

# low-quality data
par.default = SetDefaultSimulationParameter(7)
par.default$ncase = 1000
par.default$nctrl = 1000
par.default$num.covariates = 1
par.default$has.interact = FALSE
par.default$Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2)
par.default$ss.tpr = c(0.5, 0.6, 0.4, 0.55, 0.45)
par.default$bs.tpr[5] = 0.65
par.default$bs.fpr = c(0.3, 0.2, 0.35, 0.25, 0.15)
par.default$theta2.value = -1.5
par.default$theta2.pind = 0.8
par.default$Pi.seed = c(0.1, 0.7, 0.15, 0.05 - 1e-2 - 1e-3,
                        1e-2, 1e-3)

hyper.pars.list = SetVBHyperParameters(K = 5)
hyper.pars.list$aa = rep(20.6, 5)
hyper.pars.list$bb = rep(20.6, 5)
hyper.pars.list$cc = 55
hyper.pars.list$dd = 23
hyper.pars.list$rho_mu = -4
hyper.pars.list$rho_tau = 15
hyper.pars.list$theta_mu = rep(0, 5)
hyper.pars.list$theta_tau = 0.5
# -----------------------------------------------------------------------------

sim.obj = do.call(SimulatePerchData, par.default)
round(sim.obj$pars.baseline$Mu[1, ], 4)
round(sim.obj$pars.baseline$Pi[1, ], 4)

sim.obj = do.call(ReSimulateData,
                  c(par.default,
                    list(cell.prob.unique = sim.obj$cell.prob.unique)))

group0 = which(rowSums(sim.obj$X) == 1)
input.obj = list(MSS.case = sim.obj$MSS.case[group0, ],
                 MBS.case = sim.obj$MBS.case[group0, ],
                 MBS.ctrl = sim.obj$MBS.ctrl,
                 ss.available = c(), bs.available = 1:5)

# -----------------------------------------------------------------------------
par.list = SetVBInitialValues(5, input.obj, hyper.pars.list)

par.vbfit = FitVBEMnoReg(input.obj, hyper.pars.list, par.list,
                         max.iter = 450, tol = 1e-6) 

# -----------------------------------------------------------------------------
par.vbfit$mu_theta
par.vbfit$mu_rho * par.vbfit$qD

etio.info = VBEtioProbsNoReg(par.vbfit)
plot(etio.info$etio.probs.pL, type = "b", col = "blue")
lines(sim.obj$cell.prob.unique[1, ], type = "b", col = "red")
# lines(etio.info$etio.probs, type = "b", col = "black")

etio.info$ss.tpr
etio.info$bs.tpr
etio.info$bs.fpr
# round(cbind(sim.obj$Mu.unique[1, ], etio.info$etio.mean.pL), 3)