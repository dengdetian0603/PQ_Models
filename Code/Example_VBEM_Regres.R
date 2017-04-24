setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingVB.R")
source("../../BAKER/utils.R")
Rcpp::sourceCpp('VBEM_Regression.cpp')

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
hyper.pars.list$aa = rep(25.22, 5)
hyper.pars.list$bb = rep(4.56, 5)
hyper.pars.list$cc = 38.1
hyper.pars.list$dd = 2.4
hyper.pars.list$rho_mu = -4
hyper.pars.list$rho_tau = 6

# beta_parms_from_quantiles(c(0.1, 0.12), p = c(0.025, 0.975), plot = TRUE)

# low-quality data
{
par.default = SetDefaultSimulationParameter(7)
par.default$ncase = 1000
par.default$nctrl = 1000
par.default$num.covariates = 1
par.default$has.interact = FALSE
par.default$Betas = c(-0.1, -0.3, 0.4, -0.5, 0.2)
par.default$Pi.seed = c(0.1, 0.7, 0.15, 0.05 - 1e-2 - 1e-3, 1e-2, 1e-3)

par.default$theta2.value = -1.5
par.default$theta2.pind = 0.8

par.default$bs.tpr = c(0.8, 0.6, 0.7, 0.7, 0.65)
par.default$bs.fpr = c(0.45, 0.3, 0.35, 0.4, 0.35)
#par.default$ss.tpr = c(0.5, 0.6, 0.4, 0.55, 0.45)

hyper.pars.list = SetVBHyperParameters(K = 5)
# hyper.pars.list$aa = rep(11.26, 5)  # \in (0.3, 0.7)
# hyper.pars.list$bb = rep(11.26, 5)
# hyper.pars.list$aa = rep(7.6, 5)  # \in (0.05, 0.2)
# hyper.pars.list$bb = rep(59, 5)
hyper.pars.list$aa = c(412, 485, 224, 733, 344)
hyper.pars.list$bb = c(3340, 3563, 2590, 4158, 3103)
hyper.pars.list$cc = 12.7  # \in (0.5, 0.9)
hyper.pars.list$dd = 4.8
hyper.pars.list$theta_tau = 1.5

hyper.pars.list$rho_mu = -3
hyper.pars.list$rho_tau = 26

hyper.pars.list$pind_a = 13
hyper.pars.list$pind_b = 0.6
}
# -----------------------------------------------------------------------------

set.seed(603)
sim.obj = do.call(SimulatePerchData, par.default)
round(sim.obj$pars.baseline$Mu[1, ], 4)
round(sim.obj$pars.baseline$Pi[1, ], 4)

sim.dat = do.call(ReSimulateData,
                  c(par.default,
                    list(cell.prob.unique = sim.obj$cell.prob.unique)))

group0 = which(rowSums(sim.dat$X) > 0)
input.obj = list(MSS.case = sim.dat$MSS.case[group0, ],
                 MBS.case = sim.dat$MBS.case[group0, ],
                 MBS.ctrl = sim.dat$MBS.ctrl,
                 X = sim.dat$X[group0, ],
                 ss.available = 1:5, bs.available = 1:5)

# -----------------------------------------------------------------------------
hyper.pars.list$theta_mu = matrix(0, nrow = ncol(input.obj$X),
                                  ncol = ncol(input.obj$MBS.case))

par.list = SetVBInitialValues(5, input.obj, hyper.pars.list)
t0 = proc.time()
res = regVBEM(input.obj, hyper.pars.list, par.list,
              maxIter = 100, tol = 5 * 1e-5, 2)
proc.time() - t0

res$mu_theta = uniquecombs(input.obj$X) %*% res$Beta.mean
res$mu_rho = res$Rho.mean

# -----------------------------------------------------------------------------
res$mu_theta
res$mu_rho 
round(res$mu_rho * res$qD, 3)

etio.info = VBEtioProbsReg(res)
plot(etio.info[[1]]$etio.probs.pL, type = "b", col = "blue", ylim = c(0, 0.28))
lines(sim.obj$cell.prob.unique[1, ], type = "b", col = "red")
# lines(etio.info$etio.probs, type = "b", col = "black")
plot(etio.info[[2]]$etio.probs.pL, type = "b", col = "blue", ylim = c(0, 0.28))
lines(sim.obj$cell.prob.unique[2, ], type = "b", col = "red")

etio.info$ss.tpr
etio.info$bs.tpr
etio.info$bs.fpr
# round(cbind(sim.obj$Mu.unique[1, ], etio.info$etio.mean.pL), 3)
