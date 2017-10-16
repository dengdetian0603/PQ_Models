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
hyper.pars.list$aa = rep(25.22, 5)
hyper.pars.list$bb = rep(4.56, 5)
hyper.pars.list$cc = 38.1
hyper.pars.list$dd = 2.4
hyper.pars.list$rho_mu = -1
hyper.pars.list$rho_tau = 2
hyper.pars.list$theta_mu = rep(0, 5)

# beta_parms_from_quantiles(c(0.1, 0.12), p = c(0.025, 0.975), plot = TRUE)

# low-quality data
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
hyper.pars.list$aa = rep(11.26, 5)  # \in (0.3, 0.7)
hyper.pars.list$bb = rep(11.26, 5)
# hyper.pars.list$aa = rep(7.6, 5)  # \in (0.05, 0.2)
# hyper.pars.list$bb = rep(59, 5)
# hyper.pars.list$aa = c(412, 485, 224, 733, 344)
# hyper.pars.list$bb = c(3340, 3563, 2590, 4158, 3103)
hyper.pars.list$cc = 12.7  # \in (0.5, 0.9)
hyper.pars.list$dd = 4.8

hyper.pars.list$theta_mu = c(0, 0, 0, 0, 0)
hyper.pars.list$theta_tau = 2.5

hyper.pars.list$rho_mu = -3
hyper.pars.list$rho_tau = 18

hyper.pars.list$pind_a = 9
hyper.pars.list$pind_b = 3
# -----------------------------------------------------------------------------
ss.tpr = c(0.5, 0.6, 0.4, 0.55, 0.45)

bs.tpr = c(0.8, 0.6, 0.7, 0.7, 0.65) - 0.05
bs.fpr = c(0.45, 0.3, 0.35, 0.4, 0.35) - 0.1

# theta1 = c(-1, -0.5, 0, 1, 1.5) + 0.5 # v0
# theta1 = c(-1, -0.5, 0, 1, 2) + 0 # v1
theta1 = c(-1, -0.5, 0, 1, 1.5) - 0.5 # v2
theta2 = c(-0, -1, -2, -2,
           -2, -1, -2,
           -1, -5,
           -2) * 1 # v3    # v0-1: *1.5

sim.obj = SimulateNoRegData(ncase = 250, nctrl = 1000, 
                            theta1 = theta1,
                            theta2 = theta2,
                            ss.tpr = ss.tpr,
                            bs.tpr = bs.tpr,
                            bs.fpr = bs.fpr)
sim.obj$Pr.NumPathogen
sim.obj$Mu

input.obj = c(sim.obj, list(ss.available = 1:5, bs.available = 1:5))
#------------------------------------------------------------------------------
# set.seed(603)
# sim.obj = do.call(SimulatePerchData, par.default)
# round(sim.obj$pars.baseline$Mu[1, ], 4)
# round(sim.obj$pars.baseline$Pi[1, ], 4)
# 
# sim.dat = do.call(ReSimulateData,
#                   c(par.default,
#                     list(cell.prob.unique = sim.obj$cell.prob.unique)))
# 
# group0 = which(rowSums(sim.dat$X) == 1)
# input.obj = list(MSS.case = sim.dat$MSS.case[group0, ],
#                  MBS.case = sim.dat$MBS.case[group0, ],
#                  MBS.ctrl = sim.dat$MBS.ctrl,
#                  X = cbind(sim.dat$X[group0, 1]),
#                  ss.available = 1:5, bs.available = 1:5)

# -----------------------------------------------------------------------------
hyper.pars.list = SetVBHyperParameters(K = 5)
hyper.pars.list$aa = rep(11.26, 5)  # \in (0.3, 0.7)
hyper.pars.list$bb = rep(11.26, 5)

hyper.pars.list$cc = 12.7  # \in (0.5, 0.9)
hyper.pars.list$dd = 4.8

hyper.pars.list$theta_mu = c(0, 0, 0, 0, 0) + 0.5
hyper.pars.list$theta_tau = 1

hyper.pars.list$rho_mu = -3.5
hyper.pars.list$rho_tau = 15

hyper.pars.list$pind_a = 1
hyper.pars.list$pind_b = 1

par.list = SetVBInitialValues(5, input.obj, hyper.pars.list)

par.vbfit = FitVBEMnoReg(input.obj, hyper.pars.list, par.list,
                         max.iter = 450, tol = 5 * 1e-6) 

# -----------------------------------------------------------------------------
par.vbfit$mu_theta
par.vbfit$mu_rho  
round(par.vbfit$mu_rho * par.vbfit$qD, 2)

etio.info = VBEtioProbsNoReg(par.vbfit)
plot(etio.info$etio.probs.pL, type = "b", col = "blue", ylim = c(0, 0.4))
lines(sim.obj$cell.probs, type = "b", col = "red")

# lines(baker.smax3.probs , type = "b", col = "black")
# lines(baker.smax5.probs , type = "b", col = "gray")

# - - - - -    - - --   - - - - 
sum(sqrt(etio.info$etio.probs.pL * sim.obj$cell.probs))
# sum(sqrt(sim.obj$cell.probs * baker.smax3.probs))
# sum(sqrt(sim.obj$cell.probs * baker.smax5.probs))

round(cbind(sim.obj$Mu, etio.info$etio.mean.pL), 3)
round(cbind(sim.obj$Pr.NumPathogen, etio.info$etio.numPi.pL), 3)

PredNoRegGOF(5, etio.info, 5000, input.obj, 300)[3:4]
# etio.info$ss.tpr
# etio.info$bs.tpr
# etio.info$bs.fpr
# ==============================================================================
benchmark = read.csv("../WorkSpace/BenchMarks/Baker_Smax_NoReg.csv")
baker.smax3.probs = unlist(benchmark[1, -(1:12)])
baker.smax5.probs = unlist(benchmark[2, -(1:12)])

# ==============================================================================
library(dplyr)
res = read.csv("../WorkSpace/VBexperiments/VB-NoReg-HQBS+HQSS-result-v2.csv")
res = read.csv("../WorkSpace/VBexperiments/VB-NoReg-MQBS+MQSS-result-v2.csv")
res = read.csv("../WorkSpace/Manuscript2017/VB-NoReg-MQBS+MQSS-result.csv")
###
with(res, plot(BS.gof.pred, Bhattacharyya.Coef, col = rep))
with(res, plot(SS.gof.pred, Bhattacharyya.Coef, col = rep))
with(res, plot(SS.gof.pred + BS.gof.pred,
               Bhattacharyya.Coef, col = rep))
with(res, cor(0.8 * SS.gof.pred + 0.2 * BS.gof.pred, Bhattacharyya.Coef))

#### all variation
plot(sim.obj$cell.probs, type = "n", ylim = c(0, 0.5))
for (i in 1:nrow(res)) {
  points(jitter(1:32), unlist(res[i, 3:34]), pch = 16,
         col = alpha(rgb(0, 0, 0), 0.3), cex = 0.7)
}
points(sim.obj$cell.probs, col = "red", cex = 3, pch = "-")

## for Mu
plot(sim.obj$Mu[, 1], type = "n", ylim = c(0, 0.7))
for (i in 1:nrow(res)) {
  points(jitter(1:5), unlist(res[i, 35:39]), pch = 16,
         col = alpha(rgb(0, 0, 0), 0.3), cex = 0.7)
}
points(sim.obj$Mu[, 1], col = "red", cex = 3, pch = "-")


#### variation by data in best prior
plot(sim.obj$cell.probs, type = "n", ylim = c(0, 0.5))
res.best = res[res$tid == 20, ]
for (i in 1:nrow(res.best)) {
  points(jitter(1:32), unlist(res.best[i, 3:34]), pch = 16,
         col = alpha(rgb(0, 0, 0), 0.5), cex = 1)
}
points(sim.obj$cell.probs, col = "red", cex = 3, pch = "-")


### average over data set, variation by prior
res = res %>% 
  mutate(SS.gof.pred = replace(SS.gof.pred, SS.gof.pred == -Inf, NA))  
res.by.prior = res %>% group_by(tid) %>% 
               summarise_each(funs(mean(., na.rm = TRUE))) %>%
               arrange(desc(Bhattacharyya.Coef), desc(SS.gof.pred))
bc = res.by.prior$Bhattacharyya.Coef
ss.gof = scale(res.by.prior$SS.gof.pred)
bs.gof = scale(res.by.prior$BS.gof.pred)
gof = scale(res.by.prior$SS.gof.pred + res.by.prior$BS.gof.pred)
plot(gof, bc)
points(ss.gof, bc, col = "blue")

##
plot(sim.obj$cell.probs, type = "n", ylim = c(0, 0.5))
for (i in 1:nrow(res.by.prior)) {
  points(jitter(1:32), unlist(res.by.prior[i, 3:34]), pch = 16,
         col = alpha(rgb(0, 0, 0), 0.5), cex = 1.2)
}
# gof = res.by.prior$Cell.Probs.6
lines(unlist(res.by.prior[which(gof == max(gof)), 3:34]),
      type = "b", pch = "-", col = "blue", cex = 2, lwd = 2)
points(sim.obj$cell.probs, col = "red", cex = 2, pch = "-")

## for Mu
plot(sim.obj$Mu[, 1], type = "n", ylim = c(0, 0.6))
for (i in 1:nrow(res.by.prior)) {
  points(jitter(1:5), unlist(res.by.prior[i, 35:39]), pch = 16,
         col = alpha(rgb(0, 0, 0), 0.5), cex = 1.2)
}
lines(unlist(res.by.prior[which(gof == max(gof)), 35:39]),
      type = "b", pch = "-", col = "blue", cex = 2, lwd = 2)
points(sim.obj$Mu[, 1], col = "red", cex = 2, pch = "-")




#### averaged over different prior, variation by data set
res.by.dat = res %>% group_by(rep) %>% 
             summarise_each(funs(mean(., na.rm = TRUE))) %>%
             arrange(desc(Bhattacharyya.Coef), desc(SS.gof.pred))
plot(sim.obj$cell.probs, type = "n", ylim = c(0, 0.5))
for (i in 1:nrow(res.by.dat)) {
  points(jitter(1:32), unlist(res.by.dat[i, 3:34]), pch = 16,
         col = alpha(rgb(0, 0, 0), 0.5), cex = 1.2)
}
points(sim.obj$cell.probs, col = "red", cex = 2, pch = "-")
