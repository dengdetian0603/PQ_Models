setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingJAGS.R")
source("./PostProcessJAGS.R")
source("./PerchPlots.R")
source("../../BAKER/utils.R")
library(ggplot2)
# -----------------------------------------------------------------------------
top5.names = c("RSV", "RHINO", "HMPV_A_B", "PNEU", "PV_EV")

load("../WorkSpace/KenyaData_NPPCR_BCX.RData")

top5.mbs.case = data_nplcm$Mobs$MBS$NPPCR[data_nplcm$Y == 1, top5.names]
top5.mss.case = top5.mbs.case
for (i in 1:5) {
  if (top5.names[i] %in% colnames(data_nplcm$Mobs$MSS$BCX)) {
    top5.mss.case[, i] = data_nplcm$Mobs$MSS$BCX[data_nplcm$Y == 1,
                                                 top5.names[i]]
  } else {
    top5.mss.case[, i] = NA
    top5.mss.case[, i] = as.numeric(top5.mss.case[, i])
  }
}
top5.mbs.ctrl = data_nplcm$Mobs$MBS$NPPCR[data_nplcm$Y == 0, top5.names]
X.case = cbind(1, data_nplcm$X[data_nplcm$Y == 1, c("AGE", "CASECONT_SUB")])
X.case[, 3] = X.case[, 3] - 3
X.case$AGE_SEVERITY = X.case[, 2] * X.case[, 3]
colnames(X.case) = c("intercept", "AGE", "SEVERITY", "AGE*SEVERITY")
# ----------------------------------------------------------------------------
# put in real data
par.default = SetDefaultSimulationParameter(7)
sim.obj = do.call(SimulatePerchData, par.default)
sim.obj$MBS.case = top5.mbs.case
sim.obj$MSS.case = top5.mss.case
sim.obj$MBS.ctrl = top5.mbs.ctrl
sim.obj$X = X.case
# set hyper-parameter
bs.tpr.hyperpar = beta_parms_from_quantiles(c(0.5, 0.99), p = c(0.025, 0.975),
                                            plot = TRUE)
ss.tpr.hyperpar = beta_parms_from_quantiles(c(0.05, 0.2), p = c(0.025, 0.975),
                                            plot = TRUE)
hyper.pars.list = SetDefaultHyperParameters(K = 5)
hyper.pars.list$aa = ss.tpr.hyperpar$a
hyper.pars.list$bb = ss.tpr.hyperpar$b
hyper.pars.list$cc = bs.tpr.hyperpar$a
hyper.pars.list$dd = bs.tpr.hyperpar$b
hyper.pars.list$tau_theta = 0.2
hyper.pars.list$pind_a = 25
hyper.pars.list$pind_b = 15
hyper.pars.list$theta2_mu = -1.5

# specify model
model.file1 = "./jags/SparseCorr1_BSandSS_NoReg.jags"
model.file2 = "./jags/SparseCorr1_BSandSS_Reg1.jags"

# fit model
jags.result = FitJags(sim.obj,
                      c("cell_prob", "Beta", "theta2", "theta2_value",
                        "bs_tpr", "ss_tpr", "bs_fpr"),
                      model.file2, NULL, hyper.pars.list,
                      7000, 4000, 3, 1)

# collect results
coda.fit = as.mcmc(jags.result)
save(coda.fit, hyper.pars.list, sim.obj, top5.names,
     file = "../WorkSpace/Kenya_sc1.RData")

post.mean = round(summary(coda.fit)[[1]][,1], 3)
post.mean[grepl("bs_tpr", names(post.mean))]
post.mean[grepl("ss_tpr", names(post.mean))]
post.mean[grepl("theta2", names(post.mean))]
matrix(post.mean[grepl("Beta", names(post.mean))], nrow = 2)
plot(as.mcmc(coda.fit[[1]][, grepl("Beta", colnames(coda.fit[[1]]))]))

Mu.fit = ExtractMu(coda.fit[[1]], sim.obj)
plot(as.mcmc(Mu.fit[[2]]))
colMeans(Mu.fit[[1]])
colMeans(Mu.fit[[2]])

Pr.Num.Path.fit = ExtractPrNumPath(coda.fit[[1]], sim.obj)
round(colMeans(Pr.Num.Path.fit[[1]]), 3)
round(colMeans(Pr.Num.Path.fit[[2]]), 3)

ListEtiology(coda.fit[[1]], sim.obj, top5.names, FALSE)

# ----------------------------------------------------------------------------
# Compare with nplcm model
eti_log = fread("../Data/GambiaBakerPQlog.csv")
eti_top5 = eti_log[, top5.names, with = FALSE]
summary(as.mcmc(eti_top5))

PlotCompareResult(Mu.fit, eti_top5, c("SC1", "nplcm"))

# ----------------------------------------------------------------------------
Pr.NumPath.fit = ExtractPrNumPath(coda.fit[[1]], sim.obj)
Pr.NoneAbove = apply(Pr.NumPath.fit, 1, function(x) 1 - sum(x))
# Prior Informativeness evaluation
par.mat = cbind(Pr.NoneAbove, Mu.fit,
                as.data.frame(
                  coda.fit[[1]][, grepl("tpr", colnames(coda.fit[[1]]))]),
                as.data.frame(
                  coda.fit[[1]][, grepl("fpr", colnames(coda.fit[[1]]))]))
par.mat.tr = Logit(par.mat)
colnames(par.mat.tr)

par.given.ss_tpr = SweepPosterior(coda.fit[[1]], sim.obj)
pairs(par.mat.tr, pch = ".")
pairs(par.given.ss_tpr, pch = ".")

