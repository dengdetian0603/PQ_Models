setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingJAGS.R")
source("./PostProcessJAGS.R")
source("../../BAKER/utils.R")
library(ggplot2)
# -----------------------------------------------------------------------------
top5.names = c("PARA_1", "PARA_3", "RSV", "PNEU", "HMPV_A_B")

gambia = read.csv("../Data/Gambia.csv")

top5.data = GetTop5(gambia, top5.names)

# ----------------------------------------------------------------------------
# put in real data
par.default = SetDefaultSimulationParameter(4)
sim.obj = do.call(SimulatePerchData, par.default)
sim.obj$MBS.case = top5.data$top5.mbs.case
sim.obj$MSS.case = top5.data$top5.mss.case
sim.obj$MBS.ctrl = top5.data$top5.mbs.ctrl

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
hyper.pars.list$pind_a = 3
hyper.pars.list$pind_b = 2
hyper.pars.list$theta2_mu = -1.5

# specify model
model.file1 = "./jags/SparseCorr1_BSandSS_NoReg.jags"

# fit model
jags.result = FitJagsNoReg(sim.obj,
                           c("cell_prob", "theta1", "theta2",
                             "bs_tpr", "ss_tpr", "bs_fpr"),
                           model.file1, NULL, hyper.pars.list,
                           6000, 2000, 2, 1)

# collect results
coda.fit = as.mcmc(jags.result)
post.mean = round(summary(coda.fit)[[1]][,1], 3)
post.mean[grepl("theta2", names(post.mean))]
post.mean[grepl("theta1", names(post.mean))]
post.mean[grepl("tpr", names(post.mean))]
plot(coda.fit)

Mu.fit = ExtractMu(coda.fit[[1]], sim.obj)
colnames(Mu.fit) = top5.names
plot(as.mcmc(Mu.fit))
summary(as.mcmc(Mu.fit))

ListEtiology(coda.fit[[1]], sim.obj, top5.names)

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
# -----------------------------------------------------------------------------
# Calculate Covariance Matirces for Posterior Distn
## None regression version
var.par = covRob(data = par.mat.tr)$cov

var.par.reg = soft.thresholding(var.par, 0)
diag(var.par.reg) = diag(var.par)
v.post  = var.par.reg

par.names = colnames(v.post)

# Calculate Covariance Matirces for Prior Distn
## None regression version
v.prior = VarLogitPriorSC1(5, 5, hyper.pars.list)
colnames(v.prior) = colnames(v.post)

# Calculate Prior Infromativeness
p.i = PriorInfo(v.post, v.prior, match.prior = TRUE)
InfoFromLikelihood = 1 - p.i
write.csv(InfoFromLikelihood, "../Reports/2016_11_17/priorinfo.csv")
