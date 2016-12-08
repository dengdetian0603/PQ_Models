setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingJAGS.R")
source("./PostProcessJAGS.R")
source("./PerchPlots.R")

par.default = SetDefaultSimulationParameter(7)
par.default$ncase = 3000
par.default$nctrl = 1000

hyper.pars.list = SetDefaultHyperParameters(K = 5)
hyper.pars.list$aa = 7.6
hyper.pars.list$bb = 59
hyper.pars.list$pind_a = 3
hyper.pars.list$pind_b = 2
hyper.pars.list$theta2_mu = -3
hyper.pars.list$tau_theta = 0.16

# ---------------------------------------------------------
tuning.grid = expand.grid(mu.theta1 = c(-1, -0.75, -0.5),
                          mu.theta2 = c(-0.5, -2, -4),
                          pind.a = c(5, 10, 15))
tuning.grid$pind.b = 20 - tuning.grid$pind.a

gp = list()
for (i in 1:nrow(tuning.grid)) {
  hyper.pars.list$tau_theta = 0.2
  hyper.pars.list$mu_logit = rep(tuning.grid$mu.theta1[i], 5)
  hyper.pars.list$theta2_mu = tuning.grid$mu.theta2[i]
  hyper.pars.list$pind_a = tuning.grid$pind.a[i]
  hyper.pars.list$pind_b = tuning.grid$pind.b[i]
  gp[[i]] = PlotPriorPrNum(5, 5, hyper.pars.list, FALSE)
}
do.call(grid.arrange, gp[1:9])
do.call(grid.arrange, gp[10:18])
do.call(grid.arrange, gp[19:27])
# ---------------------------------------------------------

sim.obj = do.call(SimulatePerchData, par.default)
sim.obj$pars.baseline$Mu
sim.obj$pars.baseline$Pi

colSums(sim.obj$MSS.case)
colMeans(sim.obj$MSS.case)
colMeans(sim.obj$MBS.case)

# -----------------------------------------------------------------------------
model.file1 = "./jags/SparseCorr1_BSandSS_Reg1.jags"

jags.result = FitJags(sim.obj,
                      c("cell_prob", "Beta", "theta2", "theta2_value",
                        "bs_tpr", "ss_tpr", "bs_fpr"),
                      model.file1, NULL, hyper.pars.list,
                      6000, 4000, 2, 1)
save.image(file = "K5-Smax5-SC1-Reg2.RData")
coda.fit = as.mcmc(jags.result)

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
colMeans(Pr.Num.Path.fit[[1]])
colMeans(Pr.Num.Path.fit[[2]])

sim.obj$pars.baseline$Mu
sim.obj$pars.baseline$Pi
# ----------------------------------------------------------------------------
par.all = SweepPosterior(coda.fit[[1]], sim.obj, sweep.on = NULL,
                         model.type = "regression")
pairs(par.all[sample(1:1000, 200), ], pch = ".")

par.swept = SweepPosterior(coda.fit[[1]], sim.obj, sweep.on = "ss_tpr",
                           model.type = "regression")
pairs(par.swept, pch = ".")

C = cor(par.all)
corrplot.mixed(C, lower = "number", upper = "ellipse")
# -----------------------------------------------------------------------------
# result.table = rbind(
#   sim.obj$pars.baseline$Mu,
#   summary(stan.result$VB)$summary[1:5, 1],
#   apply(ExtractMu(as.mcmc(sc1.k5)[[1]], sim.obj), 2, mean),
#   apply(ExtractMu(as.mcmc(sc2.k5)[[1]], sim.obj), 2, mean),
#   apply(ExtractMu(as.mcmc(sc3.k5)[[1]], sim.obj), 2, mean)
# )
# result.table = cbind(result.table,
#                      apply(t(result.table) - result.table[1, ], 2, 
#                            function(x) sqrt(sum(x ^ 2))))
# colnames(result.table)[6] = "rmse"
# rownames(result.table) = c("true.val",
#                            "independent.vb",
#                            "sparsecorr1.mc",
#                            "sparsecorr2.mc",
#                            "sparsecorr3.mc")
# round(result.table, 3)


# --------------------------------------------------------------------------
library(Rcpp)
sourceCpp("./MakePrediction.cpp")
source("./MakePrediction.R")

pred.list = PredictCellProb(X.new = c(1, 1),
                            MBS.new = c(0, 1, 0, 0, 1),
                            MSS.new = c(0, 0, 0, 0, 1),
                            patho.names = NULL,
                            K = 5, Smax = 5,
                            X.train = sim.obj$X,
                            coda.chains = coda.fit[[1]])
pred.list$pred.Mu 
pred.list$pred.prob[1]
