setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./IndepModel.R")
source("./SparseCorrModel.R")

par.default = SetDefaultSimulationParameter(4)
par.default$ncase = 2500
par.default$nctrl = 1000
par.default$phi.iter = 500
hyper.pars.list = SetDefaultHyperParameters(K = 5)
hyper.pars.list$pind_b = 2

sim.obj = do.call(SimulatePerchData, par.default)
sim.obj$pars.baseline$Mu
sim.obj$pars.baseline$Pi
sim.obj$pars.baseline$Thetas[-(1:5)]
sim.obj$pars.baseline$Thetas[1:5]

sim.obj.withNA = sim.obj
sim.obj.withNA$MBS.case[, 1] = NA
sim.obj.withNA$MSS.case[, 5] = NA

# -----------------------------------------------------------------------------
model.file0 = "./stan/Indep_BSandSSpos_NoReg_v2.stan"
model0 = stan_model(model.file0)
stan.result = FitStanNoReg(sim.obj, c("mu", "bs_tpr", "ss_tpr"),
                           model.file0, NULL, 
                           hyper.pars.list, method = c("VB", "MAP"))
stan.result$VB
stan.result$MAP$par

# -----------------------------------------------------------------------------
model.file1 = "./jags/SparseCorr1_BSandSS_NoReg.jags"
model.file2 = "./jags/SparseCorr2_BSandSS_NoReg.jags"
model.file3 = "./jags/SparseCorr3_BSandSS_NoReg.jags"
jags.result = FitJagsNoReg(sim.obj,
                           c("cell_prob", "theta1", "theta2",
                             "bs_tpr", "ss_tpr"),
                           model.file1, NULL, hyper.pars.list,
                           6000, 2000, 2, 1)

coda.fit = as.mcmc(jags.result)
post.mean = round(summary(coda.fit)[[1]][,1], 3)
post.mean[grepl("theta2", names(post.mean))]
post.mean[grepl("theta1", names(post.mean))]
plot(coda.fit)

Mu.fit = ExtractMu(coda.fit[[1]], sim.obj)
plot(as.mcmc(Mu.fit))
summary(Mu.fit)

# -----------------------------------------------------------------------------
result.table = rbind(
    sim.obj$pars.baseline$Mu,
    summary(stan.result$VB)$summary[1:5, 1],
    apply(ExtractMu(as.mcmc(sc1.k5)[[1]], sim.obj), 2, mean),
    apply(ExtractMu(as.mcmc(sc2.k5)[[1]], sim.obj), 2, mean),
    apply(ExtractMu(as.mcmc(sc3.k5)[[1]], sim.obj), 2, mean)
)
result.table = cbind(result.table,
                     apply(t(result.table) - result.table[1, ], 2, 
                           function(x) sqrt(sum(x ^ 2))))
colnames(result.table)[6] = "rmse"
rownames(result.table) = c("true.val",
                           "independent.vb",
                           "sparsecorr1.mc",
                           "sparsecorr2.mc",
                           "sparsecorr3.mc")
round(result.table, 3)

