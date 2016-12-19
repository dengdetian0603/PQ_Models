setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingJAGS.R")
source("./PostProcessJAGS.R")
source("./PerchPlots.R")
library(coda)

ws.dir = "~/Documents/JHSPH/Research/S.Zeger/PQ_Models/WorkSpace/"

# ------------------- No regression ------------------------------------------
load(paste0(ws.dir, "2016_12_19_Top5_Singleton_NonReg.RData"))
load(paste0(ws.dir, "2016_12_19_Top5_SingletonAndPair_NonReg.RData"))

par.fit = as.mcmc(gs)
pEti.est = data.frame(Eti.Name = cause_list,
                      pEti = par.fit$mean$pEti)[order(par.fit$mean$pEti,
                                                      decreasing = TRUE), ]
# singleton
pEti.byPathogen.matrix = par.fit$sims.matrix[, paste0("pEti[", 1:6,"]")]
colnames(pEti.byPathogen.matrix) = c(cause_list[1:5], "Other")

fit1 = pEti.byPathogen.matrix

# with pair
pEti.fit = par.fit$sims.matrix[, paste0("pEti[", 1:16,"]")]
singleton.Eti = cause_list[1:5]
pEti.byPathogen.matrix = matrix(0, nrow = nrow(pEti.fit), ncol = 5)
for (i in 1:5) {
  temp = pEti.fit[, grepl(singleton.Eti[i], cause_list)]
  pEti.byPathogen.matrix[, i] = rowSums(temp)
}
pEti.byPathogen.matrix = cbind(pEti.byPathogen.matrix, pEti.fit[, 6])
colnames(pEti.byPathogen.matrix) = c(singleton.Eti, "Other")

fit2 = pEti.byPathogen.matrix

PlotCompareResult(fit1, fit2, c("Singletons", "Singletons&Pairs"))

# ----------------------- Regression ------------------------------------------
load(paste0(ws.dir, "2016_12_19_Top5_Singleton_Reg.RData"))
load(paste0(ws.dir, "2016_12_19_Top5_SingletonAndPair_Reg.RData"))

par.fit = as.mcmc(gs)

Z.unique = uniquecombs(Z_Eti00)
DIM = dim(par.fit$sims.list$pEti)

pEti.array = array(as.vector(par.fit$sims.list$pEti), c(DIM[1], DIM[3], DIM[2]))
pEti.list = list()
baker.result = list()
for (i in 1:nrow(Z.unique)) {
  temp = pEti.array[, attr(Z.unique, "index") == i, ]
  pEti.list[[i]] = matrix(0, nrow = dim(temp)[1] * dim(temp)[2],
                          ncol = dim(temp)[3])
  for (j in 1:dim(temp)[2]) {
    pEti.list[[i]][(1 + dim(temp)[1] * (j - 1)):(dim(temp)[1] * j), ] =
      temp[, j, ]
  }
  baker.result[[i]] = data.frame(EtioComb = cause_list,
                                 Probability = colMeans(pEti.list[[i]]),
                                 Prob.lower = apply(pEti.list[[i]], 2,
                                                    quantile, 0.025),
                                 Prob.upper = apply(pEti.list[[i]], 2,
                                                    quantile, 0.975))
}
baker.result[[1]]
save(pEti.list, baker.result, Z.unique,
     file = "../WorkSpace/baker_singleton.RData")

# ----------------------------------------------------------------------------
load("../WorkSpace/Kenya_sc1.RData")
load("../WorkSpace/baker_singleton_pair.RData")
load("../WorkSpace/baker_singleton.RData")
plot.obj = 
  PlotByCombination(coda.fit[[1]], sim.obj, hyper.pars.list,
                    etio.names = c("RSV", "RHINO", "HMPV_A_B", "PNEU", "PV_EV"),
                    "baker", FALSE, 16, baker.result)
do.call(grid.arrange, plot.obj)


