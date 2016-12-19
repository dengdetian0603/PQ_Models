setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingJAGS.R")
source("./PostProcessJAGS.R")
source("./PerchPlots.R")
library(coda)

ws.dir = "~/Documents/JHSPH/Research/S.Zeger/PQ_Models/WorkSpace/"

# ------------------- No regression ------------------------------------------
load(paste0(ws.dir, "2016_12_18_Top5_Singleton_NonReg.RData"))
load(paste0(ws.dir, "2016_12_18_Top5_SingletonAndPair_NonReg.RData"))

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
load(paste0(ws.dir, "2016_12_18_Top5_Singleton_Reg.RData"))
load(paste0(ws.dir, "2016_12_18_Top5_SingletonAndPair_Reg.RData"))

par.fit = as.mcmc(gs)


