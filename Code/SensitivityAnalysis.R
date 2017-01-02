setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/WorkSpace/SensitivityAnalysis")
suppressMessages(source("../../Code/PerchDataSimulation.R"))
suppressMessages(source("../../Code/PostProcessJAGS.R"))
suppressMessages(source("../../Code/PerchPlots.R"))

file.names = system("ls ./SC1-2*.csv", intern = TRUE)


f = 1
file = file.names[f]
sim.fit = read.csv(file)
load(gsub(".csv", ".RData", file))

EtioPrior = ListEtiologyPriorSC1(5, 5, hyper.pars.list,
                                 as.character(1:5), 20000)

Mu.fit = list() 
cell.prob.fit = NULL
n.strata = 4
for (i in 1:n.strata) {
  temp = sim.fit[, paste0(paste0("Mu.", i), ".", 1:5, ".")]
  colnames(temp) = paste0("Mu[", 1:5, "]")
  Mu.fit[[i]] = temp
  # sim.fit[, paste0(paste0("PrNumPath.", i), ".", 0:5, ".")]
  temp = sim.fit[, paste0(paste0("cell_prob.", i), ".", 1:32, ".")]
  colnames(temp) = paste0(paste0("cell_prob[", i), ",", 1:32, "]")
  if (length(cell.prob.fit) > 0) {
    cell.prob.fit = cbind(cell.prob.fit, temp)    
  } else {
    cell.prob.fit = temp
  }
}

PlotByPathogen(NULL, sim.obj, mu.fit = Mu.fit)
plog.obj = PlotByCombination(cell.prob.fit, sim.obj,
                             hyper.pars.list, num.keep = 16,
                             has.true.value = TRUE)