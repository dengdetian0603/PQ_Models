setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/WorkSpace/SensitivityAnalysis")
suppressMessages(source("../../Code/PerchDataSimulation.R"))
suppressMessages(source("../../Code/PostProcessJAGS.R"))
suppressMessages(source("../../Code/PerchPlots.R"))

file.names = system("ls ./SC1-2*.csv", intern = TRUE)

# -----------------------------------------------------------------------------
f = 14
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

if (grepl("Baker", file)) {
  baker.fit = ListEtiology(cell.prob.fit, sim.obj,
                           etio.names = c("A", "B", "C", "D", "E"),
                           reorder = FALSE, num.keep = 16)
} else {
  cell.prob.fit0 = cell.prob.fit
}

PlotByPathogen(NULL, sim.obj, mu.fit = Mu.fit)
plot.obj = PlotByCombination(cell.prob.fit0, sim.obj,
                             hyper.pars.list, num.keep = 16,
                             has.true.value = TRUE,
                             contrast = "baker", baker.result = baker.fit)
do.call(grid.arrange, plot.obj)
# -----------------------------------------------------------------------------
err.tab = read.csv("../SC1_Reg_SensAnalysis.csv")
pr.nonindep.mean = with(err.tab, pind.a/(pind.a + pind.b))
pr.nonindep.sd = round(with(err.tab,
                            sqrt(pind.a * pind.b / (pind.a + pind.b) ^ 2
                                 / (pind.a + pind.b + 1))), 2)
err.tab = cbind(err.tab, pr.nonindep.mean, pr.nonindep.sd)[, -c(1, 4:11)]
err.tab.sc1 = err.tab[!is.na(err.tab$theta1.mu), ]
err.tab.baker =  err.tab[is.na(err.tab$theta1.mu), ]

g = ggplot()
g + geom_point(data = err.tab.sc1, aes(x = theta2.mu, y = Bhattacharyya,
                                       col = pr.nonindep.mean,
                                       size = pr.nonindep.sd)) +
  geom_hline(data = err.tab.baker, aes(yintercept = Bhattacharyya),
             col = "red")

g + geom_point(data = err.tab.sc1, aes(x = theta2.mu, y = expKL,
                                       col = pr.nonindep.mean,
                                       size = pr.nonindep.sd)) +
  geom_hline(data = err.tab.baker, aes(yintercept = expKL),
             col = "red")
