setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/WorkSpace/")

suppressMessages(source("../Code/PerchDataSimulation.R"))
suppressMessages(source("../Code/PostProcessJAGS.R"))
suppressMessages(source("../Code/PerchPlots.R"))

library(ICC)

file.names = system("ls ./AD.SensitivityAnalysis/*.csv", intern = TRUE)
# -----------------------------------------------------------------------------
f = 3
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
err.tab = read.csv("./SC1_Reg_SensAnalysis.csv")
write.csv(round(err.tab[,c(3:5, 14)], 4), file = "temp.csv")

pr.nonindep.mean = with(err.tab, pind.a/(pind.a + pind.b))
pr.nonindep.sd = round(with(err.tab,
                            sqrt(pind.a * pind.b / (pind.a + pind.b) ^ 2
                                 / (pind.a + pind.b + 1))), 2)
err.tab = cbind(err.tab, pr.nonindep.mean, pr.nonindep.sd)[, -c(1, 4:11)]
err.tab.sc1 = err.tab[!is.na(err.tab$theta1.mu), ]
err.tab.baker =  err.tab[is.na(err.tab$theta1.mu), ]
err.tab.baker$Specification = c("Singletons", "Singletons+Pairs")


g = ggplot()
g + geom_point(data = err.tab.sc1, aes(x = theta2.mu, y = Bhattacharyya,
                                       col = pr.nonindep.mean,
                                       size = pr.nonindep.sd)) +
  labs(x = "Prior Mean of theta2") + ylim(0.75, 1) +
  scale_color_continuous(name = "Prior Mean of p") +
  scale_size_continuous(name = "Prior SD of p") +
  geom_hline(data = err.tab.baker, aes(yintercept = Bhattacharyya,
                                       linetype = Specification),
             col = "red")

# ----------------------------------------------------------------------------
n.strata = 4
prior.cat = 1
prior.name = NULL
pooled.sim.fit = NULL
for (f in 1:length(file.names)) {
  file = file.names[f]
  sim.fit = read.csv(file)
  if (!"tid" %in% colnames(sim.fit)) {
    next
  }
  sim.fit$prior.cat = prior.cat
  prior.name = c(prior.name, file)
  prior.cat = prior.cat + 1
  
  pooled.sim.fit  = rbind(pooled.sim.fit, sim.fit)
}
prior.cat.0 = pooled.sim.fit$prior.cat 
pooled.sim.fit0 = pooled.sim.fit

# pooled.sim.fit$prior.cat = prior.cat.0
# o111 = order(pooled.sim.fit[pooled.sim.fit$tid == 50, "Mu.1.1."])
# prior.cat.new = o111[pooled.sim.fit$prior.cat]
# pooled.sim.fit$prior.cat = prior.cat.new

pooled.sim.fit = pooled.sim.fit0[pooled.sim.fit0$tid < 20, ]
gp = list()
gp[[1]] = ggplot(data = pooled.sim.fit,
                 aes(x = prior.cat, y = Mu.1.1., group = tid)) +
          geom_jitter(alpha = 0.5, width = 0.1) + 
          geom_line(aes(col = tid), alpha = 0.3)
gp[[2]] = ggplot(data = pooled.sim.fit,
                 aes(x = prior.cat, y = Mu.1.2., group = tid)) +
  geom_jitter(alpha = 0.5, width = 0.1) + 
  geom_line(aes(col = tid), alpha = 0.3)
gp[[3]] = ggplot(data = pooled.sim.fit,
                 aes(x = prior.cat, y = Mu.1.3., group = tid)) +
  geom_jitter(alpha = 0.5, width = 0.1) + 
  geom_line(aes(col = tid), alpha = 0.3)
gp[[4]] = ggplot(data = pooled.sim.fit,
                 aes(x = prior.cat, y = Mu.1.4., group = tid)) +
  geom_jitter(alpha = 0.5, width = 0.1) + 
  geom_line(aes(col = tid), alpha = 0.3)
gp[[5]] = ggplot(data = pooled.sim.fit,
                 aes(x = prior.cat, y = Mu.1.5., group = tid)) +
  geom_jitter(alpha = 0.5, width = 0.1) + 
  geom_line(aes(col = tid), alpha = 0.3)

do.call(grid.arrange, gp)
# ICCest(x = tid, y = Mu.2.1., data = pooled.sim.fit0)$ICC

# ---------------------------------------------------------------------------
d1 = pooled.sim.fit[pooled.sim.fit$prior.cat == 2, c("tid", "Mu.1.1.")]
d2 = pooled.sim.fit[pooled.sim.fit$prior.cat == 4, c("tid", "Mu.1.1.")]

d1 = d1[order(d1$tid), ]
d2 = d2[order(d2$tid), ]
i.set = intersect(d1$tid, d2$tid)
plot(d1$Mu.1.1.[d1$tid %in% i.set], d2$Mu.1.1.[d2$tid %in% i.set],
     xlab = "Prior Category 1", ylab = "Prior Category 2",
     xlim = c(0, 0.6), ylim = c(0, 0.6))
abline(a = 0, b = 1, col = "red")

prior.name

