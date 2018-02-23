setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")
source("./ModelFittingVB.R")
source("../../BAKER/utils.R")

library(gtools)
library(dplyr)
library(ggplot2)
d.mat = DesignMatrixAppxQuadExp(5, 5)
L0.mat = rbind(rep(0, 5), d.mat$Lmat)

x.label = apply(L0.mat, 1, function(x) {
  label = ""
  for (xi in x) {
    if (xi > 0) {
      label = paste(label, "1")
    } else {
      label = paste(label, "0")
    }
  }
  label
})

# -----------------------------------------------------------------------------
hyper.pars.list = SetVBHyperParameters(K = 5)
hyper.pars.list$aa = rep(3.24, 5)  # \in (0.01, 0.1)
hyper.pars.list$bb = rep(71.08, 5)
# beta_parms_from_quantiles(q = c(0.01, 0.1))
hyper.pars.list$cc = 22.03  # \in (0.6, 0.9)
hyper.pars.list$dd = 6.67

# -----------------------------------------------------------------------------
ss.tpr = c(0.05, 0.06, 0.04, 0.055, 0.045)

bs.tpr = c(0.85, 0.65, 0.75, 0.7, 0.65)
bs.fpr = c(0.45, 0.3, 0.35, 0.4, 0.35) - 0.2

theta1 = c(-1, -0.5, 0, 1, 1.5) - 0.5
rho.D. = c(-0, -1, -2, -1,
           -3, -1, -2,
           -1, -2,
           -3) * 1.5
theta2 = round(rdirichlet(1, 1 - rho.D.) * (-20), 1)
theta2

sim.obj = SimulateNoRegData(ncase = 400, nctrl = 1000, 
                            theta1 = theta1,
                            theta2 = theta2,
                            ss.tpr = ss.tpr,
                            bs.tpr = bs.tpr,
                            bs.fpr = bs.fpr)
sim.obj$Pr.NumPathogen
sim.obj$Mu

input.obj = c(sim.obj, list(ss.available = 1:5, bs.available = 1:5))
#------------------------------------------------------------------------------
hyper.pars.list$theta_mu = c(0, 0, 0, 0, 0) + 1
hyper.pars.list$theta_tau = 3

hyper.pars.list$rho_mu = -6
hyper.pars.list$rho_tau = 15

hyper.pars.list$pind_a = 4 * 3
hyper.pars.list$pind_b = 1 * 3

par.list = SetVBInitialValues(5, input.obj, hyper.pars.list)

par.vbfit = FitVBEMnoReg(input.obj, hyper.pars.list, par.list,
                         max.iter = 450, tol = 1e-6) 

# -----------------------------------------------------------------------------
par.vbfit$mu_theta
par.vbfit$mu_rho  
round(par.vbfit$mu_rho * par.vbfit$qD, 3)

etio.info = VBEtioProbsNoReg(par.vbfit)

round(cbind(sim.obj$Mu, etio.info$etio.mean.pL), 3)
round(cbind(sim.obj$Pr.NumPathogen, etio.info$etio.numPi.pL), 3)

# -----------------------------------------------------------------------------
plot(sim.obj$cell.probs, type = "n", ylim = c(0, 0.5),
     xaxt = "n", xlab = "", ylab = "Probability")
axis(1, at = 1:32, labels = x.label, las = 2)
mtext(text = c("A B C D E"), side = 2, at = -0.09)

by.priors = 
  rdirichlet(100, (etio.info$etio.probs.pL[, 1] + sim.obj$cell.probs) * 70)
for (i in 1:100) {
  points(jitter(1:32), by.priors[i, ],
         pch = 16, col = alpha(rgb(0, 0, 0), 0.4), cex = 1.2)
}
lines((etio.info$etio.probs.pL[, 1] + sim.obj$cell.probs)/2,
      type = "b", pch = "-", col = "purple", cex = 4, lwd = 2)
points(sim.obj$cell.probs, col = "red", cex = 5, pch = "-")
