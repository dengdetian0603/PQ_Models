setwd("~/Documents/JHSPH/Research/S.Zeger/PQ_Models/Code")
source("./PerchDataSimulation.R")

theta1 = c(1, -0.5, 0.1)
theta2 = c(-0.5, -1, -2)

d.mat = DesignMatrixAppxQuadExp(3, 3)
lu.mat = rbind(rep(0, 6), cbind(d.mat$Lmat, d.mat$Umat))

potentials = exp(lu.mat %*% cbind(c(theta1, theta2)))
probs = potentials/sum(potentials)

prob.mu = d.mat$MuMat %*% probs[-1, ]; prob.mu
prob.pi = d.mat$PiMat %*% probs[-1, ]; prob.pi

# ------------------------------------------------------------------------------
n.sim = 5000
theta2.0 = -2
d.vec = c(0.25, 0.5, 0.99)
prob.mu.mat = matrix(NA, nrow = 3, ncol = n.sim)
prob.pi.mat = matrix(NA, nrow = 3, ncol = n.sim)

theta2.mat = cbind(rbinom(n.sim, 1, d.vec[1]),
                   rbinom(n.sim, 1, d.vec[2]),
                   rbinom(n.sim, 1, d.vec[3])) * theta2.0

x.mix = c()
for (i in 1:n.sim) {
  potentials = exp(lu.mat %*% cbind(c(theta1, theta2.mat[i, ])))
  # probs = potentials/sum(potentials)
  # prob.mu.mat[, i] = d.mat$MuMat %*% probs[-1, ]
  # prob.pi.mat[, i] = d.mat$PiMat %*% probs[-1, ]
  
  x = rmultinom(1, 1, potentials)
  prob.mu.mat[, i] = d.mat$MuMat %*% x[-1, ]
  prob.pi.mat[, i] = d.mat$PiMat %*% x[-1, ]
  
  x.mix[i] = which(x[, 1] > 0)
}
rowMeans(prob.mu.mat)
rowMeans(prob.pi.mat)

cbind(probs, table(x.mix)/n.sim)

plot(density(prob.mu.mat[, 1]))

# ==============================================================================
d.mat = DesignMatrixAppxQuadExp(3, 3)
lu.mat = rbind(rep(0, 6), cbind(d.mat$Lmat, d.mat$Umat))

#########################
theta1 = runif(3, -1, 1)
theta2 = -rgamma(3, 1, 1)
###########################

potentials = exp(lu.mat %*% cbind(c(theta1, theta2)))
probs = potentials/sum(potentials)

prob.mu = d.mat$MuMat %*% probs[-1, ]; prob.mu
prob.pi = d.mat$PiMat %*% probs[-1, ]; prob.pi

# ------------------------------------------------------------------------------
n.sim = 5000
theta2.0 = min(theta2) - 1e-3
d.vec = theta2/theta2.0
prob.mu.mat = matrix(NA, nrow = 3, ncol = n.sim)
prob.pi.mat = matrix(NA, nrow = 3, ncol = n.sim)

theta2.mat = cbind(rbinom(n.sim, 1, d.vec[1]),
                   rbinom(n.sim, 1, d.vec[2]),
                   rbinom(n.sim, 1, d.vec[3])) * theta2.0

x.mix = c()
for (i in 1:n.sim) {
  potentials = exp(lu.mat %*% cbind(c(theta1, theta2.mat[i, ])))

  x = rmultinom(1, 1, potentials)
  prob.mu.mat[, i] = d.mat$MuMat %*% x[-1, ]
  prob.pi.mat[, i] = d.mat$PiMat %*% x[-1, ]
  
  x.mix[i] = which(x[, 1] > 0)
}

cbind(probs, table(x.mix)/n.sim)
sum(sqrt(probs[, 1] * table(x.mix)/n.sim))


