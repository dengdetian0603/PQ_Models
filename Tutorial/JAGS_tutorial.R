library('R2jags')
library('rjags')

# Simulate data
N <- 1000
x <- 1:N
z <- 0.01 * x - 5
y <- sapply(1 / (1 + exp(-z)), function(p) {rbinom(1, 1, p)})
 
write.table(data.frame(X = x, Z = z, Y = y),
            file = 'example3.data',
            row.names = FALSE,
            col.names = TRUE)

sim.dat <- list('x' = x, 'y' = y, 'N' = N)

# Fit jags model
# use rjags lib
model.logistic.regression <- jags.model('./eg_logistic_regression.bug',
                                        data = sim.dat,
                                        n.chains = 3,
                                        n.adapt = 100)

# use R2jags lib
bayes.mod.params <- c("a", "b")
bayes.mod.inits <- function() {
  list("a" = rnorm(1), "b" = rnorm(1))
}
bayes.mod.fit <- jags(data = sim.dat,
                      #inits = bayes.mod.inits,
                      parameters.to.save = bayes.mod.params,
                      n.chains = 3,
                      n.iter = 2000,
                      n.burnin = 500,
                      model.file = './eg_logistic_regression.bug')
bayes.mod.fit.upd <- autojags(bayes.mod.fit)

print(bayes.mod.fit)
print(bayes.mod.fit.upd)

pdf("./bayes_trace.pdf")
traceplot(bayes.mod.fit)

dev.off()

