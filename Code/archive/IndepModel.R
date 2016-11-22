library('coda')
library('data.table')
library('R2jags')
library('rjags')
library('rstan')

# -----------------------------------------------------------------------------
FitMLEgivenGSNoReg <- function(K, sim.obj) {
# Get bs.tpr and mu estimtates assuming we know GS data
  bs.tpr.mle = rep(NA, K)
  bs.tpr.hat = rep(NA, K)
  for (j in 1:K) {
    bs.tpr.hat[j] =
        sum(sim.obj$MBS.case[, j] * sim.obj$MSS.case[, j]) /
        sum(sim.obj$MSS.case[, j])
    bs.tpr.mle[j] =
        sum(sim.obj$MBS.case[, j] * sim.obj$L[, j]) /
        sum(sim.obj$L[, j])
  }  
  mu.mle = rep(NA, K)
  mu.mle = apply(sim.obj$L, 2, mean)
  data.frame(mu.mle, bs.tpr.mle, bs.tpr.hat)
}

SetDefaultHyperParameters <- function(K = 5, logit.mu.rho = -0.22) {
  V = matrix(logit.mu.rho, nrow = K, ncol = K)
  diag(V) = 1
  V = 1 * V
  list(lambda = 1.4,
       alpha = rep(1, K),
       mu_logit = rep(-0.5, K), V = V,
       ma = 2, mb = 2, 
       aa = 1, bb = 9,
       cc = 4, dd = 2,
       ee = 1, ff = 1,
       tau_theta = 1/10,
       pind_a = 2, pind_b = 8,
       pind_dirich_par = c(5, 1, 5),
       theta2_cat = c(0, 0.5, -1))
}

# FitJagsNoReg <- function(sim.obj, parameters.to.save,
#                          model.file, bs.tpr.prefix = NULL,
#                          hyper.pars.list = SetDefaultHyperParameters(),
#                          n.iter = 6000, n.burnin = 4000,
#                          n.thin = 2, n.chains = 2) {
# # Fit JAGS non-regression models
#   sim.dat <- c(list(
#       # data
#       N_ctrl = nrow(sim.obj$MBS.ctrl),
#       K = ncol(sim.obj$MBS.ctrl),
#       mbs_ctrl = sim.obj$MBS.ctrl,
#       N_case = nrow(sim.obj$MBS.case),
#       mbs_case = sim.obj$MBS.case,
#       mss_case = sim.obj$MSS.case,
#       # hyper parameters
#       bs_tpr_prefix = bs.tpr.prefix), hyper.pars.list)
#   if (n.chains > 1) {
#     mc.fit <- jags.parallel(
#                   data = sim.dat,
#                   #inits = bayes.mod.inits,
#                   parameters.to.save = parameters.to.save,
#                   n.chains = n.chains,
#                   n.iter = n.iter,
#                   n.burnin = n.burnin,
#                   n.thin = n.thin,
#                   model.file = model.file)
#   } else {
#     mc.fit <- jags(
#                   data = sim.dat,
#                   #inits = bayes.mod.inits,
#                   parameters.to.save = parameters.to.save,
#                   n.chains = n.chains,
#                   n.iter = n.iter,
#                   n.burnin = n.burnin,
#                   n.thin = n.thin,
#                   model.file = model.file)    
#   }
#   print(mc.fit)
#   return(mc.fit)
# }

FitStanNoReg <- function(sim.obj, parameters.to.save,
                         model.file, bs.tpr.prefix = NULL,
                         hyper.pars.list = SetDefaultHyperParameters(),
                         method = c("HMC", "VB", "MAP"),
                         n.iter = 3000, n.burnin = 1000,
                         n.thin = 2, n.chains = 3, n.cores = 3) {
# Fit STAN non-regression models
#
# Example:
#   model.file = "./stan/Indep_BSandSSpos_NoReg.stan"
#   model.file = "./stan/Indep_LogitPrior_BSandSSpos_NoReg.stan"
#   par.default = SetDefaultSimulationParameter(1)
#   sim.obj = do.call(SimulatePerchData, par.default)
#   stan.result = FitStanNoReg(sim.obj, "mu", model.file, method = c("VB", "MAP"))
#
  sim.dat <- c(list(
    # data
    N_ctrl = nrow(sim.obj$MBS.ctrl),
    K = ncol(sim.obj$MBS.ctrl),
    mbs_ctrl = sim.obj$MBS.ctrl,
    N_case = nrow(sim.obj$MBS.case),
    mbs_case = sim.obj$MBS.case,
    mss_case = sim.obj$MSS.case,
    # hyper parameters
    bs_tpr_prefix = bs.tpr.prefix),
    hyper.pars.list
  )
  # Fit STAN moel
  stan.model = stan_model(model.file)
  if ("HMC" %in% method) {
    mc.fit = tryCatch(stan(
        file = model.file, data = sim.dat, pars = parameters.to.save,
        chains = n.chains, iter = n.iter, warmup = n.iter - n.burnin,
        thin = n.thin, cores = n.cores),
        error = function(c) "error in HMC sampling.")
  } else {
    mc.fit = NULL
  }
  if ("VB" %in% method) {
    vb.fit = tryCatch(vb(
        object = stan.model, data = sim.dat, pars = parameters.to.save,
        algorithm = "meanfield", tol_rel_obj = 0.0005),
        error = function(c) "error in VB approximation.")
  } else {
    vb.fit = NULL 
  }
  if ("MAP" %in% method) {
    map.fit = tryCatch(optimizing(
        object = stan.model, data = sim.dat, 
        algorithm = "BFGS"),
        error = function(c) "error in Maximizig Posterior.")
  } else {
    map.fit = NULL
  }
  return(list(HMC = mc.fit, VB = vb.fit, MAP = map.fit))
}



SimStudyNoRegJAGS <- function(sim.pars = SetDefaultSimulationParameter(1),
                              par.to.save = c("mu"),
                              model.file,
                              bs.tpr.option = 0,
                              n.iter, n.burnin, n.thin,
                              n.rep, max.core = 40) {
# Simulation study for JAGS non-regression models
# Example:
#   SimStudyNoRegJAGS(model.file = "./jags/Indep_BSandSSpos_NoReg.txt",
#                     n.iter = 6000, n.burnin = 3000, n.thin = 3, n.rep = 3)
#
  registerDoMC(min(detectCores() - 1, max.core))
  sim.obj0 = do.call(SimulatePerchData, sim.pars)
  cell.prob.unique  = sim.obj0$cell.prob.unique
  sim.pars$cell.prob.unique = cell.prob.unique
  mc.fit.all = foreach(i = 1:n.rep, .combine = rbind) %dopar% {
    print(paste0(i, "th rep: start."))
    set.seed(i*123)
    sim.obj = do.call(ReSimulateData, sim.pars)
    direct.fit = FitMLEgivenGSNoReg(sim.pars$K, sim.obj)
    if (bs.tpr.option == 0) {
      bs.tpr.prefix = sim.pars$bs.tpr
    } else {
      bs.tpr.prefix = direct.fit$bs.tpr.hat
    }
    print(paste0(i, "th rep: sampling..."))
    mc.fit = FitJagsNoReg(sim.obj, par.to.save,
                          model.file,
                          bs.tpr.prefix,
                          n.iter = n.iter, n.burnin = n.burnin,
                          n.thin = n.thin, n.chains = 1)
    coda.fit = as.mcmc(mc.fit)
    post.mean = summary(coda.fit)[[1]][,1]
    print(paste0(i, "th rep: end."))
    c(post.mean, direct.fit$bs.tpr.hat)
  }
  mc.fit.all
}


SimStudyNoRegSTAN <- function(sim.obj0 = NULL,
                              sim.pars = SetDefaultSimulationParameter(4),
                              par.to.save = c("mu"),
                              model.file,
                              method = "MAP",
                              bs.tpr.option = 0,
                              n.iter = 3000, n.burnin = 1000, n.thin = 2,
                              n.rep = 200, max.core = 40, seed0 = 123) {
  # Simulation study for JAGS non-regression models
  # Example:
  #   SimStudyNoRegSTAN(model.file = "./stan/Indep_BSandSSpos_NoReg.stan",
  #   par.to.save = c("mu", "bs_tpr"), method = "MAP", n.rep = 300) -> sim.study.obj
  registerDoMC(min(detectCores() - 1, max.core))
  if (length(sim.obj0) < 1) {
    set.seed(seed0)
    sim.obj0 = do.call(SimulatePerchData, sim.pars)
  }
  print(round(sim.obj0$pars.baseline$Mu, 3))
  cell.prob.unique  = sim.obj0$cell.prob.unique
  sim.pars$cell.prob.unique = cell.prob.unique
  stan.model = stan_model(model.file)
  stan.fit.all = 
      foreach(i = 1:n.rep, .combine = rbind, .inorder = FALSE) %dopar% {
  # ------------------------------------  
    print(paste0(i, "th rep: start."))
    set.seed(i*123)
    sim.obj = do.call(ReSimulateData, sim.pars)
    direct.fit = FitMLEgivenGSNoReg(sim.pars$K, sim.obj)
    if (bs.tpr.option == 0) {
      bs.tpr.prefix = sim.pars$bs.tpr
    } else {
      bs.tpr.prefix = direct.fit$bs.tpr.hat
    }
  # -------------------------------------------
    print(paste0(i, "th rep: model fitting..."))
    stan.result = FitStanNoReg(
        sim.obj, par.to.save, model.file, method = method,
        n.iter = 3000, n.burnin = 1000,
        n.thin = 2, n.chains = 1, n.cores = 1
        )
    if (method == "VB" & class(stan.result$VB) != "character") {
      par.hat = summary(stan.result$VB)$summary[,1]
    } else if (method == "MAP" & class(stan.result$MAP) != "character") {
      par.hat = stan.result$MAP$par
    } else {
      par.hat = rep(NA, length(par.to.save))
      names(par.hat) = par.to.save
    }
    to.save.index = 
      apply(sapply(par.to.save, grepl, names(par.hat)), 1, sum) > 0
    par.hat = par.hat[to.save.index]
    sim.result = data.frame(method = rep(method, length(par.hat)),
                            par.name = names(par.hat),
                            par.est = par.hat)
  # --------------------------------------------  
    print(paste0(i, "th rep: end."))
    sim.result
  }
  list(stan.fit.all = stan.fit.all, sim.obj0 = sim.obj0 , sim.pars = sim.pars)
}



