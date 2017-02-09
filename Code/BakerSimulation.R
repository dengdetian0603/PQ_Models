
library(baker)
library(lubridate)
library(rjags)
# ----------------------------------------------------------------------------

BakerSimulation <- function(with.pair = TRUE, is.regression = TRUE,
                            sim.obj, tid = 0,
                            n.itermcmc = 8000,
                            n.burnin   = 5000,
                            n.thin     = 3, save.img = FALSE,
                            has.MSS = c("A", "B", "C", "D", "E")) {
  #
  # 1.Clean Data:
  # specify cause list: top 5 and other
  cause_list <- c("NoA", "A", "B", "C", "D", "E")
  if (with.pair) {
    cause_list <- c("NoA", "A", "B", "C", "D", "E",
                    "A+B", "A+C", "A+D", "A+E",
                    "B+C", "B+D", "B+E",
                    "C+D", "C+E",
                    "D+E")
  }
  # specify measurements:
  patho_BrS_NPPCR = c("A", "B", "C", "D", "E")
  patho_BrS_NPPCR = intersect(patho_BrS_NPPCR, cause_list)

  patho_SS_BCX = has.MSS
  patho_SS_BCX = intersect(patho_SS_BCX, cause_list)

  ## parent directory for testing code (JHPCE):
  working_dir <- "/users/ddeng/ThesisTopic/PQModels/Experiments/Baker/"

  ### Set Sensitivity Bound
  sens_NPPCR_low_noabx = rep(0.5, 5)
  sens_NPPCR_up_noabx = rep(0.9, 5)
  sens_BCX_low_noabx = rep(0.05, 5)
  sens_BCX_up_noabx = rep(0.2, 5)

  ## clean PERCH data:
  data_nplcm.sim <- ReorgToNplcmData(sim.obj)

  #
  # Specify Model: 
  #           
  no_abx_strat <- rep(1, sum(data_nplcm.sim$Y)) # no TPR stratification

  m_opt1 <- list(likelihood = list(cause_list = cause_list,
                                   k_subclass = c(1,1),
                                   Eti_formula = ~ AGE + SEVERITY +
                                                   AGE * SEVERITY,
                                   FPR_formula = list(NPPCR = ~ 1)), 
                 use_measurements = c("BrS","SS"),
                 prior = list(Eti_prior = overall_uniform(0.1,cause_list),                     
                              TPR_prior = list(
                                BrS = list(
                                  info  = "informative",
                                  input = "match_range",
                                  grp = no_abx_strat, 
                                  val = list(
                                    NPPCR = list(
                                      up = list(sens_NPPCR_up_noabx), 
                                      low = list(sens_NPPCR_low_noabx)))),
                                SS = list(
                                  info = "informative",
                                  input = "match_range",
                                  grp = no_abx_strat,
                                  val = list(
                                    BCX = list(up = list(sens_BCX_up_noabx), 
                                               low = list(sens_BCX_low_noabx))
                                    )))))   
  model_options <- m_opt1
  if (is.regression) {
    model_options$prior$Eti_prior = 2 # for regression use
  } else {
    model_options$likelihood$Eti_formula <- ~ 1 # for non-regression use
  }

  suppressMessages(assign_model(model_options, data_nplcm.sim))
  #
  # 4. Fit model:
  #
  add.name1 = "_Singleton"
  if (with.pair) {
    add.name1 = "_SingletonAndPair"
  }
  add.name2 = "_NonReg"
  if (is.regression) {
    add.name2 = "_Reg"
  }

  # date stamp for analysis:
  Date <- gsub("-", "_", Sys.Date())

  dated_name <- paste0(working_dir, Date, "_Top5", add.name1, add.name2, "_t-", tid)
  if (!dir.exists(dated_name)) {
    dir.create(dated_name)
  }
  fullname <- dated_name

  result_folder <- fullname

  # options for MCMC chains:
  mcmc_options <- list(debugstatus = TRUE,
                       n.chains   = 1,
                       n.itermcmc = n.itermcmc,
                       n.burnin   = n.burnin,
                       n.thin     = n.thin,
                       individual.pred = !TRUE,
                       ppd             = !TRUE,
                       get.pEti        = TRUE,
                       result.folder = result_folder,
                       bugsmodel.dir = result_folder,
                       jags.dir = "",
                       use_jags = TRUE)

  rjags::load.module("glm")
  gs <- nplcm(data_nplcm.sim, model_options, mcmc_options)
  if (save.img) {
    save.image(file = paste0(fullname, "/workspace.RData"))
  }

  return(list(gs = gs, cause_list = cause_list))
}