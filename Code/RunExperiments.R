source("./PerchDataSimulation.R")
source("./IndepModel.R")

args = commandArgs(trailingOnly = TRUE)
method = 1
prefix.option = 0
if (length(args) > 1) {
  method = as.numeric(args[1])
}

if (method == 0.1) {
  par.to.save = "mu"
  model.file = "./jags/Indep_BSonly_FixPR_NoReg.txt"
} else if (method == 0.2) {
  par.to.save = "mu"
  model.file = "./jags/Indep_BSonly_FixPR_NoReg.txt"
  prefix.option = 1
} else if (method == 1) {
  par.to.save = c("mu", "bs_tpr")
  model.file = "./jags/Indep_BSandSSpos_NoReg.txt" 
} else if (method == 2) {
  par.to.save = c("mu", "bs_tpr")
  model.file = "./jags/Indep_BSandSS_NoReg.txt" 
}
print(paste0("Model: ", model.file))

registerDoMC(min(detectCores() - 1, 25))
sim.study.obj = SimStudyNoReg(par.to.save = par.to.save,
                              model.file = model.file,
                              bs.tpr.option = prefix.option,
                              n.iter = 5000, n.burnin = 2500,
                              n.thin = 2, n.rep = 100)
print(round(apply(sim.study.obj, 2, mean), 3))
save.image(file = paste0("SimStudy_NoReg_", method, ".Rdata"))