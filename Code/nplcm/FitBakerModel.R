###############################################################################
# Streamlined code for PQ regression
# Zhenke Wu | zhenkewu@gmail.com
# First version: April 27, 2016
##############################################################################
# --------------------------
# 1. Regression functionalities:
#  Change the "model_options$likelihood" below to do various regressions:
# 
# a. seasonality + age category + very severe for etiology regression,
#    seasonality + age category for false positive regression can be done by:
#
#  Eti_formula = ~ s_date_Eti(ENRLDATE,Y,basis = "ncs",5)+ as.factor(AGECAT>1)+as.factor(CXR_VS)
#  FPR_formula = list(
#   NPPCR = ~ s_date_FPR(ENRLDATE,Y,basis="ps",10)+as.factor(AGECAT>1)
#  )
#
#  Here, there is only one etiology regression among cases, but could be many 
#  FPR regressions (when using multiple BrS data), in which case you need to specify separate 
#  regression models for each measurement slice, e.g., 
#   
#   FPR_formula=list(NPPCR=~X1+X2, WBPCR=~X2+X3).
#
# b. Currently only additive models (~X1+X2, not ~X1+X2+X1*X2) are gauranteed 
#   to work. For models with interaction terms, please check again with Zhenke.
# 
# c. If you want to fit all sites together, in the raw data PQXXXXXX.csv, put a 
#  column named "include" and specify all values to be 1. Then in `clean_options`
#  specify X_strat="include", and X_strat_val=list(1).
#
# -------------------------
# 2. User-specified weights for pie-standardization.
#  When the covariates for etiology regression are all discrete, one can do user-
# specified standardization by using a vector of weights, each corresponding to 
# one level of the covariates. For example, if AGECAT>1 and AGECAT<=1 are two levels,
# then one can specify the weights as c(0.4,0.6), which could be different from
# the proportion observed in the actual sample.
#
# --------------------------
# 3. Multiple Data Sources:
#
# Can work with silver-standard data (need informative priors on sensitivities).
# Simply by constructing SS_object_1,SS_object_2, etc to feed into `clean_data_options`.
#
# -------------------------
# 4. Software
#  Regression functions currently only works with JAGS 4.x.x and R package `rjags 4-6`. 
##############################################################################

rm(list=ls())
library(baker)
library(lubridate)
library(rjags)
#
# 1.Clean Data:

# specify cause list:
cause_list <- c("BOPE","C_PNEU","M_PNEU",
                "PCP","ADENO","CMV","COR",
                "FLU_C","HBOV","HMPV_A_B",
                "FLU_A","FLU_B",
                "PARA_1","PARA_2","PARA_3","PARA_4",
                "PV_EV","RHINO","NoA",
                "HINF","MCAT","PNEU",
                "SASP","SAUR",
                "ENTRB","NFGNR","FUNGI","HAEMO","OTHSTR","NMEN","TB",#SS only 
                "RSV"
)
# specify measurements:
patho_BrS_NPPCR      <- c("BOPE","C_PNEU","M_PNEU",
                          "PCP","ADENO","CMV","COR",#WF
                          "FLU_C","HBOV","HMPV_A_B",
                          "FLU_A","FLU_B",
                          "PARA_1","PARA_2","PARA_3","PARA_4",
                          "PV_EV","RHINO","RSV",#without SS measure.
                          "HINF",
                          #"MCAT",
                          "PNEU",
                          "SASP","SAUR")
patho_BrS_WBPCR      <- c("PNEU")
patho_BrS_NPCX_VT13  <- c("PNEU_VT13","PNEU_NOVT13")
patho_SS_BCX         <- c("HINF","MCAT","PNEU",
                          "SASP","SAUR",
                          "ENTRB","NFGNR","FUNGI","HAEMO","OTHSTR","NMEN","TB")
patho_SS_LA_ADJ      <-  cause_list[cause_list!="NoA"]
patho_SS_PF_ADJ      <-  cause_list[cause_list!="NoA"]
# bronze-standard measurements:
BrS_object_1 <- make_meas_object(patho_BrS_NPPCR,"NP","PCR","BrS",cause_list)
BrS_object_2 <- make_meas_object(patho_BrS_WBPCR,"WB","PCR","BrS",cause_list)
BrS_object_3 <- make_meas_object(patho_BrS_NPCX_VT13,"NP","CX","BrS",cause_list)
# silver-standard measurements:
SS_object_1  <- make_meas_object(patho_SS_BCX,"B","CX","SS",cause_list)
SS_object_2  <- make_meas_object(patho_SS_LA_ADJ ,"LA","_ADJ","SS",cause_list)
SS_object_3  <- make_meas_object(patho_SS_PF_ADJ ,"PF","_ADJ","SS",cause_list)

## parent directory for testing code (JHPCE):
#data_dir    <- "~/data/perch/"
#working_dir <- "~/PERCH/nplcm_reg_paper/"

# parent directory for testing code (LOCAL):
data_dir    <- "/users/qshi/PQ_data/"
working_dir <- "/users/qshi/PQ_results/"

# data_dir    <- "C:/Users/ppa Shi/Dropbox (PERCH)/PQ data/"
# working_dir <- "C:/Users/ppa Shi/Dropbox (PERCH)/PQ data/"

### Read Sensitivity


sens_data=read.csv(paste0(data_dir,"sensitivity.csv"))
sens_NPPCR_low_noabx=sens_data[sens_data$pathogen %in% patho_BrS_NPPCR,"NPPCR_low_noabx"]
sens_NPPCR_up_noabx=sens_data[sens_data$pathogen %in% patho_BrS_NPPCR,"NPPCR_up_noabx"]
sens_NPPCR_low_abx=sens_data[sens_data$pathogen %in% patho_BrS_NPPCR,"NPPCR_low_abx"]
sens_NPPCR_up_abx=sens_data[sens_data$pathogen %in% patho_BrS_NPPCR,"NPPCR_up_abx"]
sens_BCX_low_noabx=sens_data[sens_data$pathogen %in% patho_SS_BCX,"BCX_low_noabx"]
sens_BCX_up_noabx=sens_data[sens_data$pathogen %in% patho_SS_BCX,"BCX_up_noabx"]
sens_BCX_low_abx=sens_data[sens_data$pathogen %in% patho_SS_BCX,"BCX_low_abx"]
sens_BCX_up_abx=sens_data[sens_data$pathogen %in% patho_SS_BCX,"BCX_up_abx"]




## clean PERCH data:
clean_options <- list (raw_meas_dir       =  file.path(data_dir,"PQ_17JUN16.csv"),                 # <-------- where to look for raw data.
                       case_def           =  "CXRFINCAT_5",                                           # <------- case definition variable.
                       case_def_val       =  1,                                                    # <---------- case definition variable's value.
                       ctrl_def           =  "CASECONT",                                           # <------control definition variable.
                       ctrl_def_val       =  2,                                                    # <------ control definition variable's value.
                       X_strat            =  c("SITE8","HIV2"),                                            # <---- focus on the stratum defined by this variable.
                       X_strat_val        =  list("02GAM",0),                                       # <---- the stratum definition value.
                       BrS_objects        =  make_list(BrS_object_1,BrS_object_2),     # <---- all bronze-standard measurements.
                       SS_objects         =  make_list(SS_object_1),      # <----- all silver-standard measurements.
                       X_extra            = c("ENRLDATE","patid","AGECAT","AGE",
                                              "HIV2","PRABXBC","CASECONT_SUB",
                                              "CXRFIN_5","CXRFINCAT_5",
                                              "ALL_VS","ELCVSPN","SITE8","SEASONCAT","PRABXNP"),                      # <----- covariates besides case/control status.
                       patho_taxo_dir     = file.path(data_dir,"pathogen_category.csv"))        # <---- where to look for pathogen taxonomy information.
# you can add: "date_formats".

data_nplcm <- clean_perch_data(clean_options)

data_nplcm$X$std_date <- dm_Rdate_FPR(data_nplcm$X$ENRLDATE,data_nplcm$Y,effect = "fixed")

#
# visualize season:
#
# data_nplcm0 <-  data_nplcm
# ind_notNA   <- which(rowSums(is.na(data_nplcm0$Mobs$MBS$NPPCR))>0)
# data_nplcm_season <- data_nplcm0
# 
# if (length(ind_notNA)>0){
#   data_nplcm_season$Mobs$MBS$NPPCR <-  data_nplcm0$Mobs$MBS$NPPCR[-ind_notNA,,drop=FALSE]
#   data_nplcm_season$Mobs$MSS$BCX <-  data_nplcm0$Mobs$MSS$BCX[-ind_notNA,,drop=FALSE]
#   data_nplcm_season$X <-  data_nplcm0$X[-ind_notNA,]
#   data_nplcm_season$Y <-  data_nplcm0$Y[-ind_notNA]
# }
# 
# # explore covariate effect: focus on seasonality.
# ind_for_season_UI <- 1:length(patho_BrS_NPPCR) # pathogen indicators:
# 
# pdf(file.path(working_dir,"seasonlity.pdf"),width=8,height=8)
# par(mfrow=c(length(ind_for_season_UI),1))
# par(mar=c(1,5,0,5),oma = c(2,2,0,0),xpd=TRUE)
# for (i in ind_for_season_UI){
#   par(mar=c(1,5,3,5))
#   visualize_season(data_nplcm_season,i,1)
# }
# dev.off()

#
# Specify Model: 
#           

#no_abx_strat <- rep(1,sum(data_nplcm$Y)) # no TPR stratification; will use the first set of 'up' and 'down'.
abx_cases_NPPCR <- data_nplcm$X$PRABXNP[data_nplcm$Y==1]+1
abx_cases_BCX <- data_nplcm$X$PRABXBC[data_nplcm$Y==1]+1 # 2 for using abx, 1 otherwise; only for cases.

m_opt1 <- list(likelihood   = list(cause_list = cause_list,                # <---- fitted causes.
                                   k_subclass = c(1,1),                      # <---- no. of subclasses.
                                   Eti_formula = ~ as.factor(ALL_VS)+as.factor(AGECAT>1),                    # <---- etiology regression formula; only for cases.
                                   FPR_formula = list(
                                     NPPCR = ~  as.factor(AGECAT>1),
                                     WBPCR = ~  1
                                   )       # <----- false positive rate regression formula.
),      #  (placeholder; currently not used)<----- true positive rate regression formula. Informative priors for the coefficients are needed.
use_measurements = c("BrS","SS"),                           # <---- which measurements to use to inform etiology
prior        = list(Eti_prior   = overall_uniform(0.1, cause_list) ,                       # <--- etiology prior.
                    TPR_prior   = list(
                      BrS  = list(info  = "informative",
                                  input = "match_range",
                                  grp=abx_cases_NPPCR ,
                                  val   = list(
                                    NPPCR = list(up = list(sens_NPPCR_up_noabx,sens_NPPCR_up_abx),
                                                 low = list(sens_NPPCR_low_noabx,sens_NPPCR_low_abx)
                                    ),
                                    WBPCR = list(up = list(0.99,0.5),
                                                 low = list(0.5,0.25)
                                    )
                                    #NPCX  = list(up = rep(0.99,length(BrS_object_3$patho)),
                                    #             low = rep(0.5,length(BrS_object_3$patho))
                                    #)
                                  )
                      ),
                      SS   = list(info  = "informative",
                                  input = "match_range",
                                  grp=abx_cases_BCX,
                                  val   = list(
                                    BCX = list(up = list(sens_BCX_up_noabx,sens_BCX_up_abx),# last one is 0.3 for SAF and 0.2 for others
                                               low = list(sens_BCX_low_noabx,sens_BCX_low_abx)
                                               
                                    )#, # <--- BCX TPR prior
                                    #LA_ADJ = list(list(up = rep(0.15,length(SS_object_2$patho)),
                                    #                   low = rep(0.05,length(SS_object_2$patho)))
                                    #), # <--- LA_ADJ TPR prior
                                    #PF_ADJ = list(list(up = rep(0.15,length(SS_object_3$patho)),
                                    #                   low = rep(0.05,length(SS_object_3$patho)))
                                    #) # <--- PF_ADJ TPR prior
                                  )
                      )
                    )  # <---- TPR prior.
)
)     

model_options <- m_opt1

## Kenya; CXR+ cases; age and seasonality regression:
#model_options$likelihood$Eti_formula      <- ~ -1+as.factor(AGE)+as.factor(ALL_VS)+as.factor(SEASONCAT)
#model_options$likelihood$FPR_formula[[1]] <- ~ -1+as.factor(AGE)+as.factor(ALL_VS)+as.factor(SEASONCAT)
#model_options$likelihood$FPR_formula[[2]] <- ~ -1+as.factor(AGE)+as.factor(ALL_VS)+as.factor(SEASONCAT)

model_options$likelihood$Eti_formula      <- ~ -1+as.factor(AGE)
model_options$likelihood$FPR_formula[[1]] <- ~ -1+as.factor(AGE)
model_options$likelihood$FPR_formula[[2]] <- ~ -1+as.factor(AGE)

## Kenya; CXR+ cases; age and seasonality regression:
#model_options$likelihood$Eti_formula      <- ~ -1+as.factor(ALL_VS)+as.factor(AGE)+as.factor(SEASONCAT)
#model_options$likelihood$FPR_formula[[1]] <- ~ -1+as.factor(ALL_VS)+as.factor(AGE)+as.factor(SEASONCAT)
#model_options$likelihood$FPR_formula[[2]] <- ~ -1+as.factor(ALL_VS)+as.factor(AGE)+as.factor(SEASONCAT)

Z_Eti00       <- stats::model.matrix(model_options$likelihood$Eti_formula,
                                   data.frame(data_nplcm$X,
                                              data_nplcm$Y)[data_nplcm$Y==1,,drop=FALSE])

str(Z_Eti00)
#tmp <- model.matrix(model_options$likelihood$Eti_formula,data = data_nplcm$X)
#tmp <- model.matrix(model_options$likelihood$FPR_formula[[1]],data = data_nplcm$X)

# Z_Eti0       <- stats::model.matrix(model_options$likelihood$Eti_formula,
#                                     data.frame(data_nplcm$X,data_nplcm$Y)[data_nplcm$Y==1,,drop=FALSE])
# a.eig        <- eigen(t(Z_Eti0)%*%Z_Eti0)
# sqrt_Z_Eti0  <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
# 
# Z_Eti  <- Z_Eti0%*%solve(sqrt_Z_Eti0)
# #Z_Eti  <- scale(Z_Eti)

assign_model(model_options,data_nplcm)

#
# 4. Fit model:
#
add.name="_A_New_5k"
# date stamp for analysis:
Date     <- gsub("-", "_", Sys.Date())
# include stratification information in file name:
dated_strat_name    <- paste(working_dir,Date,paste(c(clean_options$case_def,clean_options$case_def_val,clean_options$X_strat_val),collapse="_"),add.name,collapse="_")
# create folder
dir.create(dated_strat_name)
fullname <- dated_strat_name

# for finer scenarios, e.g., different types of analysis applicable to the
# same data set. Here we just perform one analysis:
result_folder <- fullname
dir.create(result_folder)

t1 <- Sys.time()

# options for MCMC chains:
mcmc_options <- list(debugstatus = TRUE,
                     n.chains   = 1,
                     n.itermcmc = 5000,
                     n.burnin   = 1000,
                     n.thin     = 5,
                     individual.pred = TRUE,
                     ppd             = !TRUE,
                     get.pEti        = TRUE,
                     result.folder = result_folder,
                     bugsmodel.dir = result_folder,
                     jags.dir = "",
                     use_jags = TRUE)

# Record the settings of current analysis:
#data clean options:
dput(clean_options,file.path(mcmc_options$result.folder,"data_clean_options.txt"))
dput(data_nplcm,file.path(mcmc_options$result.folder,"data_nplcm.txt")) # <-- in case covariate data X are needed for plotting.

rjags::load.module("glm")

gs <- nplcm(data_nplcm,model_options,mcmc_options)
#source("~/downloads/test_baker/regression_primary/inner_nplcm.R")

# Writing raw data into the folder.
dput(clean_options,file.path(mcmc_options$result.folder,"data_clean_options.txt"))
write.csv(data_nplcm$X,file=paste0(result_folder,"/rawdata_x.csv")) ## Save Raw Data
Y=data_nplcm$Y
mobs=cbind(do.call("cbind",data_nplcm$Mobs$MBS),do.call("cbind",data_nplcm$Mobs$MSS))
write.csv(mobs[Y==1,],file=paste0(result_folder,"/rawdata_mobs.csv"))


source(paste0(data_dir,"Util.R"))
## Traceplot
pqtrace(result_folder,param="pEti",path.name=FALSE,filename="traceplot",output.folder=result_folder)

