# This script contains many tests regarding manual computation of IV analysis. 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table", 
                   "foreign", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "fixest", "sandwich", "lmtest", "boot",
                   "ggplot2") 
# "pglm", "multiwayvcov", "clusterSEs", "alpaca", "clubSandwich",

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Load them 
lapply(neededPackages, library, character.only = TRUE)

# /!\/!\ IF renv::restore(neededPackages) FAILS TO INSTALL SOME PACKAGES /!\/!\ 

# For instance sf could cause trouble https://github.com/r-spatial/sf/issues/921 
# or magick, as a dependency of raster and rgdal. 

# FOLLOW THESE STEPS:
# 1. Remove these package names from neededPackages above, and rerun renv::restore(packages = neededPackages)
# 2. Write them in troublePackages below, uncomment, and run the following code chunk: 

# # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c() 
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/reg_results")



### PARCEL SIZE
parcel_size <- 3000

### SET NUMBER OF THREADS USED BY {FIXEST} TO ONE (TO REDUCE R SESSION CRASH)
getFixest_nthreads()




catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
island <- "Kalimantan"
outcome_variable <- "lucpfip_pixelcount_total"
dynamics <- TRUE
commo <- c("ffb", "cpo")
yoyg <- FALSE
short_run <- "unt level" # from c("unt_level", "dev", "yoyg")
imp <- 1
x_pya <- 3 # 2, 3 or 4
lag_or_not <- "_lag1" # from c("_lag1", "")
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
pya_ov <- FALSE
oneway_cluster <- ~parcel_id


# DATA
d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                              parcel_size/1000,"km_",
                              catchment_radius/1000,"CR.rds")))


# subsampling
d <- d[d$island %in% island,]

if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
  d <- d[d$any_pfc2000_total,]
}
if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
   outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
   outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
  d <- d[d$any_fc2000_30th,]
}



### Specifications
if(dynamics == FALSE){ 
  # Only the overall effect is investigated then, no short vs. long run
  # In this case, the variables have name element _Xya_ with X the number of years over which the mean has been 
  # computed, always including the contemporaneous record. See add_parcel_variables.R
  #short_run <- ""
  
  # if we omit one commodity 
  if(length(commo) == 1){
    if(yoyg == TRUE){
      regressors <- paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not)
    }else{
      regressors <- paste0("wa_",commo,"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)
    }
  }
  # if we don't omit a commodity. 
  if(length(commo) == 2){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not))
      
    }else{
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)) 
    }
  }
}

# if, on the other hand, we want to disentangle the short run effects from the long run's. 
if(dynamics == TRUE){ 
  
  if(length(commo) == 1){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                      paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
    }
    if(yoyg == FALSE){
      if(short_run == "unt level"){
        regressors <- c(paste0("wa_",commo,"_price_imp",imp,lag_or_not), # SR measure
                        paste0("wa_",commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
      }
      if(short_run == "dev"){
        regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not), # SR measure
                        paste0("wa_",commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
      }  
    }
  }
  
  if(length(commo) == 2){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                      paste0("wa_",commo[1],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not))    
    }
    if(yoyg == FALSE){
      if(short_run == "unt level"){
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,lag_or_not),# FFB SR measure
                        paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not), # FFB LR measure 
                        paste0("wa_",commo[2],"_price_imp",imp,lag_or_not),# CPO SR measure
                        paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # CPO LR measure 
      }
      if(short_run == "dev"){
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not))
      }
    }
  }
}

# add pya outcome variable 
if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}

# Formula
fe <- "parcel_id"# + district^year

fml <- as.formula(paste0(outcome_variable,
                         " ~ ",
                         paste0(regressors, collapse = "+"),
                         " + ",
                         paste0(paste0(controls,lag_or_not), collapse = "+"),
                         " | ",
                         fe))

d$fact_parcel_id <- as.factor(d$parcel_id)
fml_mfe <- as.formula(paste0(outcome_variable,
                             " ~ ",
                             paste0(regressors, collapse = "+"),
                             " + ",
                             paste0(paste0(controls,lag_or_not), collapse = "+"),
                             " + fact_parcel_id"))

fml_std <-  as.formula(paste0(outcome_variable,
                              " ~ ",
                              paste0(regressors, collapse = "+"),
                              " + ",
                              paste0(paste0(controls,lag_or_not), collapse = "+")))

# MANUAL 2SLS
# First, regress first stage using feols. Nothing special here. SE method does not matter, 
# the only purpose is computing fitted values
# Second, regress second stage with feglm with fitted values instead of price signal endogenous variable.
# Coefficients are correct, but not SE. Save coefficients only
# Run again feglm but on endogenous variable. Save residuals
# plug residuals in  

### 2SLS BY HAND ### 
# fit several models, using different packages
# with fixed effects
est_pois <- fixest::feglm(fml, data = d,
                          family = "poisson")
library(alpaca)
alp_pois <- alpaca::feglm(fml, data = d, 
                          family = poisson(link = "log"))
library(pglm)
pglm_pois <- pglm(fml_std, 
                  data = d,
                  index = c("parcel_id", "year"),
                  family = "poisson", 
                  effect = "individual", 
                  model = "within")
# (fixest, alp and pglm all yield the same coeff)
round(est_pois$coefficients,8) == round(alp_pois$coefficients,8)
round(alp_pois$coefficients, 8) == round(pglm_pois$estimate, 8)

# without fixed effects
glm_pois <- glm(fml_std,
                family = "poisson",
                data = d)

# packages that compute standard errors do not accept fitted model objects from 
# the functions we need to use to run glm-FE regressions 
# alpaca and fixest not at all, pglm is accepted by sandwich but then cluster does not work. 

# try each SE package with each regression model output. 

# sandwich
sdw_vcov_nocl <- sandwich::vcovCL(glm_pois)
sdw_vcov <- sandwich::vcovCL(glm_pois, cluster = ~parcel_id)

sdw_vcov <- sandwich::vcovCL(pglm_pois)
sdw_vcov <- sandwich::vcovCL(pglm_pois, cluster = ~parcel_id) # probably fails because number of clusters specified does not 
# correspond to numbers of grid cells actually used in the regression. 

sdw_vcov <- sandwich::vcovCL(est_pois)

sdw_vcov <- sandwich::vcovCL(alp_pois)


# multiwayvcov
library(multiwayvcov)
mw_vcov <- multiwayvcov::cluster.vcov(glm_pois, cluster = ~parcel_id)
all.equal(sdw_vcov, mw_vcov) #not very precisely equal... anyway

mw_vcov <- multiwayvcov::cluster.vcov(pglm_pois, cluster = ~parcel_id, force_posdef = T)

mw_vcov <- multiwayvcov::cluster.vcov(est_pois, cluster = ~parcel_id, force_posdef = T)

mw_vcov <- multiwayvcov::cluster.vcov(alp_pois, cluster = ~parcel_id)

# clusterSEs - nothing works
library(clusterSEs)
# does not work bc of the intercept
clSE_vcov <- clusterSEs::cluster.im.glm(glm_pois, 
                                        dat = d, cluster = ~parcel_id, return.vcv = TRUE)

clSE_vcov <- clusterSEs::cluster.im.glm(pglm_pois, 
                                        dat = d, cluster = ~parcel_id, return.vcv = TRUE, drop = TRUE)

clSE_vcov <- clusterSEs::cluster.im.glm(est_pois, 
                                        dat = d, cluster = ~parcel_id, return.vcv = TRUE, drop = TRUE)

clSE_vcov <- clusterSEs::cluster.im.glm(alp_pois, 
                                        dat = d, cluster = ~parcel_id, return.vcv = TRUE, drop = TRUE)


# clubSandwich
library(clubSandwich)
csw_vcov <- clubSandwich::vcovCR(glm_pois, cluster = d$parcel_id, type = "CR1")

csw_vcov <- clubSandwich::vcovCR(pglm_pois, cluster = d$parcel_id, type ="CR1")

csw_vcov <- clubSandwich::vcovCR(est_pois, cluster = d$parcel_id, type ="CR1")

csw_vcov <- clubSandwich::vcovCR(alp_pois, cluster = d$parcel_id, type ="CR1")

### ### ### ###
# many obs. will be removed because of NA values or because of only 0 outcome. 
# but how fixest removes NAs is not clear. 
d_nona <- dplyr::filter(d, !is.na(wa_ffb_price_imp1_lag1) & 
                          !is.na(wa_ffb_price_imp1_3pya_lag1) & 
                          !is.na(wa_cpo_price_imp1_3pya_lag1) & 
                          !is.na(wa_cpo_price_imp1_lag1) & 
                          !is.na(wa_pct_own_loc_gov_imp) & 
                          !is.na(wa_pct_own_nat_priv_imp) & 
                          !is.na(wa_pct_own_for_imp) & 
                          !is.na(n_reachable_uml))

d_nona <- d_nona[,c("parcel_id", "year", outcome_variable, regressors, controls)]
nrow(d) - nrow(d_nona) # number of rows that were removed because at least one variable is NA

est_pois <- fixest::feglm(fml, data = d,
                          family = "poisson")
# automatic package demeaning
auto_dm <- demean(d_nona[,c(outcome_variable, regressors, controls)], fe= d_nona[,"parcel_id"])
d_nona[, paste0(c(outcome_variable, regressors, controls),"_dm")] <- auto_dm
d_nona[1:30,c("parcel_id", "year", "lucpfip_pixelcount_total", "lucpfip_pixelcount_total_dm", 
              "wa_ffb_price_imp1_lag1", "wa_ffb_price_imp1_lag1_dm")]

# manual demeaning using all information, for say parcel _id 22336
# once deprived of all rows with at least one var missing, this parcel has 
d_nona[d_nona$parcel_id==22336,"wa_ffb_price_imp1_lag1"] %>% length() # 3 observations
# and a mean of 
auto_m <- d_nona[d_nona$parcel_id==22336,"wa_ffb_price_imp1_lag1"] %>% mean()
# while 
sum(!is.na(d[d$parcel_id==22336,"wa_ffb_price_imp1_lag1"])) # 8 observations are not missing for this variable, 
# and could therefore be used to demean the data... 
m <- mean(d[d$parcel_id==22336,"wa_ffb_price_imp1_lag1"],na.rm = TRUE)
m == auto_m
d_nona[d_nona$parcel_id==22336,"wa_ffb_price_imp1_lag1"] - m # is different than auto_dm 

fml_dm <-  as.formula(paste0(paste0(outcome_variable,"_dm"),
                             " ~ ",
                             paste0(paste0(regressors,"_dm"), collapse = "+"),
                             " + ",
                             paste0(paste0(controls,lag_or_not, "_dm"), collapse = "+")))

est_pois_dm <- fixest::feglm(fml_dm, data = d_nona, family = "poisson")

d_nona$lucpfip_pixelcount_total_dm
## SO WE SHOULD DEMEAN OUR DATA OURSELVES 


nrow(auto_dm)
auto_dm[1:15,]
# so fixest removes obs. as long as *ANY* of the variables is missing. And this is also the case 
# in the demeaning process, meaning that for such obs., it doesn't use the information on other 
# variables than the one missing. 

d1 <- d[-obs2remove(fml, d, "poisson"), c("parcel_id", "year", outcome_variable, regressors, controls)] 
# removes 78180 obs. with only 0 outcomes (and they might also have NA values)

d1[d1$parcel_id==225,,drop = FALSE] %>% nrow()




# Try with fixest only
obs2remove(fml, d, "poisson") %>% summary()
d_dm1 <- demean(d[,c(outcome_variable, regressors, controls)], 
                fe= d[,"parcel_id"])
d_dm <- demean(d[-obs2remove(fml, d, "poisson"),c(outcome_variable, regressors, controls)], 
               fe= d[-obs2remove(fml, d, "poisson"),"parcel_id"])

est_pois <- fixest::feglm(fml, data = d,
                          family = "poisson")

class(d_dm)
(d_dm[,outcome_variable] <= 0) %>% 
  d[,c(paste0(outcome_variable,"_dm"), paste0(regressors,"_dm"), paste0(controls,"_dm"))] <- d_dm

feglm()
#theoretical first stage 
# prepare data - keep only those that fixest will keep
# those with no NA
d_nona <- dplyr::filter(d, !is.na(wa_ffb_price_imp1_lag1) & 
                          !is.na(wa_ffb_price_imp1_3pya_lag1) & 
                          !is.na(wa_cpo_price_imp1_3pya_lag1) & 
                          !is.na(wa_cpo_price_imp1_lag1) & 
                          !is.na(wa_pct_own_loc_gov_imp) & 
                          !is.na(wa_pct_own_nat_priv_imp) & 
                          !is.na(wa_pct_own_for_imp) & 
                          !is.na(n_reachable_uml) & 
                          !is.na(wa_prex_cpo_imp1_lag1))
nrow(d) - nrow(d_nona) # number of rows that were removed because at least one variable is NA
# and those with not only zero outcome
names(d_nona)
d_clean <- d_nona[-obs2remove(lucpfip_pixelcount_total~parcel_id + year, d_nona, family = "poisson"),]
nrow(d_clean)

fststg_fml <- as.formula(paste0(regressors[1],
                                " ~ ",
                                paste0(regressors[-1], collapse = "+"),
                                " + wa_prex_cpo_imp1_lag1 +",
                                paste0(paste0(controls,lag_or_not), collapse = "+"),
                                " | parcel_id + year"))

est_fststg <- fixest::feols(fststg_fml, data = d_clean)

d_clean$fit_endo <- est_fststg$fitted.values

# theorietical second stage 
# positions in the formula matter, for matrix product later
sndstg_fml <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                "fit_endo +",
                                paste0(regressors[-1], collapse = "+"),
                                " + ",
                                paste0(paste0(controls,lag_or_not), collapse = "+"),
                                " | parcel_id + year"))

est_sndstg <- fixest::feglm(sndstg_fml, data = d_clean,
                            family = "poisson")
su_est_sndstg <- summary(est_sndstg, se="standard")

# and endogeneous model 
endo_fml <- as.formula(paste0(outcome_variable,
                              " ~ ",
                              paste0(regressors, collapse = "+"),
                              " + ",
                              paste0(paste0(controls,lag_or_not), collapse = "+"),
                              " | parcel_id + year"))

# compute correct residuals manually
cor_res <- as.matrix(d_clean[,c(regressors, controls)])%*%est_sndstg$coefficients
cor_res <- as.numeric(cor_res)
cor_est_sndstg <- est_sndstg
cor_est_sndstg$residuals <- cor_res
est_sndstg$residuals %>% class()
cor_est_sndstg$residuals %>% class()
# only rediuals are changed (no automatic recomputation of vcov arguments, normal)
all.equal(est_sndstg$cov.unscaled,cor_est_sndstg$cov.unscaled)
all.equal(vcov(est_sndstg), vcov(cor_est_sndstg))
all.equal(bread(est_sndstg),bread(cor_est_sndstg))
all.equal(est_sndstg$residuals, cor_est_sndstg$residuals)
head(est_sndstg$residuals)
head(cor_est_sndstg$residuals)
# and no recomputation in calls of summary either
cor_su_est_sndstg <- summary(cor_est_sndstg, se="cluster")
etable(su_est_sndstg, cor_su_est_sndstg)

# try removing cov.unscaled args --> does NOT prevent SE computation
rm(cor_su_est_sndstg)
cor_est_sndstg <- est_sndstg
cor_est_sndstg$residuals <- cor_res
cor_est_sndstg$cov.unscaled <- NULL
cor_su_est_sndstg <- summary(cor_est_sndstg, se="cluster")
etable(su_est_sndstg, cor_su_est_sndstg)
# try replacing cov.unscaled args --> CHANGES SE computation
rm(cor_su_est_sndstg)
cor_est_sndstg <- est_sndstg
cor_res <- as.matrix(d_clean[,c(regressors, controls)])%*%est_sndstg$coefficients
cor_res <- as.numeric(cor_res)
cor_est_sndstg$residuals <- cor_res
rep_mat <- matrix(data = rnorm(length(su_est_sndstg$coefficients)),
                  nrow = length(su_est_sndstg$coefficients), 
                  ncol = length(su_est_sndstg$coefficients))
colnames(rep_mat) <- attr(cor_est_sndstg$cov.unscaled, "dimnames")[[1]]
rownames(rep_mat) <- attr(cor_est_sndstg$cov.unscaled, "dimnames")[[1]]

cor_est_sndstg$cov.unscaled <- rep_mat

cor_su_est_sndstg <- summary(cor_est_sndstg, se="standard")
etable(su_est_sndstg, cor_su_est_sndstg)

# maintenant on calcule notre vcov unscaled basé sur nos correct residuals. 
# on pourra vérifier qu'on a pas fait de la merde si on peut reproduire la cov.unscaled produite par la fonction 
vcov <- cov(as.matrix(d_clean[,c("fit_endo", regressors[-1], controls)]))
cov(est_sndstg$coefficients, est_sndstg$coefficients)
round(est_sndstg$cov.unscaled, 5) == round(solve(est_sndstg$hessian), 5)
# we need to pass from the Hessian. 
# compute it from the model with the endogeneous variable
rm(cor_su_est_sndstg)
cor_est_sndstg <- est_sndstg
cor_res <- as.matrix(d_clean[,c(regressors, controls)])%*%est_sndstg$coefficients
cor_res <- as.numeric(cor_res)
cor_est_sndstg$residuals <- cor_res
est_endo <- fixest::feglm(endo_fml, data = d_clean,
                          family = "poisson")

new_cov.unscaled <- solve(est_endo$hessian)
cor_est_sndstg$cov.unscaled <- new_cov.unscaled
cor_su_est_sndstg <- summary(cor_est_sndstg, se="cluster")
etable(su_est_sndstg, cor_su_est_sndstg)

# Now verify that our "hessian-unscaled vcov" method works. 
# i.e. verify it produces the same standard errors as a 2sls facility, in a non-glm case, and in a non fe case. 

# no FE, linear, test
# automatic se calculation
library(AER)
simple_aer_fml <- as.formula(paste0(outcome_variable,
                                    " ~ ",
                                    paste0(regressors, collapse = "+"),
                                    " + ",
                                    paste0(paste0(controls,lag_or_not), collapse = "+"),
                                    " | wa_prex_cpo_imp1_lag1 +",
                                    paste0(regressors[-1], collapse = "+"),
                                    " + ",
                                    paste0(paste0(controls,lag_or_not), collapse = "+")))

est_simple_aer <- ivreg(simple_aer_fml, data = d_clean)
summary(est_simple_aer)

# manual se calculation

fststg_fml <- as.formula(paste0(regressors[1],
                                " ~ wa_prex_cpo_imp1_lag1 +",
                                paste0(regressors[-1], collapse = "+"),
                                " + ",
                                paste0(paste0(controls,lag_or_not), collapse = "+")))

est_fststg <- fixest::feglm(fststg_fml, data = d_clean, family = "gaussian")

d_clean$fit_endo <- est_fststg$fitted.values

# theorietical second stage 
# positions in the formula matter, for matrix product later
sndstg_fml <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                "fit_endo +",
                                paste0(regressors[-1], collapse = "+"),
                                " + ",
                                paste0(paste0(controls,lag_or_not), collapse = "+")))

est_sndstg <- fixest::feglm(sndstg_fml, data = d_clean, family = "gaussian")


su_est_sndstg <- summary(est_sndstg, se="standard")

# and endogeneous model 
endo_fml <- as.formula(paste0(outcome_variable,
                              " ~ ",
                              paste0(regressors, collapse = "+"),
                              " + ",
                              paste0(paste0(controls,lag_or_not), collapse = "+")))

cor_est_sndstg <- est_sndstg

# residuals
cor_fit <- as.matrix(cbind(rep(1, nrow(d_clean)),
                           d_clean[,c(regressors, controls)]))%*%est_sndstg$coefficients
cor_res <- as.matrix(d_clean[,outcome_variable]) - cor_fit
d_clean[14,outcome_variable] - cor_fit[14,] == cor_res[14,]

cor_res <- as.numeric(cor_res)
cor_est_sndstg$residuals <- cor_res

est_endo <- fixest::feglm(endo_fml, data = d_clean, family = "gaussian")

new_cov.unscaled <- solve(est_endo$hessian)
cor_est_sndstg$cov.unscaled <- new_cov.unscaled
cor_su_est_sndstg <- summary(cor_est_sndstg, se = "standard")
#etable(su_est_sndstg, cor_su_est_sndstg)

summary(cor_est_sndstg, se = "standard", dof = dof(
  adj = TRUE,
  fixef.K = "nested",
  cluster.adj = FALSE#,
  # cluster.df = "conventional",
  # t.df = "conventional",
  # fixef.force_exact = FALSE
))

summary(est_simple_aer, vcov = sandwich)
summary(est_sndstg, se = "white")
vcov(est_simple_aer)
vcov(cor_su_est_sndstg)
# well it does not work. (and note that the manual method yields different se if we use feglm gaussian or feols)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# so taking the Hessian of the endogenous model did not work. 
# now try to compute the correct vcov. 
# recall that replacing cov.unscaled args --> CHANGES SE computation
# do it with simpler formulas : no "regressor" object, just one endo variable and the controls


# data first 
d_nona <- dplyr::filter(d, !is.na(wa_ffb_price_imp1_lag1) &
                          !is.na(n_reachable_uml_lag1) & 
                          !is.na(wa_prex_cpo_imp1_lag1))
nrow(d) - nrow(d_nona) # number of rows that were removed because at least one variable is NA
# and those with not only zero outcome
names(d_nona)
d_clean <- d_nona[-obs2remove(lucpfip_pixelcount_total~parcel_id + year, d_nona, family = "poisson"),]
nrow(d_clean)

outcome_variable 
endo_var <- "wa_ffb_price_imp1_lag1"
exo_var <- "n_reachable_uml_lag1"
instr_var <- "wa_prex_cpo_imp1_lag1"
d_clean <- dplyr::select(d_clean, parcel_id, year, outcome_variable, endo_var, exo_var, instr_var)
dim(d_clean)
names(d_clean)
write.dta(d_clean, "temp_data/temp_sample.dta")
# theoretical first stage
fststg_fml <- as.formula(paste0(endo_var," ~ ",
                                exo_var, " + ", instr_var))

est_fststg <- fixest::feols(fststg_fml, data = d_clean)

d_clean$fit_endo <- est_fststg$fitted.values

# theorietical second stage 
# positions in the formula matter, for matrix product later
sndstg_fml <- as.formula(paste0(outcome_variable," ~ ",
                                "fit_endo +", exo_var))

est_sndstg <- fixest::feols(sndstg_fml, data = d_clean)


# compute the correct residuals: the diff between the true outcome and the fitted model with the 
# true endogenous variable and not the fitted endogenous variable. 
# here we add a 1 column to the rhs variables to account for the presence of the intercept. 
# with fixed effects we will have to use predict? 
# use of predict requires that the true endo is named after the fitted endo exactly.
predi_data <- dplyr::select(d_clean, -fit_endo)
names(predi_data)[names(predi_data)==endo_var] <- "fit_endo"

cor_fit <- predict(est_sndstg, 
                   newdata = predi_data, 
                   type = "response")

cor_res <- d_clean[,outcome_variable] - cor_fit
cor_res_sq <- cor_res^2
d_clean[,outcome_variable] %>% head()
cor_fit %>% head()
cor_res %>% head()
cor_res_sq %>% head()
dfk <- est_sndstg$nobs - est_sndstg$nparams
cor_rmse <- sqrt(sum(cor_res_sq)/dfk)

manual_rmse <- sqrt(sum((est_sndstg$residuals)^2)/dfk)
manual_vce <- est_sndstg$cov.unscaled 
class(manual_vce)
cor_vce <- ((cor_rmse/manual_rmse)^2)*manual_vce
class(cor_vce)

cor_est_sndstg <- est_sndstg
cor_est_sndstg$cov.unscaled <- cor_vce
#etable(est_sndstg, cor_est_sndstg)


# and with the automatic iv 
library(AER)
simple_aer_fml <- as.formula(paste0(outcome_variable, " ~ ", 
                                    endo_var, " + ", exo_var, 
                                    " | ", instr_var, " + ", exo_var))

est_simple_aer <- ivreg(simple_aer_fml, data = d_clean)
summary(est_simple_aer)
summary(cor_est_sndstg)

# SE ARE THE SAME ! 

# now try to rescale the vcov to compute white SE. 
summary(cor_est_sndstg, se = "white")
# and clustered
summary(cor_est_sndstg, cluster = ~parcel_id)

# they both yield LOWER SE... (in stata too, with ivreghdfe with option cluster or not)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###  ### ### ### ### ### ###

# redo simply the same thing, on the manual side, but with fixed effects. Compare with ivreghdfe in Stata
fe <- "parcel_id"

# data first 
d_nona <- dplyr::filter(d, !is.na(wa_ffb_price_imp1_lag1) &
                          !is.na(n_reachable_uml_lag1) & 
                          !is.na(wa_prex_cpo_imp1_lag1))
nrow(d) - nrow(d_nona) # number of rows that were removed because at least one variable is NA
# and those with not only zero outcome
names(d_nona)
d_clean <- d_nona[-obs2remove(lucpfip_pixelcount_total~parcel_id + year, d_nona, family = "poisson"),]
nrow(d_clean)

outcome_variable 
endo_var <- "wa_ffb_price_imp1_lag1"
exo_var <- "n_reachable_uml_lag1"
instr_var <- "wa_prex_cpo_imp1_lag1"
d_clean <- dplyr::select(d_clean, parcel_id, year, outcome_variable, endo_var, exo_var, instr_var)
dim(d_clean)
names(d_clean)
write.dta(d_clean, "temp_data/temp_sample.dta")
# theoretical first stage
fststg_fml <- as.formula(paste0(endo_var," ~ ",
                                exo_var, " + ", instr_var,
                                " | ", fe))

est_fststg <- fixest::feols(fststg_fml, data = d_clean)

d_clean$fit_endo <- est_fststg$fitted.values

# theorietical second stage 
# positions in the formula matter, for matrix product later
sndstg_fml <- as.formula(paste0(outcome_variable," ~ ",
                                "fit_endo +", exo_var,
                                " | ", fe))

est_sndstg <- fixest::feols(sndstg_fml, data = d_clean)

# compute the correct residuals: the diff between the true outcome and the fitted model with the 
# true endogenous variable and not the fitted endogenous variable. 
# here we add a 1 columne to the rhs variables to account for the presence of the intercept. 
# with fixed effects we will have to use predict? 
# use of predict requires that the true endo is named after the fitted endo exactly.
predi_data <- dplyr::select(d_clean, -fit_endo)
names(predi_data)[names(predi_data)==endo_var] <- "fit_endo"

cor_fit <- predict(est_sndstg, 
                   newdata = predi_data, 
                   type = "response")

cor_res <- d_clean[,outcome_variable] - cor_fit
cor_res_sq <- cor_res^2
d_clean[,outcome_variable] %>% head()
cor_fit %>% head()
cor_res %>% head()
cor_res_sq %>% head()
dfk <- est_sndstg$nobs - est_sndstg$nparams

cor_rmse <- sqrt(sum(cor_res_sq)/dfk)

manual_rmse <- sqrt(sum((est_sndstg$residuals)^2)/dfk)
manual_vce <- est_sndstg$cov.unscaled 
class(manual_vce)
cor_vce <- ((cor_rmse/manual_rmse)^2)*manual_vce
class(cor_vce)

cor_est_sndstg <- est_sndstg
cor_est_sndstg$cov.unscaled <- cor_vce
summary(cor_est_sndstg, se = "cluster")



# degrees of freedom 
dof_rep <- dof(t.df= "min")
dof()
summary(cor_est_sndstg, se = "standard", dof = dof(adj = 0, 
                                                   fixef.K = "nested",
                                                   fixef.exact = TRUE, 
                                                   cluster.adj = FALSE))#

etable(est_sndstg, cor_est_sndstg, se = "cluster")




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Now design the CONTROL FUNCTION APPROACH WITH BOOTSTRAP
ctrl_fun_endo <- function(myfun_data, indices, fixed_effects){
  
  d_clean_sampled <- d_clean[indices,]
  fststg_fml <- as.formula(paste0(endo_var," ~ ",
                                  exo_var, " + ", instr_var,
                                  " | ", fixed_effects))
  
  est_fststg <- fixest::feols(fststg_fml, data = d_clean_sampled)
  
  d_clean_sampled$fst_res <- est_fststg$residuals
  head(d_clean$fst_res)
  # residuals are the same as in Stata 
  # reghdfe wa_ffb_price_imp1_lag1 wa_prex_cpo_imp1_lag1 n_reachable_uml_lag1, absorb(parcel_id) residuals(v2h_fe)
  
  # theorietical second stage 
  # positions in the formula matter, for matrix product later
  sndstg_fml <- as.formula(paste0(outcome_variable," ~ ",
                                  endo_var, "+ fst_res +", exo_var,
                                  " | ", fixed_effects))
  
  est_sndstg <- fixest::feglm(sndstg_fml, data = d_clean_sampled, family = "poisson")
  #su_est <- summary(est_sndstg, se = "white")
  return(est_sndstg$coefficients)
}



### Design the bootstrap sampling function 

original_data <- d_clean
cluster_var <- "parcel_id"
# get the different cluster sizeS
sizes <- sapply(cluster_names, FUN = function(name){
                                        sum(as.character(d_clean[,cluster_var]) == name)
                                       })
# this is a vector of same length as cluster_names. Values are the number of observations (the size) per cluster. 

ran.gen_cluster <- function(original_data, cluster_var, sizes){
  cluster_names <- original_data[,cluster_var] %>% unique() %>% as.character()
 
  # sample_cl %>% length()
  # sample_cl %>% unique() %>% length()
  # a_duplic_cl <- sample_cl[duplicated(sample_cl)][1]
  # bb <- table(sample_cl)
  # bb[a_duplic_cl]
  # g_name <- as.integer(a_duplic_cl)
  
  cl_boot_dat <- NULL
  for(s in unique(sizes)){
    
    # names and numbers of clusters of size s
    cluster_names_s <- names(sizes[sizes == s])
    G_s <- length(cluster_names_s)
    
    # resample cluster names 
    sample_cl_s <- sample(cluster_names_s, G_s, replace = TRUE) 
    
    for(g in 1:G_s){
      # name of the g_th cluster that got resampled
      g_name <- sample_cl_s[g]
      
      # select the s observations from cluster g 
      # if this cluster name was picked more than once, then its obs. will be added again for another value of g. 
      cc <- original_data[as.character(original_data[,cluster_var]) == g_name,]
      
      # we need to give a new cluster identifier during the resampling, otherwise a cluster sampled more than once 
      # will be "incorrectly treated as one large cluster rather than two distinct cluster" (by the fixed effects) (Cameron 2015)
      cc[,cluster_var] <- paste0(s,"_",g)
      
      cl_boot_dat <- rbind(cl_boot_dat, cc)
    }
  }
  # test that the returned data are the same dimension as input
  # dim(cl_boot_dat) 
  # dim(original_data)
  # dim(cl_boot_dat) == dim(d_clean)
  # 
  # # test new clusters are not duplicated
  # anyDuplicated(cl_boot_dat[,c(cluster_var,"year")])
  # 
  
  # design in apply  
  # needs two lower level functions imbricated
  # one that samples 
  # one, lower, that extracts individual clusters 
  # extract_cl_d <- function(g_name){
  #          # for(g in 1:arg_list[["number_clusters"]][[s]]){
  #               # name of the g_th cluster that got resampled
  #               # g_name <- sample_cl_s[g]
  #               
  #               # select the s observations from cluster g 
  #               # if this cluster name was picked more than once, then its obs. will be added again for another value of g. 
  #               cc <- original_data[as.character(original_data[,arg_list[["cluster_variable"]]]) == g_name,]
  #               # we need to give a new cluster identifier during the resampling, otherwise a cluster sampled more than once 
  #               # will be "incorrectly treated as one large cluster rather than two distinct cluster" (by the fixed effects) (Cameron 2015)
  #               cc[,arg_list[["cluster_variable"]]] <- paste0(s,"_",g)
  #               
  #               return(cc)
  # }
  # bd_s <- lapply(X = sample_cl_s, FUN = extract_cl_d) %>% bind_rows()
  
   
  return(cl_boot_dat)
}

bootdat_test <- ran.gen_cluster(original_data = d_clean, clusters = d_clean$parcel_id)
  
  
  
fixed_effects <- "parcel_id + year"
#d_clean_sampled <- d_clean
#data <- d_clean
clusters = d_clean$parcel_id
reps <- 1000
reg1 <- est_sndstg
i <- 1
# cluster bootstrap function
clusbootreg <- function(formula, data, cluster, reps=1000){
  reg1 <- lm(formula, data)
  clusters <- names(table(cluster))
  sterrs <- matrix(NA, nrow=reps, ncol=length(coef(reg1)))
  for(i in 1:reps){
    # this picks ("samples") integers in 1:G, G times 
    index <- sample(1:length(clusters), length(clusters), replace=TRUE)
    aa <- clusters[index]
    bb <- table(aa)
    bootdat <- NULL
    # max bb is the maximum amount of times a same cluster has been sampled. 
    # j is the number of times a same cluster has been sampled. 
    for(j in 1:max(bb)){
      # so this selects rows of the data, that correspond to the clusters which names have 
      # been sampled j times. 
      cc <- data[cluster %in% names(bb[bb %in% j]),]
      for(k in 1:j){
        # and this binds these rows j times. 
        bootdat <- rbind(bootdat, cc)
      }
    }
    # stores the coeffs estimated on the bootstrap sample 
    sterrs[i,] <- coef(lm(formula, bootdat))
  }
  val <- cbind(coef(reg1),apply(sterrs,2,sd))
  colnames(val) <- c("Estimate","Std. Error")
  return(val)
}
# the coefficients and SE are equal to those in stata using 
# reghdfe wa_ffb_price_imp1_lag1 wa_prex_cpo_imp1_lag1 n_reachable_uml_lag1, absorb(parcel_id) residuals(v2h_fe)
# reghdfe lucpfip_pixelcount_total wa_ffb_price_imp1_lag1 v2h_fe n_reachable_uml_lag1, absorb(parcel_id) vce(robust)
# (they don't incorporate the constant in fe and therefore it is displayed, but that does not change anything.)

# Now see if we can replicate the boostrap SE. 
bootstraped <- boot(data = d_clean, 
                    statistic = ctrl_fun_endo,
                    fixed_effects = "parcel_id",
                    sim = "ordinary",
                    R = 200)


# now try to device the boot fucntion with ran.gen function that tells to 
# sample on clusters and not on single obs.  


norm.ci(bootstraped, index = 1) # equivalent to : 
norm.ci(t0 = bootstraped$t0[1],
        t = bootstraped$t[,1])

# replicating norm.ci (bias-corrected)
2*bootstraped$t0[1]-mean(bootstraped$t[,1])-qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,1]))
2*bootstraped$t0[1]-mean(bootstraped$t[,1])+qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,1]))

# non bias corrected CI
bootstraped$t0[1] - qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,1]))
bootstraped$t0[1] + qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,1]))

# replicating Stata (non bias-corrected): 
bootstraped$t0[1] - qnorm((1+0.95)/2)*0.0281059 
bootstraped$t0[1] + qnorm((1+0.95)/2)*0.0281059 








# nothing works, so try glm with dummy fixed effects
# takes > 10min 
glm_pois_fe <- glm(fml_mfe, 
                   family = "poisson", 
                   data = d)
# but yields the same coeff as feglm
round(glm_pois_fe$coefficients[2:9], 8) == round(est_pois$coefficients,8) 
# residuals (and fitted values are not the same length, because fixest removes fixed effects that have only
# zero outcomes - see NOTES after regression)
glm_pois_fe$residuals %>% length()
est_pois$residuals %>% length()
glm_pois_fe$fitted.values %>% length()
est_pois$fitted.values %>% length()

# we can identify these observations with 
only_o <- obs2remove(fml, data = d, family = "poisson")
attr(only_o, "fixef")$parcel_id %>% head()
only_o %>% head()
d$parcel_id %>% unique()%>% head()
length(only_o)
nrow(d)-nrow(d2)
d2 <- d[-only_o,]
est_pois2 <- fixest::feglm(fml, data = d,
                           family = "poisson")
only_o2 <- obs2remove(fml, data = d2, family = "poisson")
glmfe_sdw_vcov <- sandwich::vcovCL(glm_pois_fe, cluster = ~parcel_id)
coeftest(glm_pois_fe, glmfe_sdw_vcov) 
# wa_ffb_price_imp1_lag1        1.2816e-02  2.1179e-03  6.0511 1.438e-09 ***
# wa_ffb_price_imp1_3pya_lag1   7.1283e-03  6.1041e-03  1.1678  0.242893    
# wa_cpo_price_imp1_lag1        9.2404e-04  4.9445e-04  1.8688  0.061648 .  
# wa_cpo_price_imp1_3pya_lag1  -8.9259e-03  1.2863e-03 -6.9394 3.937e-12 ***


# now retrieve standard errors
pglm_vcov <- vcovCL(pglm_pois, cluster = d$parcel_id)
pglm_coeft<- coeftest(glm_pois, vcov. = pglm_vcov)

vcovCL(est_pois$coefficients, cluster = ~parcel_id)
etable(est_pois, se = "cluster")

cluster.vcov(pglm_pois)

# with clusterSEs
t <- cluster.im.glm(glm_pois, dat = d, cluster = ~parcel_id, drop = TRUE, return.vcv = TRUE, report = T)
### ### ### ###
class(est_pois) <- class(glm_pois)
names(glm_pois)%in%names(est_pois)

vcovCL(est_pois, cluster = ~parcel_id)
summary(est_pois) %>% class()
est_pois 
pglm_pois
pglm_pois %>% class()

vcov_t <- vcovCL(pglm_pois)

ft <- coeftest(pglm_pois, vcov. = vcov_t)


