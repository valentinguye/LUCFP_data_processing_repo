# reproduce script from 31012021 but with removal of indus-sm overlaps (straight in this script, pour ne pas tout cbhamlbouler avant )


##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table",
                   "foreign", "readstata13", "readxl",
                   "raster", "rgdal",  "sp", "spdep", "sf",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "msm", "car", "fixest", "sandwich", "lmtest", "boot", "multcomp",
                   "ggplot2", "leaflet", "htmltools")
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


# # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("leaflet", "leaflet.providers", "png")
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/reg_results")

### SET NUMBER OF THREADS USED BY {FIXEST} TO ONE (TO REDUCE R SESSION CRASH)
getFixest_nthreads()



### INDONESIAN CRS 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE
parcel_size <- 3000

# PIXEL AREA
# to rescale the average partial effects and the predictions from pixel counts to hectares
pixel_area_ha <- (27.8*27.6)/(1e4)


# OK SO NOW FIXED EFFECT VARIABLES ARE DEFINED IN add_parcel_variables.R AS geographic_year 
# the only change implied is that now district^year can also be specified as district_year.

### READ ALL POSSIBLE DATASETS HERE
# They are outputs of merge_lhs_rhs_parcels.R


d_50 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                 parcel_size/1000,"km_",
                                 "50CR.rds")))
# Split them into islands of interest
d_50_suma <- d_50[d_50$island == "Sumatra",]
d_50_kali <- d_50[d_50$island == "Kalimantan",]
#d_50_papu <- d_50[d_50$island == "Papua",]


d_30 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                 parcel_size/1000,"km_",
                                 "30CR.rds")))
d_30 <- d_30[,names(d_50)]
# Split them into islands of interest
d_30_suma <- d_30[d_30$island == "Sumatra",]
d_30_kali <- d_30[d_30$island == "Kalimantan",]
#d_30_papu <- d_30[d_30$island == "Papua",]

rm(d_30, d_50)

# d_30_all <- d_30[d_30$island %in% c("Sumatra", "Kalimantan", "Papua"),]
# d_50_all <- d_50[d_50$island %in% c("Sumatra", "Kalimantan", "Papua"),]
# 
# d_30_CA2 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
#                                  parcel_size/1000,"km_",
#                                  "30CR_CA2.rds")))
# # Split them into islands of interest
# d_30_CA2_suma <- d_30_CA2[d_30_CA2$island == "Sumatra",]
# d_30_CA2_kali <- d_30_CA2[d_30_CA2$island == "Kalimantan",]
# #d_30_papu <- d_30[d_30$island == "Papua",]
# 
# rm(d_30_CA2)
# 
# d_50_CA2 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
#                                  parcel_size/1000,"km_",
#                                  "50CR_CA2.rds")))
# # Split them into islands of interest
# d_50_CA2_suma <- d_50_CA2[d_50_CA2$island == "Sumatra",]
# d_50_CA2_kali <- d_50_CA2[d_50_CA2$island == "Kalimantan",]
# #d_50_papu <- d_50[d_50$island == "Papua",]
# 
# rm(d_50_CA2)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### REGRESSION FUNCTION ##### 

# Prefered specifications are passed as default arguments. 
# Currently, the returned object is a fixest object, of *ONE* regression.  
# d <- d_30_kali
catchment = "CR"
outcome_variable = "lucpfap_pixelcount"
island = "both"
start_year = 2002
end_year = 2014
alt_cr = FALSE
nearest_mill = FALSE
margin = "both"
restr_marg_def = TRUE
commo = c("cpo")
x_pya = 3
dynamics = FALSE
log_prices = TRUE
yoyg = FALSE
only_sr = FALSE
short_run = "full"
imp = 1
distribution = "quasipoisson"
fe = "parcel_id + district_year"#
offset = FALSE
lag_or_not = "_lag1"
controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml")#, "wa_prex_cpo_imp1""wa_pct_own_loc_gov_imp",
remaining_forest = FALSE
interaction_terms = NULL # "illegal2"  #c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1")
interact_regressors = TRUE
interacted = "regressors"
pya_ov = FALSE
illegal = "all"# "ill2" #
weights = FALSE
min_forest_2000 = 0
min_coverage = 0
output_full = FALSE
# 
# rm(catchment,outcome_variable,island,alt_cr,commo,x_pya,dynamics,log_prices,yoyg,short_run,imp,distribution,fe,remaining_forest,offset,lag_or_not,controls,interaction_terms ,interacted,pya_ov,illegal, nearest_mill, weights)

make_base_reg <- function(island,
                          start_year = 2002, 
                          end_year = 2014, 
                          outcome_variable = "lucpfip_pixelcount", # LHS. One of "lucfip_pixelcount", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount"p
                          catchment = "CR",  
                          alt_cr = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                          nearest_mill = FALSE, # whether the ibs variables should be attributed to parcels as from the nearest mill or inverse distance weighted average.  
                          margin = "both", # "both"
                          restr_marg_def = TRUE, # should the restrictive definition of extensive margin be used, if margin is either "intensive" of "extensive"
                          commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo"), commodities the price signals of which should be included in the RHS
                          x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                          dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                          yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                          only_sr = FALSE,
                          log_prices = TRUE, # Logical, should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
                          short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                          imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                          distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                          fe = "parcel_id + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                          offset = FALSE, # Logical. Should the log of the remaining forest be added as an offset.  
                          lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                          controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml"), # , "wa_prex_cpo_imp1"character vectors of names of control variables (don't specify lags in their names)
                          remaining_forest = FALSE, # Logical. If TRUE, the remaining forest is added as a control
                          interaction_terms = NULL, # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1"), # may be one or several of the controls specified above. 
                          interacted = "regressors",
                          interact_regressors = TRUE, # if there are two regressors (e.g. ffb and cpo), should their interaction be included in the model? 
                          pya_ov = FALSE, # logical, whether the lagged (by one year) outcome_variable should be added in controls
                          illegal = "all", # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                          min_forest_2000 = 0, 
                          min_coverage = 0, # fraction, from 0 to 1. Minimum share of reachable IBS over all reachable (UML), for an obs.to be included in sample. 
                          weights = FALSE, # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
                          output_full = FALSE
){
  
  ### BASE DATA 
  
  if(catchment == "CR"){
    # these are preloaded, because CR is the preferred data set
    # Catchment radius
    if(island == "Sumatra" & alt_cr == FALSE){
      d <- d_30_suma
    }
    if(island == "Sumatra" & alt_cr == TRUE){
      d <- d_50_suma
    }
    if(island == "Kalimantan" & alt_cr == FALSE){
      d <- d_50_kali
    }
    if(island == "Kalimantan" & alt_cr == TRUE){
      d <- d_30_kali
    }
    if(island == "both" & alt_cr == FALSE){
      d <- rbind(d_30_suma, d_50_kali)#, d_50_papu
    }
    if(island == "both" & alt_cr == TRUE){
      d <- rbind(d_50_suma, d_30_kali)#, d_30_papu
    }
  } else if(catchment == "CA"){
    # here we load it, it slows the function but we will use it only once or twice.
    d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                  parcel_size/1000,"km_",
                                  "2h_CA.rds")))
    if(island == "Sumatra"){
      d <- d[d$island == "Sumatra",]
    }
    if(island == "Kalimantan"){
      d <- d[d$island == "Kalimantan",]
    }
    
  }
  
  
  # ////!!!!\\\\\ remove overlap of indus and sm lucfp the same year 
  # 
  # d$lucp_i_and_sm_bar <- (d$lucpfip_pixelcount == 0 | d$lucpfsmp_pixelcount == 0)
  # d$luc_i_and_sm_bar <- (d$lucfip_pixelcount == 0 | d$lucfsmp_pixelcount == 0)

  ## make a variable that counts lucfp events on all types of plantations
  # d <- mutate(d, lucpfap_pixelcount = (lucpfip_pixelcount + lucpfsmp_pixelcount)*lucp_i_and_sm_bar)
  # d <- mutate(d, lucfap_pixelcount = (lucfip_pixelcount + lucfsmp_pixelcount)*luc_i_and_sm_bar)
  
  ### SPECIFICATIONS  
  
  if(offset & grepl("lucpf", outcome_variable)){offset_fml <- ~log(remain_pf_pixelcount)}
  #if(offset & grepl("lucf", outcome_variable)){offset_fml <- ~log(remain_f30th_pixelcount)}
  
  ## VARIABLES OF INTEREST (called regressors here)
  if(dynamics == FALSE){ 
    # Only the overall effect is investigated then, no short vs. medium run
    # In this case, the variables have name element _Xya_ with X the number of years over which the mean has been 
    # computed, always including the contemporaneous record. See add_parcel_variables.R
    #short_run <- ""
    
    # if we omit one commodity 
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- paste0(commo,"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not)
      }else{
        regressors <- paste0(commo,"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)
      }
    }
    # if we don't omit a commodity. 
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not))
        
      }else{
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)) 
      }
    }
  }
  
  # if, on the other hand, we want to identify the short run effects, controlling for the medium run's. 
  if(dynamics == TRUE){ 
    
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                        paste0(commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
      }
      if(yoyg == FALSE){
        if(short_run == "full"){
          regressors <- c(paste0(commo,"_price_imp",imp,lag_or_not), # SR measure
                          paste0(commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
        }
        if(short_run == "dev"){
          regressors <- c(paste0(commo,"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not), # SR measure
                          paste0(commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
        }  
      }
    }
    
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0(commo[1],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not))    
      }
      if(yoyg == FALSE){
        if(short_run == "full"){
          regressors <- c(paste0(commo[1],"_price_imp",imp,lag_or_not),# FFB SR measure
                          paste0(commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not), # FFB LR measure 
                          paste0(commo[2],"_price_imp",imp,lag_or_not),# CPO SR measure
                          paste0(commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # CPO LR measure 
        }
        if(short_run == "dev"){
          regressors <- c(paste0(commo[1],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                          paste0(commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not),
                          paste0(commo[2],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                          paste0(commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not))
        }
      }
    }
  }
  
  if(only_sr){
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- paste0(commo,"_price_imp",imp,"_yoyg",lag_or_not)
      }else{
        regressors <- paste0(commo,"_price_imp",imp,lag_or_not)
      }
    }
    # if we don't omit a commodity. 
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg",lag_or_not))
        
      }else{
        regressors <- c(paste0(commo[1],"_price_imp",imp,lag_or_not),
                        paste0(commo[2],"_price_imp",imp,lag_or_not)) 
      }
    }
  }
  
  if(nearest_mill==FALSE){
    regressors <- paste0("wa_",regressors)
  }
  
  # Logarithms
  if(log_prices){
    for(reg in regressors){
      d[,paste0("ln_",reg)] <- log(d[,reg])
    }
    regressors <- paste0("ln_", regressors)
    rm(reg)
  }
  
  ## CONTROLS
  # add lagged outcome variable 
  #if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}
  # for 1 year lag
  if(any(grepl(pattern = paste0(outcome_variable,"_lag1"), x = controls))){
    d <- dplyr::arrange(d, parcel_id, year)
    d <- DataCombine::slide(d,
                            Var = outcome_variable,
                            TimeVar = "year",
                            GroupVar = "parcel_id",
                            NewVar = paste0(outcome_variable,"_lag1"), # the name passed in controls should correspond. 
                            slideBy = -1,
                            keepInvalid = TRUE)
    d <- dplyr::arrange(d, parcel_id, year)
  }
  # for 4 year lag
  if(any(grepl(pattern = paste0(outcome_variable,"_lag4"), x = controls))){
    d <- dplyr::arrange(d, parcel_id, year)
    d <- DataCombine::slide(d,
                            Var = outcome_variable,
                            TimeVar = "year",
                            GroupVar = "parcel_id",
                            NewVar = paste0(outcome_variable,"_lag4"), # the name passed in controls should correspond. 
                            slideBy = -4,
                            keepInvalid = TRUE)
    d <- dplyr::arrange(d, parcel_id, year)
  }
  
  
  # lag controls that are from IBS 
  select_ibs_controls <- grepl("pct_own", controls) | grepl("prex_", controls)
  if(lag_or_not=="_lag1" & length(controls)>0 & any(select_ibs_controls)){
    # find them
    controls[select_ibs_controls] <- sapply(controls[select_ibs_controls], FUN = paste0, lag_or_not)
  }
  
  # add remaining forest or not
  if("remaining_forest" %in% controls){
    if(grepl("lucpf", outcome_variable)){
      controls[controls == "remaining_forest"] <- "remain_pf_pixelcount"
    }
    if(grepl("lucf", outcome_variable)){
      controls[controls == "remaining_forest"] <- "remain_f30th_pixelcount"
    }
    
    offset <- FALSE
  }
  
  # build baseline forest trend if it is featured in the controls. 
  if("baseline_forest_trend" %in% controls){
    if(grepl("lucpf", outcome_variable)){
      d[,"baseline_forest_trend"] <- d$pfc2000_total_pixelcount*(d$year-2000)
    }
    if(grepl("lucf", outcome_variable)){
      d[,"baseline_forest_trend"] <- d$fc2000_30th_pixelcount*(d$year-2000)
    }
  }
  
  if(nearest_mill==TRUE){
    controls <- gsub("wa_", "", controls)
  }
  
  #d$reachable_year <- d$reachable*(d$year - 2000)
  
  # ### WEIGHTS
  # if(weights){
  #   d$sample_coverage <- d$n_reachable_ibs/d$n_reachable_uml
  # }
  
  ### SELECT DATA FOR REGRESSION
  
  ## group all the variables necessary in the regression
  # important to do that after outcome_variable, regressors controls etc. have been (re)defined. 
  # (interactions do not need to be in there as they are fully built from the used_vars)
  used_vars <- c(outcome_variable, regressors, controls,
                 "parcel_id",  "year", "lat", "lon", 
                 "village", "subdistrict", "district", "province", "island", #"illegal2",
                 "village_year", "subdistrict_year", "district_year", "province_year")#,"reachable", "reachable_year"
  #"n_reachable_ibsuml_lag1", "sample_coverage_lag1", #"pfc2000_total_ha", 
  #"remain_f30th_pixelcount","remain_pf_pixelcount"
  if(nearest_mill){
    used_vars <- c(used_vars, "nearest_firm_id")
  }
  
  # if we use an offset (for now only remaining forest possible) then it should be in the necessary variables, but not among the controls. 
  if(offset){
    if(grepl("lucpf", outcome_variable)){
      used_vars <- c(used_vars, "remain_pf_pixelcount")
    }
    if(grepl("lucf", outcome_variable)){
      used_vars <- c(used_vars, "remain_f30th_pixelcount")
    }
  } 
  
  ### KEEP OBSERVATIONS THAT: 
  
  # - are in study period 
  d <- dplyr::filter(d, year >= start_year)
  d <- dplyr::filter(d, year <= end_year)
  
  # - are covered with some of the studied forest 
  # in 2000
  if(grepl("lucpf", outcome_variable)){
    d <- dplyr::filter(d, pfc2000_total_pixelcount*pixel_area_ha/900 > min_forest_2000)
    #d <- d[d$pfc2000_total_pixelcount*pixel_area_ha/900 > min_forest_2000,]
  }
  if(grepl("lucf", outcome_variable)){
    d <- dplyr::filter(d, fc2000_30th_pixelcount*pixel_area_ha/900 > min_forest_2000)
    # d <- d[d$fc2000_30th_pixelcount*pixel_area_ha/900 > min_forest_2000,]
  }
  
  # # in years after grid cells get completely deforested)
  if(grepl("lucpf", outcome_variable)){
    d <- dplyr::filter(d, remain_pf_pixelcount > 0,)
    #d <- d[d$remain_pf_pixelcount > 0,]
  }
  if(grepl("lucf", outcome_variable)){
    d <- dplyr::filter(d, remain_f30th_pixelcount > 0,)
    # d <- d[d$remain_f30th_pixelcount > 0,]
  }
  
  # Do not experience industrial and smallholder deforestation at the same time.
  if(grepl("lucpf", outcome_variable)){
    d <- d[d$lucp_i_and_sm_bar,]
  }
  if(grepl("lucf", outcome_variable)){
    d <- d[d$luc_i_and_sm_bar,]
  }
  
  
  # - are in years when the outcome can actually be observed
  if(grepl("_slow_",outcome_variable)){
    d <- dplyr::filter(d, year < 2011)
    d <- d[d$year<2011,]
  }
  
  # remove year 2015 as it is only available for industrial plantations
  # d <- dplyr::filter(d, year < 2015)
  #  d <- d[d$year<2015,]
  # if(grepl("fsmp_", outcome_variable) | grepl("fsp_", outcome_variable) | grepl("fmp_", outcome_variable)){
  #   d <- d[d$year<2015,]
  # }
  
  # - are in intensive margin expansion areas
  if(restr_marg_def){
    if(margin == "intensive"){
      d <- dplyr::filter(d, extensive_restr == FALSE)
      # d <- d[d$extensive_restr == FALSE,]
    }
    if(margin == "extensive"){
      d <- dplyr::filter(d, extensive_restr == TRUE)
      # d <- d[d$extensive_restr == TRUE,]
    }
  } else {
    if(margin == "intensive"){
      d <- dplyr::filter(d, extensive == FALSE)
      # d <- d[d$extensive == FALSE,]
    }
    if(margin == "extensive"){
      d <- dplyr::filter(d, extensive == TRUE)
      # d <- d[d$extensive == TRUE,]
    }
  }
  
  # - have a minimum reachable mill sample coverage (and a positive amount of reachable mills)
  d <- dplyr::filter(d, n_reachable_ibs > 0 &
                       n_reachable_uml > 0 & 
                       n_reachable_ibs/n_reachable_uml >= min_coverage)
  
  # d <- d[d$n_reachable_ibs > 0 &
  #        d$n_reachable_uml > 0 & 
  #        d$n_reachable_ibs/d$n_reachable_uml >= min_coverage,]
  
  
  # - are not in an intended RSPO-certified supply base
  d <- dplyr::filter(d, rspo_cert == FALSE)
  #d <- d[d$rspo_cert==FALSE,]
  
  # are legal, illegal, or neither
  if(illegal == "no_ill1"){
    d <- d[!is.na(d$illegal1) & d$illegal1 == FALSE, ]
  }
  if(illegal == "no_ill2"){
    d <- d[!is.na(d$illegal1) & d$illegal2 == FALSE, ]
  }
  if(illegal == "ill1"){
    d <- d[!is.na(d$illegal2) & d$illegal1 == TRUE, ]
  }  
  if(illegal == "ill2"){
    d <- d[!is.na(d$illegal2) & d$illegal2 == TRUE, ]
  }
  
  # - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
  usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
  # used_vars[!(used_vars%in%names(d))]
  names(usable) <- used_vars            
  usable <- bind_cols(usable)
  filter_vec <- base::rowSums(usable)
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d[filter_vec, c(used_vars)]
  if(anyNA(d_nona)){stop()}
  rm(filter_vec, usable)
  
  # # - sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
  # # probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
  # # just remove them if any, so that there is no bug. 
  # # hm no it's rather panel obs. that don't have any reachable mill YET. So they are important, but they should be
  # if(weights){
  #   d_nona <- d_nona[d_nona$sample_coverage!=0,]
  #   d_nona[d_nona$sample_coverage>1,"sample_coverage"] <- 1
  # }
  
  # OUTLIERS # 
  # otl <- boxplot.stats(d_nona[,regressors[1]])$out
  # d_nona <- d_nona[!(d_nona[,regressors[1]] %in% otl), ]
  
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_clean <- d_nona[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                                d_nona, 
                                family = "poisson"),]
  
  # note that obs2remove has to be the last filtering operation on data, otherwise some units (along fe dimension)
  # may become "only-zero-outcome" units after other data removal.  
  
  
  
  # make the interaction variables in the data
  
  ## INTERACTIONS
  interaction_vars <- c()
  
  # here we produce the names of the actual interaction variables
  if(length(interaction_terms)>0){
    # unless variables of interest to be interacted are specified, they are all the presently defined regressors
    if(interacted == "regressors"){interacted <- regressors}
    # this function selects the actual controls (lagged if necessary) into the interaction terms
    make_int_term <- function(x){int_term <- controls[grepl(x,controls)] %>% paste0("X",interacted)}
    interaction_vars <- sapply(interaction_terms, FUN = make_int_term)
  }else{interacted <- NULL}
  
  # add the interaction between the regressors 
  if(length(regressors) == 2 & interact_regressors){
    interaction_vars <- c(interaction_vars, paste0(regressors[1],"X",regressors[2]))
  }  
  
  
  # and here we build the interactions 
  if(length(interaction_terms)>0){
    for(con in interaction_terms){
      actual_ <- controls[grepl(con, controls)]
      for(reg in interacted){
        d_clean[,paste0(actual_,"X",reg)] <- d_clean[,actual_]*d_clean[,reg]
      }
    }
  }
  
  # and add the interaction between the regressors 
  if(length(regressors) == 2 & interact_regressors){
    d_clean[,paste0(regressors[1],"X",regressors[2])] <- d_clean[,regressors[1]]*d_clean[,regressors[2]]
    #d_clean[,paste0(regressors[2],"X",regressors[1])] <- d_clean[,regressors[2]]*d_clean[,regressors[1]]
  }
  
  
  ### REGRESSIONS
  
  # Model specification
  fe_model <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                paste0(regressors, collapse = "+"),
                                " + ",
                                paste0(controls, collapse = "+"),
                                " | ",
                                fe))
  
  if(length(interaction_terms)>0 | (interact_regressors & length(regressors) == 2)){
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ ",
                                  paste0(regressors, collapse = "+"),
                                  " + ",
                                  paste0(controls, collapse = "+"),
                                  " + ",
                                  paste0(interaction_vars, collapse = "+"),
                                  " | ",
                                  fe))
    # 
    # fe_model <- as.formula(paste0(outcome_variable,
    #                               " ~ ",
    #                               paste0(regressors, collapse = "+"),
    #                               " + ",
    #                               "interact(",regressors,",",interaction_terms,")",#",","FALSE",
    #                               " + ",
    #                               paste0(controls, collapse = "+"),
    #                               " | ",
    #                               fe))
  }
  
  
  # Estimation
  if(offset == TRUE){
    if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
      if(weights == TRUE){
        var_weights <- d_clean$sample_coverage
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution,
                                 offset = offset_fml,
                                 glm.iter = 200,
                                 notes = TRUE, 
                                 weights = var_weights)
        
      }else{
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution,
                                 offset = offset_fml,
                                 glm.iter = 200,
                                 notes = TRUE)
      }
    }else{ # no weights allowed in negative binomial
      reg_res <- fixest::fenegbin(fe_model,
                                  data = d_clean,
                                  offset = offset_fml,
                                  notes = TRUE)
    }
  }  
  
  if(offset == FALSE){
    if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
      if(weights == TRUE){
        var_weights <- d_clean$sample_coverage
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution, 
                                 glm.iter = 200,
                                 #fixef.iter = 100000,
                                 notes = TRUE, 
                                 weights = var_weights)
        
      }else{
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution, 
                                 glm.iter = 100,
                                 #fixef.iter = 100000,
                                 notes = TRUE)
      }
    }else{ # no weights allowed in negative binomial
      reg_res <- fixest::fenegbin(fe_model,
                                  data = d_clean,
                                  #fixef.iter = 100000,
                                  notes = TRUE)
    }
  }
  
  
  
  if(output_full){
    toreturn <- list(reg_res, d_clean, d)
  }else{
    toreturn <- list(reg_res, d_clean)
  }
  
  rm(d, d_nona)
  return(toreturn)
  
}

#### APE FUNCTIONS ####
# helper function that transforms the list of results into a data frame of average partial effects (APEs) and their standard errors (SEs), 
# for the K first regressors in the models fitted by make_base_reg 
# If there are interactions in the models, the APEs (and SEs) of the interaction effects are computed (may not work if K > 1 then)


# # If we want to compute the partial effect for a one standard deviation (or the equivalent relative change in price)
# d_clean <- left_join(d_clean, d[,c("parcel_id", "year", "wa_cpo_price_imp1_4ya_lag1")], by = c("parcel_id", "year")) 
# # this removes fixed effect variations before computing standard deviation. 
# reg_sd <- fixest::feols(fml = as.formula("wa_cpo_price_imp1_4ya_lag1 ~ 1 |parcel_id + district_year"), #  
#                           data = d_clean)
# reg_sd_ln <- fixest::feols(fml = as.formula("ln_wa_cpo_price_imp1_4ya_lag1 ~ 1 |parcel_id + district_year"), #  
#                         data = d_clean)
# sd <- sd(reg_sd$residuals)
# sd_ln <- sd(reg_sd_ln$residuals)
# m <- mean(d_clean$wa_cpo_price_imp1_4ya_lag1)
# m_ln <- mean(d_clean$ln_wa_cpo_price_imp1_4ya_lag1)
# sd/m
# sd_ln/m_ln

# res_data <- res_data_list_full[[1]]
# k = 1
# K=1
# controls_pe = F
# SE = "cluster"
# CLUSTER = "subdistrict"
# rel_price_change = 0.01 # sd/m #
# abs_price_change = 1
# stddev = FALSE
# rounding = 2

make_APEs <- function(res_data, K=1, 
                      controls_pe = FALSE, 
                      #SE = "cluster", 
                      CLUSTER = "subdistrict",
                      stddev = FALSE,
                      rel_price_change = 0.01, 
                      abs_price_change = 1, 
                      rounding = 2){
  
  # store APEs and their deltaMethod statistics in this list 
  dM_ape_roi_list <- list()
  
  # get estimation results and data 
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]]
  
  ## Redefine changes in regressor to one standard deviation if asked 
  if(stddev){
    # remove fixed effect variations from the regressor
    reg_sd <- fixest::feols(fml = as.formula(paste0(
      names(reg_res$coefficients)[1],
      " ~ 1 | ", 
      paste0(reg_res$fixef_vars, collapse = "+"))),
      data = d_clean)
    
    # and take the remaining standard deviation
    rel_price_change <- sd(reg_sd$residuals)
    abs_price_change <- sd(reg_sd$residuals)
  }
  
  
  # identify the nature of different variables
  coeff_names <- names(coef(reg_res))
  # define actual interaction terms (possibly lagged)
  interaction_effects <- coeff_names[grepl(pattern = "X", names(coef(reg_res)))]
  others <- coeff_names[!(grepl(pattern = "X", coeff_names))]
  interaction_terms <- others[paste0(others,"X",coeff_names[1]) %in% interaction_effects]
  
  
  # data POINTS that deltaMethod will grab 
  
  # averages of the interaction terms -the controls that we interacted with the regressor of interest) 
  int_term_avg <- list()
  if(length(interaction_terms) >0){
    for(i in 1:length(interaction_terms)){
      int_term_avg[[i]] <- mean(d_clean[,interaction_terms[i]])
    }
  }
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values) 
  
  
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  
  ## repeat the following for the K regressors of interest 
  linear_ape_fml_list <- list() # to store some results, see at end of loop
  for(k in 1:K){
    
    # add interaction term and average of the other regressor of interest if there are several (K == 2 only for now)
    if(K == 1){
      interaction_terms_k <- interaction_terms
      int_term_avg_k <- int_term_avg
    }
    
    if(K == 2 & k == 1){
      interaction_terms_k <- c(interaction_terms, coeff_names[2])
      int_term_avg_k <- c(int_term_avg, mean(d_clean[,coeff_names[2]]))
    }
    if(K == 2 & k == 2){
      interaction_terms_k <- c(interaction_terms, coeff_names[1])
      int_term_avg_k <- c(int_term_avg, mean(d_clean[,coeff_names[1]]))
    }
    
    
    linear_ape_fml <- paste0(coeff_names[k])
    
    # if there are interactions, the partial effect formula needs to be appended.
    i_t <- 1
    while(i_t<=length(interaction_terms_k)){
      linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[grepl(coeff_names[k],interaction_effects)][i_t],"*",int_term_avg_k[[i_t]])
      i_t <- i_t +1
    }
    rm(i_t)
    
    # the final formula is different depending on the regressor of interest being in the log scale or not. 
    if(grepl("ln_",coeff_names[k])){
      ape_fml_roi <- paste0("((",1+rel_price_change,")^(",linear_ape_fml,") - 1)*100")#*fv_bar*",pixel_area_ha)
    } else{
      ape_fml_roi <- paste0("(exp(",linear_ape_fml,"*",abs_price_change,") - 1)*100")#*fv_bar*",pixel_area_ha)
    } 
    
    dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                              vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                              g. = ape_fml_roi, 
                              rhs = 0)
    
    row.names(dM_ape_roi) <- NULL
    dM_ape_roi <- as.matrix(dM_ape_roi)
    dM_ape_roi <- dM_ape_roi[,c("Estimate","2.5 %","97.5 %")]
    
    dM_ape_roi_list[[k]] <- dM_ape_roi  
    
    
    ### COMPUTE PARTIAL EFFECTS OF INTERACTION TERMS 
    dM_ape_roi_list[[K + k]] <- list()
    
    if(length(interaction_terms) >0){
      # average of the regressor of interest
      reg_bar <- mean(d_clean[,coeff_names[k]]) 
      
      for(i in 1:length(interaction_terms)){
        
        # there are several parts 
        # 1. the coeff of the interaction effect of interest
        iei <- paste0(interaction_terms[i], "X", coeff_names[k])
        ape_fml1 <- paste0(iei," + ") 
        
        # 2. the linear expression derivative wrt. the regressor of interest 
        ape_fml2 <- paste0("(",linear_ape_fml,")")
        
        # 3. the derivative of the linear expression wrt. the interaction term; that multiplies the previous linear expression
        ape_fml3 <- paste0("*(",
                           coeff_names[coeff_names==interaction_terms[i]]," + ",
                           iei,"*",reg_bar,")")#
        
        # paste everything together; divide by 100 to approximate the log transformation of regressor of interest if it's in the log scale
        if(grepl("ln_", coeff_names[k])){
          ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*0.01")#,"*",fv_bar,"*",pixel_area_ha)
          
        } else{
          ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")")#,"*",fv_bar,"*",pixel_area_ha)
        }
        
        dM <- deltaMethod(object = coef(reg_res), 
                          vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE,
                          g. = ape_fml_it, 
                          rhs = 0)
        
        row.names(dM) <- NULL
        dM <- as.matrix(dM)
        dM <- dM[,c("Estimate","2.5 %","97.5 %")]
        
        dM_ape_roi_list[[K + k]][[i]] <- dM
        
      } # closes loop on i (interaction terms)
    } # closes condition on length(interaction_terms) > 0
    
    # save the linear derivative formula wrt. the regressor of interest
    linear_ape_fml_list[[k]] <- linear_ape_fml
  }# closes loop over k. 
  
  ### COMPUTE THE PARTIAL EFFECTS OF CONTROLS
  
  if(controls_pe & K == 1 & length(interaction_terms) == 0){ # we provide controls' partial effects only in the simplest case 
    for(coi in others[!grepl("price",others)]){
      # make the linear formula 
      ape_fml_coi <- paste0("(exp(",coi,"*",1,") - 1)*100")
      
      dM_ape_coi <- deltaMethod(object = coef(reg_res), 
                                vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                                g. = ape_fml_coi, 
                                rhs = 0)   
      
      row.names(dM_ape_coi) <- NULL
      dM_ape_coi <- as.matrix(dM_ape_coi)
      dM_ape_coi <- dM_ape_coi[,c("Estimate","2.5 %","97.5 %")]
      
      dM_ape_roi_list[[length(dM_ape_roi_list) + 1]] <- dM_ape_coi    
    }
  }
  
  # FINALLY, COMPUTE THE APE OF THE INTERACTION BETWEEN TWO REGRESSORS OF INTEREST
  if(K == 2){
    # there are several parts 
    # 1. the coeff of the interaction effect of interest
    iei <- paste0(coeff_names[1], "X", coeff_names[2]) # the order is important, as in the regression function, the interaction name follows this order with the two first regressors
    ape_fml1 <- paste0(iei," + ") 
    
    # 2. the linear expression derivative wrt. the regressor of interest 
    # linear_ape_fml will be the last computed one in the loop over k above, 
    # so if K = 2, it is computed as coeff_names[2] being the regressor with respect to which the expression is first derivate 
    ape_fml2 <- paste0("(",linear_ape_fml_list[[1]],")")
    
    # 3. and therefore, the derivative in this part should be wrt. coeff_names[1] that multiplies the previous linear expression
    ape_fml3 <- paste0("(",linear_ape_fml_list[[2]],")")
    
    # paste everything together; divide by 100 to approximate the log transformation of regressor of interest if it's in the log scale
    if(any(grepl("ln_", coeff_names[c(1,2)]))){
      ape_fml_it <- paste0("0.01*(",ape_fml1,ape_fml2, "*", ape_fml3, ")")#,"*",fv_bar,"*",pixel_area_ha)
      
    } else{
      ape_fml_it <- paste0("(",ape_fml1,ape_fml2, "*", ape_fml3, ")")#,"*",fv_bar,"*",pixel_area_ha)
    }
    
    dM <- deltaMethod(object = coef(reg_res), 
                      vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                      g. = ape_fml_it, 
                      rhs = 0)
    
    row.names(dM) <- NULL
    dM <- as.matrix(dM)
    dM <- dM[,c("Estimate","2.5 %","97.5 %")]
    dM_ape_roi_list[[length(dM_ape_roi_list) + 1]] <- dM
  }
  
  
  # make a one column matrix with all computed APEs' estimates, LB and HB values. 
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(dM_ape_roi_list)), 
                data = unlist(dM_ape_roi_list))  
  
  
  ## Shape the matrix of APE stats
  
  # round and coerce to character
  # mat[1,] <- mat[1,] %>% round(digits = 2)
  # mat[c(2,3),] <- mat[c(2,3),] %>% round(digits = 2)
  # mat[4,] <- mat[4,] %>% formatC(digits = 0, format = "f")
  # mat[5,] <- mat[5,] %>% formatC(digits = 0, format = "f")
  # 
  # mat[2,] <- paste0("[",mat[2,],"; ",mat[3,],"]")
  # mat <- mat[-3,] %>% as.matrix()
  
  mat <- round(mat, digits = rounding)
  
  k  <- 1
  while(k < nrow(mat)){
    mat[k+1,] <- paste0("[",mat[k+1,],"; ",mat[k+2,],"]")
    k <- k + 3
  } 
  row.names(mat) <- c(rep(c("Estimate","CI","delete"), nrow(mat)/3))
  mat <- mat[row.names(mat)!="delete",] %>% as.matrix()
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Observations"
  mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of clusters
  mat <- rbind(mat, length(unique(d_clean[,CLUSTER])))
  row.names(mat)[nrow(mat)] <- "Clusters"
  mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  # # add a row with the size of the variation 
  # if(grepl("ln_", coeff_names[k])){
  #   mat <- rbind(mat, rel_price_change)
  # } else{
  #   mat <- rbind(mat, abs_price_change)
  # }  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  
  
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, ape_fml_it, dM_ape_roi_list, linear_ape_fml_list)
  return(mat)
}


## REGRESSIONS
# infrastructure to store results
res_data_list_full <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")


for(SIZE in size_list){
  for(ILL in ill_status){
    res_data_list_full[[elm]] <- make_base_reg(island = ISL,
                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                               illegal = ILL,
                                               offset = FALSE)
    names(res_data_list_full)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }
}

## PARTIAL EFFECTS
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list_full, FUN = make_APEs)) %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat