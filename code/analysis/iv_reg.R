
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
                   "msm", "car", "fixest", "sandwich", "lmtest", "boot", "multcomp", "urca",
                   "ggplot2")#,"leaflet", "htmltools"
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

### SET NUMBER OF THREADS USED BY {FIXEST} TO 1 (might be neccessary to avoid R session crash) 
getFixest_nthreads()

### INDONESIAN CRS necessary for some spatial operations
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE (data available only for this size currently)
parcel_size <- 3000

# PIXEL AREA
# to rescale the average partial effects and the predictions from pixel counts to hectares
pixel_area_ha <- (27.8*27.6)/(1e4)

### Prepare polygons of three Indonesian islands of interest 
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"
island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

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


##### REGRESSION FUNCTION ##### 
# Commented out below are the arguments of the regression making function. 
# They may be useful to run parts of the operations within the function. 
island = "both"
start_year = 2002 
end_year = 2014 
outcome_variable = "lucpfsmp_pixelcount" # LHS. One of "lucfip_pixelcount" "lucfip_pixelcount_60th" "lucfip_pixelcount_90th" "lucpfip_pixelcount_intact" "lucpfip_pixelcount_degraded" "lucpfip_pixelcount"p
catchment = "CR"  
alt_cr = FALSE # logical if TRUE Sumatra's catchment radius is 50000 meters and Kalimantan's is 30000. 
nearest_mill = FALSE # whether the ibs variables should be attributed to parcels as from the nearest mill or inverse distance weighted average.  
margin = "both" # "both"
restr_marg_def = TRUE # should the restrictive definition of extensive margin be used if margin is either "intensive" of "extensive"
commo = "cpo" # either "ffb" "cpo" or c("ffb" "cpo") commodities the price signals of which should be included in the RHS
x_pya = 3 # either 2 3 or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
dynamics = FALSE # Logical should the total price signal(s) be split into current year and x_pya past year average. 
annual = FALSE # should every annual price signal be specified in the RHS (for years from current or lag1 if lag_or_not == "_lag1" to x_pya+1...)
price_variation = FALSE # should the regressors be price variation over the past years or average (the default)
yoyg = FALSE # logical should the price variables be computed in year-on-year growth rate instead of level.
only_sr = FALSE
log_prices = TRUE # Logical should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
short_run = "full" # either "full" or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
imp = 1 # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
distribution = "quasipoisson" # either "poisson" "quasipoisson" or "negbin"
fe = "reachable + district_year" # fixed-effects interactions should not be specified in {fixest} synthax with fe1^fe2
offset = FALSE # Logical. Should the log of the remaining forest be added as an offset.  
lag_or_not = "_lag1" # either "_lag1" or  "" should the 
controls = c("wa_avg_prex_cpo_imp1") #, "n_reachable_uml"  "wa_prex_cpo_imp1"character vectors of names of control variables (don't specify lags in their names)
remaining_forest = FALSE # Logical. If TRUE the remaining forest is added as a control
interaction_terms = NULL # c("wa_pct_own_nat_priv_imp""wa_pct_own_for_imp""n_reachable_uml" "wa_prex_cpo_imp1") # may be one or several of the controls specified above. 
interacted = "regressors"
interact_regressors = TRUE # if there are two regressors (e.g. ffb and cpo) should their interaction be included in the model? 
pya_ov = FALSE # logical whether the lagged (by one year) outcome_variable should be added in controls
illegal = "no_ill2" # the default "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions see add_parcel_variables.R)
min_forest_2000 = 0 
min_coverage = 0 # fraction from 0 to 1. Minimum share of reachable IBS over all reachable (UML) for an obs.to be included in sample. 
weights = FALSE # logical should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
output_full = FALSE # if TRUE the larger dataset (not only the one prepared foranalysis within the function) is returned as a third element of the list output 

IV_REG = TRUE # there is the option to set IV_REG to false just to test similarity with make_main_reg
instru_share = "avged"  
# APE parameters that were in APE function in standard final script but with IV reg APE and its SE needs to be computed within bootstrap estimation
CLUSTER = "reachable"#"subdistrict"
stddev = FALSE # if TRUE the PEs are computed for a one standard deviation (after removing variation in the fixed-effect dimensions)
rel_price_change = 0.01 
abs_price_change = 1 
rounding = 2 
boot_rep = 2
# 
# rm(catchment,outcome_variable,island,alt_cr,commo,x_pya,dynamics,log_prices,yoyg,short_run,imp,distribution,fe,remaining_forest,offset,lag_or_not,controls,interaction_terms ,interacted,pya_ov,illegal, nearest_mill, weights)


make_IV_reg <- function(island,
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
                        annual = TRUE, # should every annual price signal be specified in the RHS (for years from current, or lag1 if lag_or_not == "_lag1", to x_pya+1...)
                        price_variation = FALSE, # should the regressors be price variation over the past years, or average (the default)
                        yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                        only_sr = FALSE,
                        log_prices = TRUE, # Logical, should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
                        short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                        imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                        distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                        fe = "lonlat + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
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
                        output_full = FALSE, # if TRUE, the larger dataset (not only the one prepared foranalysis within the function) is returned as a third element of the list output 
                        
                        IV_REG = TRUE, # there is the option to set IV_REG to false, just to compare with regression without IV
                        instru_share = "avged",  
                        preclean_level = "FE", # on what dimensions should precleaning occur? Either same as fixed effects (FE) or another string value with "+" btwn dimension vars.
                        # APE parameters, that were in APE function in standard final script, but with IV reg, APE and its SE needs to be computed within bootstrap estimation
                        CLUSTER = "reachable",#"subdistrict",
                        stddev = FALSE, # if TRUE, the PEs are computed for a one standard deviation (after removing variation in the fixed-effect dimensions)
                        rel_price_change = 0.01, 
                        abs_price_change = 1, 
                        rounding = 2, 
                        boot_rep = 500
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
    
    # this is just for convenience, because this variable was not created for the CA stream, but is required to exist below 
    # (although wont be effectively used, see in specification charts). 
    d$reachable <- rep(1, nrow(d))
  }
  
  ### MAKE THE lonlat 
  d$lonlat <- d$lonlat
  
  # This would also be possible if done here, once we have stacked/merged all datasets to the final one. 
  # But changes nothing and is more intricate... 
  # d <- dplyr::arrange(d, lonlat, year)
  # uni_lonlat <- unique(d$lonlat)
  # d <- mutate(d,
  #             lonlat = match(lonlat, uni_lonlat))
  # rm(uni_lonlat)
  
  
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
  
  # if we want every annual price to be a regressor
  if(annual){
    regressors <- paste0(commo,"_price_imp",imp,"_lag",c(1:(x_pya+1)))
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
  
  if(price_variation){
    if(length(commo) == 1){
      regressors <- paste0(commo,"_price_imp",imp,"_",x_pya+1,"yv",lag_or_not)
    }
    if(length(commo) == 2){
      regressors <- c(paste0(commo[[1]],"_price_imp",imp,"_",x_pya+1,"yv",lag_or_not), 
                      paste0(commo[[2]],"_price_imp",imp,"_",x_pya+1,"yv",lag_or_not))
    }
  }
  
  
  
  if(nearest_mill==FALSE){
    regressors <- paste0("wa_",regressors)
  }
  
  # Logarithms
  if(log_prices & !price_variation){
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
    d <- dplyr::arrange(d, lonlat, year)
    d <- DataCombine::slide(d,
                            Var = outcome_variable,
                            TimeVar = "year",
                            GroupVar = "lonlat",
                            NewVar = paste0(outcome_variable,"_lag1"), # the name passed in controls should correspond. 
                            slideBy = -1,
                            keepInvalid = TRUE)
    d <- dplyr::arrange(d, lonlat, year)
  }
  # for 4 year lag
  if(any(grepl(pattern = paste0(outcome_variable,"_lag4"), x = controls))){
    d <- dplyr::arrange(d, lonlat, year)
    d <- DataCombine::slide(d,
                            Var = outcome_variable,
                            TimeVar = "year",
                            GroupVar = "lonlat",
                            NewVar = paste0(outcome_variable,"_lag4"), # the name passed in controls should correspond. 
                            slideBy = -4,
                            keepInvalid = TRUE)
    d <- dplyr::arrange(d, lonlat, year)
  }
  
  
  # lag controls that are from IBS 
  select_ibs_controls <- grepl("wa", controls) & !grepl("wa_avg_", controls) # | grepl("prex_", controls)
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
  
  ### INSTRUMENTS ### 
  if(IV_REG){
    
    instruments <- paste0("iv_",instru_share,"_imp",imp, "_lag", c(1:(x_pya+1)))
    
    controls <- c(controls, instruments)
  }
  
  # Log price variables that are among controls (in case of independent/integrated investigations)
  if(log_prices & !price_variation){
    for(reg in controls[grepl(pattern = "price", controls)]){ 
      d[,paste0("ln_",reg)] <- log(d[,reg])
      controls[match(reg, controls)] <- paste0("ln_", reg)
    }
    rm(reg)
  }
  
  
  # ### WEIGHTS
  # if(weights){
  #   d$sample_coverage <- d$n_reachable_ibs/d$n_reachable_uml
  # }
  
  ### SELECT DATA FOR REGRESSION
  
  ## group all the variables necessary in the regression
  # important to do that after outcome_variable, regressors controls etc. have been (re)defined. 
  # (interactions do not need to be in there as they are fully built from the used_vars)
  used_vars <- c(outcome_variable, regressors, controls,
                 "lonlat",  "year", "lat", "lon", 
                 "village", "subdistrict", "district", "province", "island", "reachable", # "illegal2",
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
  
  # # Do not experience industrial and smallholder deforestation at the same time. 
  # if(grepl("lucpf", outcome_variable)){
  #   d <- d[d$lucp_i_and_sm_bar,]
  # }
  # if(grepl("lucf", outcome_variable)){
  #   d <- d[d$luc_i_and_sm_bar,]
  # }
  
  
  # - are in years when the outcome can actually be observed
  # if(grepl("_slow_",outcome_variable)){
  #   d <- dplyr::filter(d, year < 2011)
  #   d <- d[d$year<2011,]
  # }
  
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
  
  # PRE - CLEAN 
  # After lots of consideration, we DO pre-clean the data, based on second stage feglm criteria, 
  # i.e. remove always zero a FE dimension, BEFORE the first stage. 
  # note that this has to be the last filtering operation on data, otherwise some units may become "only-zero-outcome" units after other data removal. 
  if(preclean_level == "FE"){
    preclean_level <- fe
  }
  temp_est <- feglm(fml = as.formula(paste0(outcome_variable, " ~ 1 | ", preclean_level)),
                    data = d_nona,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    d_precleaned <- d_nona[unlist(temp_est$obs_selection),]
  }  else { 
    d_precleaned <- d_nona
  }
  
  # save these for later
  avg_defo_ha <- (fitstat(temp_est, "my")[[1]]*pixel_area_ha) %>% round(1)
  
  G <- length(unique(d_precleaned[,CLUSTER]))
  
  
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
        d_precleaned[,paste0(actual_,"X",reg)] <- d_precleaned[,actual_]*d_precleaned[,reg]
      }
    }
  }
  
  # and add the interaction between the regressors 
  if(length(regressors) == 2 & interact_regressors){
    d_precleaned[,paste0(regressors[1],"X",regressors[2])] <- d_precleaned[,regressors[1]]*d_precleaned[,regressors[2]]
    #d_precleaned[,paste0(regressors[2],"X",regressors[1])] <- d_precleaned[,regressors[2]]*d_precleaned[,regressors[1]]
  }
  
  
  
  ### FORMULAE ###
  
  # INTRUMENTS ARE IN CONTROLS   
  # store first stage formulae in this list 
  fml_1st_list <- list()
  classic_iv_fixest_fml <- list()
  for(f in 1:length(regressors)){
    fml_1st_list[[f]] <- as.formula(paste0(regressors[f],
                                           " ~ ", 
                                           paste0(controls, collapse = "+"), #  instruments are in controls
                                           " | ", 
                                           fe))
    
    # this is not for within bootstraps, but to fit in feols-iv in fixest package to extract 1st stage stats readily
    classic_iv_fixest_fml[[f]] <- as.formula(paste0(outcome_variable,
                                                    " ~ ", 
                                                    paste0(controls[!(controls %in% instruments)], collapse = "+"), # only "included" controls, i.e. z1, i.e. not instruments.
                                                    " | ",
                                                    fe, 
                                                    " | ",
                                                    regressors[f], 
                                                    " ~ ",
                                                    paste0(instruments, collapse = "+")))
    
  }
  
  # need to specify the names of first stages' residuals now
  errors_1stg <- paste0("est_res_1st_", c(1:length(regressors)))
  
  fml_2nd <- as.formula(paste0(outcome_variable,
                               " ~ ", 
                               paste0(regressors, collapse = "+"),
                               " + ", 
                               paste0(errors_1stg, collapse = "+"),
                               " + ",
                               paste0(controls[!(controls %in% instruments)], collapse = "+"), # only "included" controls, i.e. z1, i.e. not instruments.
                               " | ",
                               fe))
  
  # ESTIMATION 
  if(IV_REG){
    
    # We run the first stage outside the bootstrapping process, in order to extract related information
    est_1st_list <- list()
    Ftest_list <- list()
    Waldtest_list <- list()
    for(f in 1:length(regressors)){
      est_1st <- fixest::feols(classic_iv_fixest_fml[[f]], 
                               data = d_precleaned)
      
      # Extract first stage information
      est_1st_list[[f]] <- summary(est_1st, .vcov = vcov(est_1st, cluster = CLUSTER) )#se = "cluster", cluster = CLUSTER
      
      Ftest_list[[f]] <- fitstat(est_1st, type = "ivf1")
      Waldtest_list[[f]] <- fitstat(est_1st, type = "ivwald1")
      
      # save residuals 
      # d_precleaned[,errors_1stg[f]] <- summary(est_1st, stage = 1)$residuals 
    }
    
    toreturn <- list(firstg_summary = est_1st_list, 
                     # fitstat not available in this version of fixest. Might want to compute a F-test by hand ... wald_jointnull_1st = fitstat(est_1st_list[[f]], "wald"), # F-test not available for GLM. 
                     # sum_fv_1st = sum(est_1st_list$fitted.values), # leave it in bare unit and convert/scale post estimation 
                     boot_info = NA,
                     APE_mat = NA) # this will be filled below
    
    Fstat <- unlist(Ftest_list)[[1]] %>% round(2)
    Waldstat <- unlist(Waldtest_list)[[1]] %>% as.numeric() %>% round(2)
    
    
    # SECOND STAGE ESTIMATION   
    # est_2nd <- fixest::feglm(fml_2nd,
    #                          data = d_precleaned,
    #                          family = distribution,
    #                          glm.iter = 100,
    #                          notes = TRUE)
    # summary(est_2nd, se = "cluster", cluster = CLUSTER)  
    
    
    ## Redefine changes in endogenous variable (first stage outcome) to one standard deviation if asked 
    if(stddev){
      # remove fixed effect variations from the regressor
      reg_sd <- fixest::feols(fml = as.formula(paste0(
        regressors[1],
        " ~ 1 | ", 
        paste0(est_1st$fixef_vars, collapse = "+"))),
        data = d_precleaned)
      
      # and take the remaining standard deviation
      rel_lu_change <- sd(reg_sd$residuals)
      abs_lu_change <- sd(reg_sd$residuals)
    }
    
    # function applied to each bootstrap replicate
    # myfun_data <- d_precleaned
    # fsf_list <- fml_1st_list
    # ssf <- fml_2nd
    ctrl_fun_endo <- function(myfun_data, fsf_list, ssf){
      
      ape <- try({
        # As many first stages as there are endogenous variables
        for(f in 1:length(regressors)){
          BS_est_1st <- fixest::feols(fsf_list[[f]], 
                                      data = myfun_data)
          
          # save residuals 
          myfun_data[,errors_1stg[f]] <- BS_est_1st$residuals
        }
        
        # 2nd stage
        BS_est_2nd <- fixest::feglm(ssf, 
                                    data = myfun_data, 
                                    family = distribution, 
                                    glm.iter = 100,
                                    notes = TRUE)
        
        ## MAKE APE 
        # Unlike in the non-IV case, APEs are not estimated with the delta Method, because their SEs are only computed asymptotically, by bootstrap. 
        # we don't need to compute inference statistics on them, as inference is handled by bootstrap (on APE).
        
        # select coefficient of interest 
        roi <- BS_est_2nd$coefficients[regressors]
        annual_ape <- c()
        i <- 1
        for(r in roi){
          # the APE is different depending on the regressor of interest being in the log scale or not. 
          if(grepl("ln_",regressors[1])){
            annual_ape[i] <- ((1+rel_price_change)^(r) - 1)*100
          } else{
            annual_ape[i] <- (exp(r*abs_price_change) - 1)*100 
          } 
          i <- i+1
        }  
        
        # return both cumulative and annual APE as statistics to bootstrap (i.e. don't force to decide here which is intersting, because it's as cotsly to compute both or only one)
        ape <- c(sum(annual_ape), annual_ape)
        ape
      })
      
      # if fixest returns an error in either stage in this boot sample   
      if(!is.numeric(ape)){
        ape <- NA
      }
      
      # statistics we want to evaluate the variance of:
      return(ape)
    }
    
    ## BOOTSTRAP
    
    # get the different cluster sizeS. This is necessary to cluster bootstrapping with clusters of different sizes. 
    sizes <- table(d_precleaned[,CLUSTER])
    u_sizes <- sort(unique(sizes))
    
    # names and numbers of clusters of every sizes
    cl_names <- list()
    n_clusters <- list()
    for(s in u_sizes){
      cl_names[[s]] <- names(sizes[sizes == s]) # these names are unique, by construction of using table() for sizes
      n_clusters[[s]] <- length(cl_names[[s]])
    }
    
    par_list <- list(unique_sizes = u_sizes,
                     cluster_names = cl_names,
                     number_clusters = n_clusters)
    
    # helper function
    # original_data <- d_precleaned
    # arg_list <- par_list
    ran.gen_cluster <- function(original_data, arg_list){
      
      # to store 
      cl_boot_dat <- NULL
      
      # non-unique names of clusters (repeated when there is more than one obs. in a cluster) 
      nu_cl_names <- as.character(original_data[,CLUSTER]) 
      
      for(s in arg_list[["unique_sizes"]]){
        # sample, in the vector of names of clusters of size s, as many draws as there are clusters of that size, with replacement
        sample_cl_s <- sample(arg_list[["cluster_names"]][[s]], 
                              arg_list[["number_clusters"]][[s]], 
                              replace = TRUE) 
        
        # because of replacement, some names are sampled more than once
        sample_cl_s_tab <- table(sample_cl_s)
        
        # because of replacement, some names are sampled more than once
        # we need to give them a new cluster identifier, otherwise a cluster sampled more than once 
        # will be "incorrectly treated as one large cluster rather than two distinct clusters" (by the fixed effects) (Cameron and Miller, 2015)    
        # here we do not necessarily need to bother with this, as clustering at the set of reachable mills is not a FE dimension. 
        sample_cl_tab <- table(sample_cl_s_tab)
        
        for(n in 1:max(sample_cl_s_tab)){ # from 1 to the max number of times a name was sampled bc of replacement
          # vector to select obs. that are within the sampled clusters. 
          names_n <- names(sample_cl_s_tab[sample_cl_s_tab == n])
          sel <- nu_cl_names %in% names_n
          
          # select data accordingly to the cluster sampling (duplicating n times observations from clusters sampled n times)
          clda <- original_data[sel,][rep(seq_len(sum(sel)), n), ]
          
          #identify row names without periods, and add ".0" 
          row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
          
          # add the suffix due to the repetition after the existing cluster identifier. 
          clda[,CLUSTER] <- paste0(clda[,CLUSTER], sub(".*\\.","_",row.names(clda)))
          
          # stack the bootstrap samples iteratively 
          cl_boot_dat <- rbind(cl_boot_dat, clda)       
        }
      }  
      return(cl_boot_dat)
    }
    
    # # test that 
    # test_boot_d <- ran.gen_cluster(original_data = d_precleaned,
    #                                arg_list = par_list)
    # dim(test_boot_d)
    # dim(d_precleaned)
    
    ## RUN BOOTSTRAP PROCEDURE
    
    # bootstrap on d_precleaned, no specific subset, because bootstrapping WILL lead to different data sets that have different 
    # patterns of always zero units. 
    set.seed(8888)
    bootstraped_1 <- boot(data = d_precleaned, 
                          statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                          # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                          fsf_list = fml_1st_list,
                          ssf = fml_2nd,
                          ran.gen = ran.gen_cluster,
                          mle = par_list,
                          sim = "parametric", # Note, from ?boot : "Use of sim = "parametric" with a suitable ran.gen allows the user to implement any types of nonparametric resampling which are not supported directly."
                          R = boot_rep)
    
    # bootstraped_1$t is a matrix with as many rows as bootrstrap replicates, and one column for each bootstrap statistic. 
    
    # store final information and bootstrap information (for checks) separately
    APE_estimand <- matrix(nrow = ncol(bootstraped_1$t), ncol = 4)
    colnames(APE_estimand) <- c("Estimate", "SE", "2.5 %", "97.5 %")
    output_list <- list()
    # Bootstrap info 
    boot_info <-  matrix(nrow = ncol(bootstraped_1$t), ncol = 5)
    colnames(boot_info) <- c("Estimate", "bootstrap_bias", "boot bias to coeff", "SE", "p_value")
    
    for(n_reg in 1:ncol(bootstraped_1$t)){ # 
      SE <- sd(bootstraped_1$t[,n_reg])
      
      beta_hat <- bootstraped_1$t0[n_reg]
      
      tval <- (beta_hat - 0)/SE
      
      # find degrees of freedom with "min" method: the number of clusters G, minus one.  
      # this follows the default in fixest package, as of version 0.7.0. See https://lrberge.github.io/fixest/articles/standard_errors.html
      # and from further version of the package, I verified that degrees_freedom(est_object, type = "t") indeed returns G-1
      boot_info[n_reg, "p_value"] <- (2*pt(abs(tval),
                                           lower.tail = FALSE,
                                           df = G-1))
      boot_info[n_reg, "bootstrap_bias"] <- mean(bootstraped_1$t[,n_reg]) - bootstraped_1$t0[n_reg] 
      
      boot_info[n_reg, "Estimate"] <- beta_hat
      boot_info[n_reg, "boot bias to coeff"] <- boot_info[n_reg, "bootstrap_bias"] / beta_hat
      boot_info[n_reg, "SE"] <- SE
      
      
      APE_estimand[n_reg, "Estimate"] <- beta_hat
      APE_estimand[n_reg, "SE"] <- SE
      APE_estimand[n_reg, "2.5 %"] <- beta_hat - qt(0.975, lower.tail=T, df=(G-1)) * SE
      APE_estimand[n_reg, "97.5 %"] <- beta_hat + qt(0.975, lower.tail=T, df=(G-1)) * SE
      
      output_list[[n_reg]] <- APE_estimand[n_reg, ]
    }
    
    # just because it is doubled when there is only one regressor
    if(length(regressors)==1){ 
      boot_info <- boot_info[1,]
      output_list <- output_list[1]
    }
    
    toreturn[["boot_info"]] <- boot_info
  }
  
  # this is just for comparison purpose
  if(!IV_REG){
    # Model specification
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ ",
                                  paste0(regressors, collapse = "+"),
                                  " + ",
                                  paste0(controls[!(controls %in% instruments)], collapse = "+"),
                                  " | ",
                                  fe))
    
    
    
    reg_res <- fixest::feglm(fe_model,
                             data = d_precleaned, 
                             family = distribution, 
                             glm.iter = 100,
                             #fixef.iter = 100000,
                             notes = TRUE)
    # select coefficient of interest 
    roi <- names(reg_res$coefficients[regressors])
    annual_ape_list <- list()
    i <- 1
    for(r in roi){
      # the final formula is different depending on the regressor of interest being in the log scale or not. 
      if(grepl("ln_",r)){
        ape_fml_roi <- paste0("((",1+rel_price_change,")^(",r,") - 1)*100")#*fv_bar*",pixel_area_ha)
      } else{
        ape_fml_roi <- paste0("(exp(",r,"*",abs_price_change,") - 1)*100")#*fv_bar*",pixel_area_ha)
      } 
      
      dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                                vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                                g. = ape_fml_roi, 
                                rhs = 0)
      
      row.names(dM_ape_roi) <- NULL
      dM_ape_roi <- as.matrix(dM_ape_roi)
      output_list[[match(r, roi)]] <- dM_ape_roi[,c("Estimate", "SE", "2.5 %","97.5 %")]
    }
    
    Fstat <- ""
    Waldstat <- ""
  } 
  ### OUTPUT ### 
  
  # make a one column matrix with all computed APEs' estimates, LB and HB values. 
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(output_list)), 
                data = unlist(output_list))  
  
  mat <- round(mat, digits = rounding)
  
  k  <- 1
  while(k < nrow(mat)){
    mat[k+2,] <- paste0("[",mat[k+2,],"; ",mat[k+3,],"]")
    k <- k + 4
  } 
  row.names(mat) <- c(rep(c("Estimate","SE", "CI","delete"), nrow(mat)/4))
  mat <- mat[row.names(mat)!="delete",] %>% as.matrix()
  
  # Average dependent variable 
  mat <- rbind(mat, avg_defo_ha)
  row.names(mat)[nrow(mat)] <- "Average deforestation (ha)"
  
  # add a row with the number of clusters
  mat <- rbind(mat, G)
  row.names(mat)[nrow(mat)] <- "Clusters"
  mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of observations
  mat <- rbind(mat, temp_est$nobs)
  row.names(mat)[nrow(mat)] <- "Observations"
  mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # First stage info
  mat <- rbind(mat, Fstat)
  row.names(mat)[nrow(mat)] <- "1st stage F-stat"
  mat <- rbind(mat, Waldstat)
  row.names(mat)[nrow(mat)] <- "1st stage joint Wald-stat"
  
  toreturn[["APE_mat"]] <- mat
  
  
  rm(d, d_nona, d_precleaned)
  return(toreturn)
  
}

### REGRESSIONS ### 
# infrastructure to store results
res_iv_avged_indep <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")


for(SIZE in size_list){
  for(ILL in ill_status){
    res_iv_avged_indep[[elm]] <- make_IV_reg(island = ISL,
                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                               illegal = ILL,
                                               commo = "cpo",
                                               controls = c("wa_avg_prex_cpo_imp1", "wa_cpo_price_imp1_4ya_lag1", 
                                                            "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp", "n_reachable_uml"), 
                                               instru_share = "avged",
                                               annual = TRUE, # if annual = FALSE, everything works well, just the APE_estimand dataframe outputed will have identical first two rows.
                                               boot_rep = 1,
                                               offset = FALSE)
    names(res_iv_avged_indep)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }
}
ape <- bind_cols(lapply(res_iv_avged_indep, FUN = function(x){x[["APE_mat"]]})) %>% as.matrix()

ape_mat <- bind_cols(lapply(res_iv_lagged, FUN = function(x){x[["APE_mat"]]})) %>% as.matrix()

boot_info <- lapply(res_iv_lagged, FUN = function(x){x[["boot_info"]]})

first_stage_info <- lapply(res_iv_lagged, FUN = function(x){x[[1]]})
res_iv_lagged[["both_i_all"]][["APE_mat"]]

ape_mat # this is instru share = lagged ; imp1

bind_cols(lapply(res_iv_lagged, FUN = function(x){x[["APE_mat"]]})) %>% as.matrix() # this is instru share = lagged ; imp2

bind_cols(lapply(res_data_list_full, FUN = function(x){x[["APE_mat"]]})) %>% as.matrix() # this is instru share = avged ; imp1

bind_cols(lapply(res_iv_avged_imp2, FUN = function(x){x[["APE_mat"]]})) %>% as.matrix() # this is instru share = avged ; imp2

### FFB
bind_cols(lapply(res_iv_avged_imp1_ffb, FUN = function(x){x[["APE_mat"]]})) %>% as.matrix() 

lapply(res_iv_avged_imp1_ffb, FUN = function(x){x[[1]]})



# -------------------------------------------------------------------------------------------------------------------------------------

#### APE FUNCTIONS ####
# helper function that transforms the list of results into a data frame of average partial effects (APEs) and their standard errors (SEs), 
# for the K first regressors in the models fitted by make_base_reg 
# If there are interactions in the models, the APEs (and SEs) of the interaction effects are computed (may not work if K > 1 then)

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
                      CLUSTER = "reachable",#"subdistrict",
                      stddev = FALSE, # if TRUE, the PEs are computed for a one standard deviation (after removing variation in the fixed-effect dimensions)
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
  mat <- rbind(mat, G)
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

# note that this function yields a warning that is not important: 
# In rm(coeff_names, interaction_effects, others, interaction_terms,  :
# objet 'ape_fml_it' introuvable

# this is a simpler version of make_APEs function above, that does not handle interaction terms. 
make_APEs_1regr <- function(res_data, 
                            #SE = "cluster", 
                            stddev = TRUE,
                            CLUSTER = "reachable",# "subdistrict", 
                            rel_price_change = 0.01, 
                            abs_price_change = 1){
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
  
  linear_ape_fml <- paste0(coeff_names[1])
  i_t <- 1
  while(i_t<=length(interaction_terms)){
    linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[grepl(coeff_names[1],interaction_effects)][i_t],"*",int_term_avg[[i_t]])
    i_t <- i_t +1
  }
  
  # the final formula is different depending on the regressor of interest being in the log scale or not. 
  if(grepl("ln_",coeff_names[1])){
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
  dM_ape_roi <- dM_ape_roi[,c("Estimate","SE","Pr(>|z|)")]
  
  dM_ape_list <- list()
  dM_ape_list[[1]] <- dM_ape_roi
  
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(dM_ape_list)), 
                data = unlist(dM_ape_list))  
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, dM_ape_list)
  return(mat)
}

#### PLANTATION DESCRIPTIVE STATISTICS WITH AND WITHOUT NAs #### 
make_desstats <- function(sample_1, sample_2){
  
  # get estimation results and data 
  # this is the exact set of observations used for estimation
  # #reg_res <- res_data[[1]]
  # parcels <- sample_1
  # d <- res_data[[3]]
  # sample_1 <- res_data_list_des[["both_i"]][[2]]
  # sample_2 <- res_data_list_des[["both_sm"]][[2]]
  
  
  # add variables that were not used in the regressions
  d_all <- rbind(d_30_suma, d_50_kali)
  # determine dependent variable
  dep_var_1 <- names(sample_1)[grepl("pixelcount",names(sample_1))]
  dep_var_2 <- names(sample_2)[grepl("pixelcount",names(sample_2))]
  
  # it is important that both have the same name
  # dep_var_1 is always length 1, but dep_var_2 is not if sample_2 is d
  if(length(dep_var_2)>1){
    dep_var_2 <- dep_var_1
  }
  
  
  
  
  
  # group public ownership
  #d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  
  # process a bit sample_2, in case it's the with-NA sample. 
  # identify grid cells with null lucfp every year
  sample_2 <- sample_2[-obs2remove(fml = as.formula(paste0(dep_var_2, " ~ lonlat + district_year")),
                                   sample_2[!is.na(sample_2[,dep_var_2]),], 
                                   family = "poisson"),]
  
  variables_1 <- c(dep_var_1,
                   "wa_cpo_price_imp1_4ya_lag1",
                   "wa_pct_own_gov_imp_lag1", "wa_pct_own_nat_priv_imp_lag1", "wa_pct_own_for_imp_lag1",
                   "n_reachable_uml") 
  variables_2 <- c(dep_var_2,
                   "wa_cpo_price_imp1_4ya_lag1",
                   "wa_pct_own_gov_imp_lag1", "wa_pct_own_nat_priv_imp_lag1", "wa_pct_own_for_imp_lag1",
                   "n_reachable_uml") 
  
  sample_1 <- left_join(sample_1[,c("lonlat","year")], d_all[,c("lonlat","year",variables_1)], by = c("lonlat","year"))#
  sample_2 <- left_join(sample_2[,c("lonlat","year")], d_all[,c("lonlat","year",variables_2)], by = c("lonlat","year"))
  
  # pixelcount to hectares
  sample_1[,dep_var_1] <- sample_1[,dep_var_1]*pixel_area_ha
  sample_2[,dep_var_2] <- sample_2[,dep_var_2]*pixel_area_ha
  
  
  statistics <- c("mean", "std.dev", "median", "min", "max")
  
  ## Des. stats. for sample_1
  des_sample_1 <- matrix(NA, nrow = length(variables_1), ncol = length(statistics))
  row.names(des_sample_1) <- variables_1
  colnames(des_sample_1) <- c(statistics)
  
  for(var in variables_1){
    des_sample_1[var,statistics] <- summarise(sample_1,
                                              mean = mean(get(var), na.rm=TRUE),
                                              std.dev = sd(get(var), na.rm= TRUE),
                                              median = median(get(var), na.rm= TRUE), 
                                              min = min(get(var), na.rm= TRUE),
                                              max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    des_sample_1[var, "median"] <- paste0(des_sample_1[var,"median"]," [",des_sample_1[var, "min"],"; ",des_sample_1[var,"max"],"]")
  }
  
  des_sample_1 <- des_sample_1[,c("mean", "std.dev", "median")]
  
  length(unique(sample_1$lonlat)) %>% paste0(" number of grid cells in sample") %>%  print()
  nrow(sample_1) %>% paste0(" number of observations in sample") %>%  print()
  
  
  ## Des. stats for sample_2
  des_sample_2 <- matrix(NA, nrow = length(variables_2), ncol = length(statistics))
  row.names(des_sample_2) <- variables_2
  colnames(des_sample_2) <- c(statistics)
  
  for(var in variables_2){
    des_sample_2[var,statistics] <- summarise(sample_2,
                                              mean = mean(get(var), na.rm=TRUE),
                                              std.dev = sd(get(var), na.rm= TRUE),
                                              median = median(get(var), na.rm= TRUE), 
                                              min = min(get(var), na.rm= TRUE),
                                              max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    des_sample_2[var, "median"] <- paste0(des_sample_2[var,"median"]," [",des_sample_2[var, "min"],"; ",des_sample_2[var,"max"],"]")
  }
  
  des_sample_2 <- des_sample_2[,c("mean", "std.dev", "median")]
  
  length(unique(sample_2$lonlat)) %>% paste0(" total number of grid cells") %>%  print()
  nrow(sample_2) %>% paste0(" total number of observations") %>%  print()
  
  
  ### add the t-test column
  t_test <- matrix(NA, nrow = length(variables_1), ncol = 1)
  row.names(t_test) <- variables_1
  colnames(t_test) <- "t test"
  
  for(var in variables_1[-1]){ # do not include the first one as it would not be the right one for sample_2
    test <- t.test(x = sample_1[,var],
                   y = sample_2[,var], 
                   alternative = "two.sided", 
                   mu = 0, 
                   paired = F, 
                   var.equal = FALSE)
    
    t_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  # t test on dependent variables
  test <- t.test(x = sample_1[,dep_var_1],
                 y = sample_2[,dep_var_2], 
                 alternative = "two.sided", 
                 mu = 0, 
                 paired = F, 
                 var.equal = FALSE)
  
  t_test[dep_var_1,] <- test$p.value %>% formatC(digits = 3, format = "f")
  # interpretation: if the 95% CI includes 0, then the difference in means between two samples 
  # is not statistically different from 0. Hence, the two samples are "similar" in means 
  # with respect to the variable tested. 
  # (In other words, we cannot reject the null hypothesis that the difference is null
  # -i.e. the two samples are alike - with 95% confidence) 
  
  ### add the Smirnov test
  ks_test <- matrix(NA, nrow = length(variables_1), ncol = 1)
  row.names(ks_test) <- variables_1
  colnames(ks_test) <- "KS test"
  
  for(var in variables_1[-1]){
    test <- ks.test(x = sample_1[,var],
                    y = sample_2[,var], 
                    alternative = "two.sided", 
                    exact = FALSE)
    
    ks_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  test <- ks.test(x = sample_1[,dep_var_1],
                  y = sample_2[,dep_var_2], 
                  alternative = "two.sided", 
                  exact = FALSE)
  
  ks_test[dep_var_1,] <- test$p.value %>% formatC(digits = 3, format = "f")
  # intepretation: if p-value < 0.05 we cannot reject with 95% confidence that 
  # two distributions are equal, implying that they are different. 
  
  
  ### Bind all 
  des_table <- cbind(des_sample_1, des_sample_2, t_test, ks_test)
  # row names
  row.names(des_table) <- c("Deforestation (ha)",
                            "Price signal ($/tCPO)", 
                            "Public ownership (%)", 
                            "Domestic private ownership (%)", 
                            "Foreign ownership (%)", 
                            #"Competition", 
                            "# reachable mills")
  
  return(des_table)
}


### OVERALL 
res_data_list_des <- list()
res_data_list_des[["both_a"]] <- make_base_reg(island = "both",
                                               outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                               output_full = TRUE)

# # export for analysis on other plateform
# d_clean_sf <- st_as_sf(d_clean, coords = c("lon","lat"), crs = 4326)
# st_write(d_clean_sf, file.path("temp_data/actual_analysis_data_Kalimantan_a"), driver = "ESRI Shapefile", delete_dsn = TRUE)



des_table <- make_desstats(sample_1 = res_data_list_des[["both_a"]][[2]],
                           sample_2 = res_data_list_des[["both_a"]][[3]])
des_table
colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Estimation sample - descriptive statistics") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "p-value" = 1, "p-value" = 1),
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "# grid cells = 4757 \n # grid cell-year = 31721" = 3,
                     "# grid cells = 8309 \n # grid cell-year = 87004" = 3, 
                     " " = 2),
                   align = "c",
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Without missing values" = 3,
                     "With missing values" = 3, 
                     "t test" = 1,
                     "KS test" = 1), 
                   bold = TRUE,
                   align = "c",
                   strikeout = F) 
# pack_rows("Ownership signals (%)", 3, 5, 
#           italic = TRUE)  %>%
# column_spec(column = c(2,3,5,6,8,9),
#             width = "3em") %>%
# column_spec(column = c(4,7,10),
#             width = "9em") %>%



#### Table of deforestation in different catchment radius / sample ####

### TOTAL DEFORESTATION 
rows <- c("Sumatra", "Kalimantan", "both")
cols <- c("sample", "30km", "50km", "total")

accu_lucpfp <- matrix(ncol = length(cols), nrow = length(rows)) # 4 cols for sample, 30km, 50km and total. 3 rows for Sumatra, Kalimantan, and both. 
row.names(accu_lucpfp) <- rows
colnames(accu_lucpfp) <- cols



#island <- "Sumatra"
for(island in c("Sumatra", "Kalimantan")){
  # LUCPFIP
  defo_dyna <- list()
  for(DYNA in c("rapid", "slow")){
    brick_lucpfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",DYNA,"_",island,"_",parcel_size/1000,"km_total.tif")))
    
    # select 2002-2014 layers 
    layer_names <- paste0("parcel_lucpfip_",DYNA,"_",island,"_",parcel_size/1000,"km_total.",c(2:14))
    
    brick_lucpfip <- raster::subset(brick_lucpfip, layer_names)
    
    # Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
    r_accu_lucfp <- calc(brick_lucpfip, fun = sum, na.rm = TRUE)
    
    
    
    
    defo_dyna[[match(DYNA, c("rapid", "slow"))]] <- raster::extract(x = r_accu_lucfp, 
                                                                    y = island_sf_prj[island_sf_prj$shape_des==island,] %>% st_geometry() %>% as("Spatial"), 
                                                                    fun = sum, 
                                                                    na.rm = TRUE) 
  }
  accu_lucpfp[island, "total"] <- sum(unlist(defo_dyna))
  
  # LUCPSMP
  defo_sm <- list()
  for(size in c("s", "m")){
    brick_lucpfsmp <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpf",size,"p_",island,"_",parcel_size/1000,"km_total.tif")))
    
    # select 2002-2014 layers 
    layer_names <- paste0("parcel_lucpf",size,"p_",island,"_",parcel_size/1000,"km_total.",c(2:14))
    
    brick_lucpfsmp <- raster::subset(brick_lucpfsmp, layer_names)
    
    # Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
    r_accu_lucfp <- calc(brick_lucpfsmp, fun = sum, na.rm = TRUE)
    
    defo_sm[[match(size, c("s", "m"))]] <- raster::extract(x = r_accu_lucfp, 
                                                           y = island_sf_prj[island_sf_prj$shape_des==island,] %>% st_geometry() %>% as("Spatial"), 
                                                           fun = sum, 
                                                           na.rm = TRUE) 
  }
  
  # Add indus and smallholder deforestation
  accu_lucpfp[island, "total"] <- sum(unlist(defo_dyna)) + sum(unlist(defo_sm)) 
}

### Now restrict to full CR. 

# sum pixelcounts over parcels and years, under the condition that n_reachable is not null
accu_lucpfp["Sumatra", "30km"] <- dplyr::filter(d_30_suma, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Sumatra", "50km"] <- dplyr::filter(d_50_suma, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Kalimantan", "30km"] <- dplyr::filter(d_30_kali, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Kalimantan", "50km"] <- dplyr::filter(d_50_kali, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()

### FOr analysis sample
for(ISL in c("Sumatra", "Kalimantan")){
  res_data_list  <- make_base_reg(island = ISL,
                                  outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                  output_full = FALSE)
  
  d_clean <- res_data_list[[2]]
  accu_lucpfp[ISL, "sample"] <- sum(d_clean$lucpfap_pixelcount)
}


# Fill the row of total over islands
for(col in c(1:ncol(accu_lucpfp))){
  accu_lucpfp["both", col] <- sum(accu_lucpfp["Sumatra", col], 
                                  accu_lucpfp["Kalimantan", col])
}

accu_lucpfp <- accu_lucpfp*pixel_area_ha/1000 # convert to *kilo* hectares
accu_lucpfp <- accu_lucpfp %>% round(digits = 2)

row.names(accu_lucpfp)[3] <- "Both"

accu_lucpfp_save <- accu_lucpfp

# Look at how much land is concerned with immediate and transitional deforestation
res_data_list  <- make_base_reg(island = "both",
                                #end_year = 2010,
                                outcome_variable = paste0("lucpfip_rapid_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                output_full = FALSE)

d_clean <- res_data_list[[2]]
sum(d_clean$lucpfip_rapid_pixelcount)*pixel_area_ha/1000 %>% round(digits = 2)

res_data_list  <- make_base_reg(island = "both",
                                end_year = 2010,
                                outcome_variable = paste0("lucpfip_slow_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                output_full = FALSE)

d_clean <- res_data_list[[2]]
sum(d_clean$lucpfip_slow_pixelcount)*pixel_area_ha/1000 %>% round(digits = 2)
rm(d_clean)

# Finally, add figures from Austin et al. 2017 (SI) - no because not comparable enough because we add smallholders without breaking down here. 
# austin <- matrix(nrow = 2, ncol = ncol(accu_lucpfp))
# austin["Sumatra", 4] <- "(448)"
# austin["Kalimantan", 4] <- "(1021)"

### Print the LateX table code ###

accu_lucpfp
colnames(accu_lucpfp) <- NULL

options(knitr.table.format = "latex") 
kable(accu_lucpfp, booktabs = T, align = "c", 
      caption = "Deforestation accumulated over 2002-2014, in kha.") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "Sample" = 1, 
                     "30km from sample mill" = 1, 
                     "50km from sample mill" = 1,
                     "Total" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  column_spec(column = 1,
              width = "5em") %>%
  column_spec(column = c(2:ncol(accu_lucpfp)),
              width = "5em") 


### INDUSTRIAL VS SMALLHOLDERS
# res_data_list_des <- list()
# for(SIZE in c("i", "sm")){
#   res_data_list_des[[paste0("both_",SIZE)]] <- make_base_reg(island = "both",
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                                                  output_full = FALSE)
# }
# # export for analysis on other plateform
# d_clean_sf <- st_as_sf(d_clean, coords = c("lon","lat"), crs = 4326)
# st_write(d_clean_sf, file.path("temp_data/actual_analysis_data_Kalimantan_a"), driver = "ESRI Shapefile", delete_dsn = TRUE)



# des_table <- make_desstats(sample_1 = res_data_list_des[["both_i"]][[2]],
#                            sample_2 = res_data_list_des[["both_sm"]][[2]])
# colnames(des_table) <- NULL
# 
# options(knitr.table.format = "latex") 
# kable(des_table, booktabs = T, align = "c", 
#       caption = "Estimation sample - descriptive statistics") %>% 
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
#   add_header_above(c(" " = 1, 
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
#                      "p-value" = 1, "p-value" = 1),
#                    align = "c", 
#                    strikeout = F) %>% 
#   add_header_above(c(" " = 1, 
#                      "# grid cells = 4757 \n # grid cell-year = 31721" = 3,
#                      "# grid cells = 8368 \n # grid cell-year = 91620" = 3, 
#                      " " = 2),
#                    align = "c",
#                    strikeout = F) %>% 
#   add_header_above(c(" " = 1,
#                      "Without missing values" = 3,
#                      "With missing values" = 3, 
#                      "t test" = 1,
#                      "KS test" = 1), 
#                    bold = TRUE,
#                    align = "c",
#                    strikeout = F) %>% 
#   # pack_rows("Ownership signals (%)", 3, 5, 
#   #           italic = TRUE)  %>%
#   # column_spec(column = c(2,3,5,6,8,9),
#   #             width = "3em") %>%
#   # column_spec(column = c(4,7,10),
#   #             width = "9em") %>%
#   
#   footnote(general = c("Note"),
#            threeparttable = TRUE, 
#            escape = TRUE) 


### Descriptive IBS ###

# all IBS oil palm related establishments. 
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))

length(unique(ibs$firm_id))
length(unique(ibs$firm_id[ibs$uml_matched_sample==1])) 
length(unique(ibs$firm_id[ibs$is_mill==1])) 
length(unique(ibs$firm_id[ibs$analysis_sample==1]))
# ibs[ibs$geo_sample & !ibs$is_mill,1:40]
# out of 1473 establishments initially extracted from IBS,
# 1004 are involved at least one year in some FFB processing or CPO or PKO producing, but 930 are not in Java nor Bali 
# out of which 468 have been geolocalized (and 2 more are in Java)
# 2 additional IBS firms have been matched with UML and geolocalized but have no sign of FFB, PKO or CPO activity. 
# (firm_id 2036 and 55630)


# Now we also remove those that have some FFB-CPO-PKO activity but have been identified as refineries. 
# But why would we remove refineries? They may be different but if they input FFB, they have an impact on proximate LUC. 
# Because we excluded refineries from geo_sample, and here the purpose is to compare this sample to the population of mills. 
# let's first see the comparative stat des without removing refineries. 



# We want to produce, for a set of ibs variables, the mean, median, std.dev. min and max statistics, 
# for the geo_sample that we are going to use, and the is_mill sample. 

make_des_table_ibs <- function(ibs_isl){
  
  variables <- c("min_year", 
                 "ffb_price_imp1", "in_ton_ffb_imp1",
                 "cpo_price_imp1", "out_ton_cpo_imp1",
                 "pko_price_imp1", "out_ton_pko_imp1",
                 "prex_cpo_imp1", 
                 "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp")
  
  statistics <- c("mean", "std.dev", "median", "min", "max")
  
  
  ## Matrix for analysis_sample
  rhs_des_sample <- matrix(NA, nrow = length(variables), ncol = length(statistics))
  row.names(rhs_des_sample) <- variables
  colnames(rhs_des_sample) <- c(statistics)
  
  for(var in variables){
    rhs_des_sample[var,statistics] <- summarise(dplyr::filter(ibs_isl, analysis_sample == TRUE),
                                                mean = mean(get(var), na.rm=T),
                                                std.dev = sd(get(var), na.rm= TRUE),
                                                median = median(get(var), na.rm= TRUE), 
                                                min = min(get(var), na.rm= TRUE),
                                                max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    rhs_des_sample[var, "median"] <- paste0(rhs_des_sample[var,"median"]," [",rhs_des_sample[var, "min"],"; ",rhs_des_sample[var,"max"],"]")
  }
  
  rhs_des_sample <- rhs_des_sample[,c("mean", "std.dev", "median")]
  
  # number of establishments in sub-group
  # N_sample_row <- c(length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"])), 
  #                           rep("", length(statistics))) 
  # rhs_des_sample <- rbind(N_sample_row, rhs_des_sample)
  
  ## Matrix for all ibs mills
  rhs_des_pop <- matrix(NA, nrow = length(variables), ncol = length(statistics))
  row.names(rhs_des_pop) <- variables
  colnames(rhs_des_pop) <- c(statistics)
  
  for(var in variables){
    rhs_des_pop[var,statistics] <- summarise(dplyr::filter(ibs_isl, is_mill == TRUE),
                                             mean = mean(get(var), na.rm=T),
                                             std.dev = sd(get(var), na.rm= TRUE),
                                             median = median(get(var), na.rm= TRUE), 
                                             min = min(get(var), na.rm= TRUE),
                                             max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    rhs_des_pop[var, "median"] <- paste0(rhs_des_pop[var,"median"]," [",rhs_des_pop[var, "min"],"; ",rhs_des_pop[var,"max"],"]")
  }
  
  rhs_des_pop <- rhs_des_pop[,c("mean", "std.dev", "median")]
  # # number of establishments in sub-group
  # N_pop_row <- c(length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"])), 
  #                   rep("", length(statistics))) 
  # rhs_des_pop <- rbind(N_pop_row, rhs_des_pop)
  
  # bind two groups together
  rhs_des <- cbind(rhs_des_sample, rhs_des_pop)    
  
  # add the t-test column
  t_test <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(t_test) <- variables
  colnames(t_test) <- "t test"
  
  for(var in variables){
    test <- t.test(x = ibs_isl[ibs_isl$analysis_sample==TRUE,var],
                   y = ibs_isl[ibs_isl$is_mill==TRUE,var], 
                   alternative = "two.sided", 
                   mu = 0, 
                   paired = F, 
                   var.equal = FALSE)
    
    t_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  # interpretation: if the 95% CI includes 0, then the difference in means between two samples 
  # is not statistically different from 0. Hence, the two samples are "similar" in means 
  # with respect to the variable tested. 
  # (In other words, we cannot reject the null hypothesis that the difference is null
  # -i.e. the two samples are alike - with 95% confidence) 
  
  # add the Smirnov test
  ks_test <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(ks_test) <- variables
  colnames(ks_test) <- "KS test"
  
  for(var in variables){
    test <- ks.test(x = ibs_isl[ibs_isl$analysis_sample==TRUE,var],
                    y = ibs_isl[ibs_isl$is_mill==TRUE,var], 
                    alternative = "two.sided", 
                    exact = FALSE)
    
    ks_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  
  # intepretation: if p-value < 0.05 we cannot reject with 95% confidence that 
  # two distributions are equal, implying that they are different. 
  
  rhs_des <- cbind(rhs_des, t_test, ks_test)
  
  
  
  # row names
  row.names(rhs_des) <- c("First year in IBS", "FFB muv (USD/ton)", "FFB input (ton)", 
                          "CPO muv (USD/ton", "CPO output (ton)", 
                          "PKO muv (USD/ton)", "PKO output (ton)", 
                          "CPO export share (%)", 
                          "Central government ownership (%)", 
                          "Local government ownership (%)", 
                          "National private ownership (%)", 
                          "Foreign ownership (%)")
  
  return(rhs_des)
}

#### Print the LateX table code ALL #### 
des_table <- make_des_table_ibs(ibs) 
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, all islands") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "p-value" = 1,
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 587 mills" = 3,
                     "All IBS palm oil mills \n n = 930 mills" = 3, 
                     "t-test" = 1,
                     "KS test" = 1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "4em") %>% 
  column_spec(column = c(4,7),
              width = "9em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 


#### VARIATIONS IN PRICES WITHIN DISTRICT YEAR #### 
# Restrict to analysis sample 
as <- ibs[ibs$analysis_sample==TRUE,]
# remove district-year variations 
rm_fevar <- fixest::feols(fml = as.formula("cpo_price_imp1 ~ 1 | district^year"),
                          data = as)

# and take the standard deviation of remaining variation
sd(rm_fevar$residuals)

#### DESCRIPTIVE MAP #####
res_data <- make_base_reg(island = "both",
                          outcome_variable = paste0("lucpfap_pixelcount"))

d_clean <- res_data[[2]]
d_clean_cs <- d_clean[!duplicated(d_clean$lonlat),]

d_cs <- st_as_sf(d_clean_cs, coords = c("lon", "lat"), crs = 4326)
d_cs <- st_transform(d_cs, crs = indonesian_crs)

d_cs <- st_buffer(d_cs, dist = 1500)
st_geometry(d_cs) <- sapply(st_geometry(d_cs), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
d_geo <- st_union(st_geometry(d_cs))
d_geo <- st_transform(d_geo, crs = 4326)


# d_cs <- ddply(d_clean, "lonlat", summarise, 
#               accu_lucfp = sum(lucpfap_pixelcount))
# 
# d_cs <- left_join(d_cs, d_clean_cs[,c("lonlat", "lon", "lat")], by = "lonlat")





### MILLs
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# keep only those that are used in analysis. 
ibs <- ibs[ibs$analysis_sample ==1,]
# make it a cross section
ibs <- ibs[!duplicated(ibs$firm_id),]
ibs <- ibs[!is.na(ibs$lat),]
ibs <- st_as_sf(ibs, coords = c("lon", "lat"), remove = FALSE, crs = 4326)

# we do not plot uml actually, it's too much. 
# uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
# uml <- uml[!is.na(uml$lat),]
# uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = FALSE, crs = 4326)
# # remove those matched with ibs
# uml <- uml[!(uml$trase_code %in% ibs$trase_code),]

# prepare backgroud layers with other countries
countries <- st_read(file.path("input_data/Global_LSIB_Polygons_Detailed"))
countries <- countries[countries$COUNTRY_NA == "Indonesia" | 
                         countries$COUNTRY_NA == "Malaysia" | 
                         countries$COUNTRY_NA == "Brunei", "geometry"]
# these two lines to speed up mapping
countries <- st_transform(countries, crs = indonesian_crs) %>% st_simplify(dTolerance = 1000)
countries <- st_transform(countries, crs = 4326)

ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_sf(data = countries, fill = FALSE) +
  geom_sf(data = ibs, color = "black", size = 0.05) + 
  geom_sf(data=st_geometry(d_geo), color=alpha("grey",0.2))+
  coord_sf(xlim = c(90, 120), ylim = c(-7, 7), expand = FALSE) 

# leaflet() %>% 
#   # we do not add these tiles if the figure is purposed for print in white and black. 
#   # addTiles()%>%
#   # addProviderTiles(providers$Esri.WorldGrayCanvas) %>% # CartoDB.PositronNoLabels 
#   addPolygons(data = countries, stroke = TRUE, color = "black", weight = 1,
#               fill = FALSE) %>%
#   setView(lat = -0.493, 
#           lng = 107.367, 
#           zoom = 4.5) %>% 
#   #addPolygons(data = island_sf, stroke = FALSE, fill = FALSE, fillColor = "grey", fillOpacity = 0) %>% 
#   addPolygons(data = st_geometry(d_geo), 
#               stroke = TRUE, opacity = 0.25,  weight = 2, color = "green") %>%   #~cb_ind(d_cs$accu_lucfp)
#   addCircleMarkers(data = ibs, radius = 0.001, fillOpacity = 1, fillColor = "black", stroke = FALSE, weight = 0) 
# # exported in width = 1150 and height = 560 and zoom once 


rm(d_clean_cs, d_cs, d_geo, ibs)

#### DESCRIPTIVE PARTIAL AUTOCORRELATION FUNCTION OF PRICES #### 
# read this panel, as it still features the annual price observations from 1998 (final data only from 2001)
# 50km CR because it is the most general
RHS_50 <-  readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                    parcel_size/1000,"km_",
                                    "50CR.rds")))

prices <- RHS_50[,c("lonlat", "year", "cpo_price_imp1")]

nrow(prices)/length(unique(prices$year))

# Apply descriptive statistics to a sample of plantations, because neighboring plantations can have spatially correlated price time series, 
# that would conflate the inference performed here.
sample_ids <- prices$lonlat %>% sample(size = 700) # 700 this is roughly 1% of the grid cells. 
prices <- prices[prices$lonlat %in% sample_ids,]

fp <- prices[,c("lonlat", "year")]

# separate each time series (of each plantation) by an "empty" (NA) time series, to isolate them from each others. 
length_panel <- length(unique(prices$year))
fp <- mutate(fp, year = year + length_panel)
fp$cpo_price_imp1 <- NA
prices_extd <- rbind(prices, fp)
prices_extd <- dplyr::arrange(prices_extd, lonlat, year)
prices_ts <- prices_extd[,c("cpo_price_imp1")]
prices_ts <- ts(prices_ts, frequency = length_panel*2)#

ur.kpss(prices_ts) %>% summary() # --> data need to be differenced
ur.kpss(diff(prices_ts, differences = 1)) %>% summary() # --> first differencing suffices

arima(prices_ts, order = c(length_panel - 1,1,0)) 


# this does not handle properly the specific structure of our time series (with NA sequences separating grid cells' respective time series)
# pacf_14 <- pacf(diff(prices_ts, differences = 1), na.action = na.exclude, lag.max = length_panel - 1)

rm(RHS_50)

#### MAIN : REGRESSIONS DES STATS AND EQUALITY TESTS OVER INDUS, SM, LEGAL AND ILLEGAL,  ####

### REGRESSIONS ### 
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
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Price elasticities of deforestation across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3, 
                     "All" = 3),
                   align = "c",
                   strikeout = F) %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

rm(ape_mat)


### PARTIAL EFFECTS WITH ALL CONTROLS 
ape_mat <- bind_cols(lapply(res_data_list_full, FUN = make_APEs, controls_pe = TRUE)) %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Price elasticity  and partial effects of control variables on deforestation across Indonesian oil palm sectors") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3, 
                     "All" = 3),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("Price elasticity", 1, 2, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Partial effects of:", 3, 8,
            italic = FALSE, bold = TRUE) %>%
  pack_rows("Domestic private \n mill ownership", 3, 4,
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Foreign mill \n ownership", 5, 6,
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("# reachable mills", 7, 8,
            italic = TRUE, bold = TRUE)  %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = ncol(ape_mat)+1,
              bold = TRUE)


### DESCRIPTIVE STATISTICS ACROSS PLANTATION TYPES ---------------------------------

# We make three tables, for readability purpose; one for industrial, one for smallholder, and one for all plantations. 
# Each one has the legal/illegal/all breakdown. Each features same quantities as feature in the main des stat table above. 
make_desstats_simple <- function(sample_1){
  
  # add variables that were not used in the regressions
  d_all <- rbind(d_30_suma, d_50_kali)
  # determine dependent variable
  dep_var_1 <- names(sample_1)[grepl("pixelcount",names(sample_1))]
  
  # group public ownership
  #d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  
  variables_1 <- c(dep_var_1,
                   "wa_cpo_price_imp1_4ya_lag1",
                   "wa_pct_own_gov_imp_lag1", "wa_pct_own_nat_priv_imp_lag1", "wa_pct_own_for_imp_lag1",
                   "n_reachable_uml") 
  sample_1 <- left_join(sample_1[,c("lonlat","year")], d_all[,c("lonlat","year",variables_1)], by = c("lonlat","year"))#
  
  # pixelcount to hectares
  sample_1[,dep_var_1] <- sample_1[,dep_var_1]*pixel_area_ha
  
  statistics <- c("mean", "std.dev", "median", "min", "max")
  
  ## Des. stats. for sample_1
  des_sample_1 <- matrix(NA, nrow = length(variables_1), ncol = length(statistics))
  row.names(des_sample_1) <- variables_1
  colnames(des_sample_1) <- c(statistics)
  
  for(var in variables_1){
    des_sample_1[var,statistics] <- summarise(sample_1,
                                              mean = mean(get(var), na.rm=TRUE),
                                              std.dev = sd(get(var), na.rm= TRUE),
                                              median = median(get(var), na.rm= TRUE), 
                                              min = min(get(var), na.rm= TRUE),
                                              max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    des_sample_1[var, "median"] <- paste0(des_sample_1[var,"median"]," [",des_sample_1[var, "min"],"; ",des_sample_1[var,"max"],"]")
  }
  
  des_sample_1 <- des_sample_1[,c("mean", "std.dev", "median")]
  
  length(unique(sample_1$lonlat)) %>% paste0(" number of grid cells in sample") %>%  print()
  nrow(sample_1) %>% paste0(" number of observations in sample") %>%  print()
  
  
  
  
  return(des_sample_1)
}


# ## Descriptive statistics across legal/illegal, for all plantation types
# # template to store 
# list_desstat_all <- list()
# i <- 1
# for(ILL in ill_status){
#   list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[[paste0("both_a_",ILL)]][[2]])
#   i <- i +1
# }
# 
# des_table <- bind_cols(list_desstat_all) %>% as.matrix()
# # row names
# row.names(des_table) <- c("Deforestation (ha)",
#                              "Price signal ($/tCPO)", 
#                              "Public ownership (%)", 
#                              "Domestic private ownership (%)", 
#                              "Foreign ownership (%)", 
#                              #"Competition", 
#                              "# reachable mills")
# des_table
# colnames(des_table) <- NULL
# 
# options(knitr.table.format = "latex") 
# kable(des_table, booktabs = T, align = "c", 
#       caption = "Estimation sample for all plantation types - descriptive statistics") %>% 
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
#   add_header_above(c(" " = 1, 
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1),
#                    align = "c", 
#                    strikeout = F) %>% 
#   add_header_above(c(" " = 1, 
#                      "# grid cells = 2245 \n # grid cell-year = 15139" = 3,
#                      "# grid cells = 1216 \n # grid cell-year = 7848" = 3, 
#                      "# grid cells = 4747 \n # grid cell-year = 31650" = 3, 
#                      " " = 2),
#                    align = "c",
#                    strikeout = F) %>% 
#   add_header_above(c(" " = 1,
#                      "Legal" = 3,
#                      "Illegal" = 3, 
#                      "All" = 3), 
#                    bold = FALSE,
#                    align = "c",
#                    strikeout = F) 

## Descriptive statistics across legal/illegal, for industrial plantations
# template to store 
list_desstat_all <- list()
i <- 1
for(ILL in ill_status){
  list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[[paste0("both_i_",ILL)]][[2]])
  i <- i +1
}

des_table <- bind_cols(list_desstat_all) %>% as.matrix()
# row names
row.names(des_table) <- c("Deforestation (ha)",
                          "Price signal ($/tCPO)", 
                          "Public ownership (%)", 
                          "Domestic private ownership (%)", 
                          "Foreign ownership (%)", 
                          #"Competition", 
                          "# reachable mills")
des_table
colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Estimation sample for industrial plantations - descriptive statistics") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1),
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "# grid cells = 1979 \n # grid cell-year = 13081" = 3,
                     "# grid cells = 781 \n # grid cell-year = 4951" = 3, 
                     "# grid cells = 3853 \n # grid cell-year = 25249" = 3, 
                     " " = 2),
                   align = "c",
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Legal" = 3,
                     "Illegal" = 3, 
                     "All" = 3), 
                   bold = FALSE,
                   align = "c",
                   strikeout = F) 

## Descriptive statistics across legal/illegal, for smallholder plantations
# template to store 
list_desstat_all <- list()
i <- 1
for(ILL in ill_status){
  list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[[paste0("both_sm_",ILL)]][[2]])
  i <- i +1
}

des_table <- bind_cols(list_desstat_all) %>% as.matrix()
# row names
row.names(des_table) <- c("Deforestation (ha)",
                          "Price signal ($/tCPO)", 
                          "Public ownership (%)", 
                          "Domestic private ownership (%)", 
                          "Foreign ownership (%)", 
                          #"Competition", 
                          "# reachable mills")
des_table
colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Estimation sample for smallholder plantations - descriptive statistics") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1),
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "# grid cells = 385 \n # grid cell-year = 2971" = 3,
                     "# grid cells = 522 \n # grid cell-year = 3412" = 3, 
                     "# grid cells = 1189 \n # grid cell-year = 8611" = 3, 
                     " " = 2),
                   align = "c",
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Legal" = 3,
                     "Illegal" = 3, 
                     "All" = 3), 
                   bold = FALSE,
                   align = "c",
                   strikeout = F) 


### SPATIAL BREAKDOWN ----------------------------------------------------
ape_mat_list <- list()
ape_elm <- 1

isl_list <- list("Sumatra", "Kalimantan") 
for(ISL in isl_list){
  res_data_list_bd <- list()
  elm <- 1
  
  
  size_list <- list("i","sm", "a")
  
  # legality definition
  ill_def <- 2
  ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")
  
  
  for(SIZE in size_list){
    for(ILL in ill_status){
      if(!(SIZE == "sm" & ISL == "Kalimantan" & ILL == "ill2")){# because in this case there is not enough variation
        
        res_data_list_bd[[elm]] <- make_base_reg(island = ISL,
                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                                 illegal = ILL,
                                                 offset = FALSE)
        names(res_data_list_bd)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
        elm <- elm + 1
      }
    }
  }
  
  ## PARTIAL EFFECTS
  rm(ape_mat)
  ape_mat <- bind_cols(lapply(res_data_list_bd, FUN = make_APEs)) %>% as.matrix()
  
  row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
  ape_mat
  colnames(ape_mat) <- NULL
  
  if(ISL == "Kalimantan"){
    ape_mat <- cbind(ape_mat[,1:4], rep("", 4), ape_mat[,5:8])} # fill the column so that it aligns with Sumatra
  
  ape_mat_list[[ape_elm]] <- ape_mat
  ape_elm <- ape_elm + 1
}

stacked_ape_mat <- rbind(ape_mat_list[[1]], ape_mat_list[[2]])

options(knitr.table.format = "latex")
kable(stacked_ape_mat, booktabs = T, align = "r",
      caption = "Price elasticities of deforestation across the oil palm sector, by island") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3, 
                     "All" = 3),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("Sumatra", 1, 4,
            italic = TRUE, bold = TRUE) %>% 
  pack_rows("Kalimantan", 5, 8,
            italic = TRUE, bold = TRUE) %>% 
  # pack_rows("Price elasticity", 1, 2, 
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows(start_row =  nrow(stacked_ape_mat)-1, end_row = nrow(stacked_ape_mat),  latex_gap_space = "1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(stacked_ape_mat))),
              width = "7em",
              latex_valign = "b") 



### TEMPORAL BREAKDOWN -----------------------------------------------------------
# ISL <- "both"
# size_list <- list("i","sm", "a")
# 
# # legality definition
# ill_def <- 2
# ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")
# 
# res_data_list_bd <- list()
# elm <- 1
# for(SIZE in size_list){
#   for(ILL in ill_status){
#     res_data_list_bd[[elm]] <- make_base_reg(island = ISL,
#                                              start_year = 2002,
#                                              end_year = 2008,
#                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                                              illegal = ILL,
#                                              offset = FALSE)
#     names(res_data_list_bd)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#     elm <- elm + 1
#   }
# }
# ## PARTIAL EFFECTS
# rm(ape_mat)
# ape_mat <- bind_cols(lapply(res_data_list_bd, FUN = make_APEs)) %>% as.matrix()
# 
# row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters")
# ape_mat
# colnames(ape_mat) <- NULL
# 
# ape_mat_list[[ape_elm]] <- ape_mat
# ape_elm <- ape_elm + 1
# 
# # repeat for 2009-2014
# res_data_list_bd <- list()
# elm <- 1
# for(SIZE in size_list){
#   for(ILL in ill_status){
#     res_data_list_bd[[elm]] <- make_base_reg(island = ISL,
#                                              start_year = 2009,
#                                              end_year = 2014,
#                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                                              illegal = ILL,
#                                              offset = FALSE)
#     names(res_data_list_bd)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#     elm <- elm + 1
#   }
# }
# ## PARTIAL EFFECTS
# rm(ape_mat)
# ape_mat <- bind_cols(lapply(res_data_list_bd, FUN = make_APEs)) %>% as.matrix()
# 
# row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters")
# ape_mat
# colnames(ape_mat) <- NULL
# 
# ape_mat_list[[ape_elm]] <- ape_mat
# ape_elm <- ape_elm + 1
# 
# 
# stacked_ape_mat <- rbind(ape_mat_list[[1]], ape_mat_list[[2]])
# 
# options(knitr.table.format = "latex")
# kable(stacked_ape_mat, booktabs = T, align = "r",
#       caption = "Price elasticities of deforestation across the Indonesian oil palm sector, by time period") %>% #of 1 percentage change in medium-run price signal
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c(" " = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1),
#                    bold = F,
#                    align = "c") %>%
#   add_header_above(c(" " = 1,
#                      "Industrial plantations" = 3,
#                      "Smallholder plantations" = 3,
#                      "All" = 3),
#                    align = "c",
#                    strikeout = F) %>%
#   pack_rows("2002-2008", 1, 4,
#             italic = TRUE, bold = TRUE) %>%
#   pack_rows("2009-2014", 5, 8,
#             italic = TRUE, bold = TRUE) %>%
#   column_spec(column = 1,
#               width = "7em",
#               latex_valign = "b") %>%
#   column_spec(column = c(2:(ncol(stacked_ape_mat))),
#               width = "7em",
#               latex_valign = "b") %>%
#   column_spec(column = ncol(stacked_ape_mat)+1,
#               bold = TRUE)

# rm(ape_mat)

### SECONDARY FORESTS -----------------------

# infrastructure to store results
res_data_list_full_2ndry <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")


for(SIZE in size_list){
  for(ILL in ill_status){
    res_data_list_full_2ndry[[elm]] <- make_base_reg(island = ISL,
                                                     outcome_variable = paste0("lucf",SIZE,"p_pixelcount"), # note the absence of "p": it is not primary forest data. 
                                                     illegal = ILL,
                                                     offset = FALSE)
    names(res_data_list_full_2ndry)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }
}

## PARTIAL EFFECTS
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list_full_2ndry, FUN = make_APEs)) %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Price elasticities of deforestation, in secondary forest, across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3, 
                     "All" = 3),
                   align = "c",
                   strikeout = F) %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

rm(ape_mat)


### COMPARE GROUPS --------------------------------------------------------------------------------

compare_APEs_across_groups <- function(group1, group2, m0 = 0, alternative = "two.sided") { 
  # important that it's selecting the "1" (resp. 2) row, and not the "Estimate" (resp. "SE") row in cases where there are several 
  # rows with the same names, i.e. if APEs of interaction effects have been computed as well. 
  ape1 <- ape_mat[1,group1] 
  ape2 <- ape_mat[1,group2]
  # # n1 <- ape_mat["Observations","Sumatra industrial"]
  # # n2 <- ape_mat["Observations","Sumatra smallholders"]
  sigma1 <- ape_mat[2,group1]^2
  sigma2 <- ape_mat[2,group2]^2
  
  statistic <- (ape1 - ape2 - m0) / sqrt(sigma1 + sigma2)
  
  pval <- if (alternative == "two.sided") { 
    2 * pnorm(abs(statistic), lower.tail = FALSE) 
  } else if (alternative == "less") { 
    pnorm(statistic, lower.tail = TRUE) 
  } else { 
    pnorm(statistic, lower.tail = FALSE) 
  } 
  # LCL <- (M1 - M2 - S * qnorm(1 - alpha / 2)) UCL <- (M1 - M2 + S * qnorm(1 - alpha / 2)) value <- list(mean1 = M1, mean2 = M2, m0 = m0, sigma1 = sigma1, sigma2 = sigma2, S = S, statistic = statistic, p.value = p, LCL = LCL, UCL = UCL, alternative = alternative) 
  # print(sprintf("P-value = %g",p)) # print(sprintf("Lower %.2f%% Confidence Limit = %g", 
  # alpha, LCL)) # print(sprintf("Upper %.2f%% Confidence Limit = %g", # alpha, UCL)) return(value) } test <- t.test_knownvar(dat1$sample1, dat1$sample2, V1 = 1, V2 = 1 )
  return(pval)
}


## infrastructure to store results
ape_mat_list <- list()
elm <- 1

groups <- c("both_a_no_ill2", "both_a_ill2", "both_a_all",
            "both_i_all",
            "both_sm_all")

comparisons <- c("industrial = smallholders",
                 "legal = illegal")#,"primary forest = broadly def. forest"

comp_ape_mat <- matrix(ncol = length(groups), nrow = length(comparisons), data = NA)
colnames(comp_ape_mat) <- groups
row.names(comp_ape_mat) <- comparisons

ISL <- "both"
for(SIZE in size_list){
  for(ILL in ill_status){
    
    # make the APE
    ape_mat_list[[elm]] <- make_APEs_1regr(res_data = res_data_list_full[[paste0(ISL,"_",SIZE, "_",ILL)]])
    names(ape_mat_list)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }
}

rm(ape_mat)
ape_mat <- bind_cols(ape_mat_list) %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","SE","p-value"), (nrow(ape_mat)-1)/3), "Observations")
# keep ape_mat like this for later comparisons between APEs

# fill the matrix of p values of equality tests BY ROW

for(ILL in ill_status){
  comp_ape_mat["industrial = smallholders",paste0("both_a_",ILL)] <- compare_APEs_across_groups(group1 = paste0("both_i_",ILL), 
                                                                                                group2 = paste0("both_sm_",ILL)) 
}

for(SIZE in size_list){
  comp_ape_mat["legal = illegal",paste0("both_",SIZE,"_all")] <- compare_APEs_across_groups(group1 = paste0("both_",SIZE,"_no_ill2"), 
                                                                                            group2 = paste0("both_",SIZE,"_ill2")) 
}

comp_ape_mat <- comp_ape_mat %>% formatC(digits = 4, format = "f")
comp_ape_mat[comp_ape_mat=="   NA"] <- "" 
comp_ape_mat
colnames(comp_ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(comp_ape_mat, booktabs = T, align = "r",
      caption = "p-values from equality tests of price elasticities") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Ho" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     " " = 1,
                     " " = 1),
                   bold = F,
                   align = "c") %>%  
  add_header_above(c(" " = 1,
                     "All plantations" = 3,
                     "Industrial plantations" = 1,
                     "Smallholder plantations" = 1),
                   bold = F,
                   align = "c") %>%
  column_spec(column = 1,
              width = "13em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(comp_ape_mat))),
              width = "5em",
              latex_valign = "m")




rm(res_data_list_full)


### WITH ALL INTERACTIONS -------------------------------------------------------------
res_data_list_interact <- list()
elm <- 1


# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  for(ILL in ill_status){
    res_data_list_interact[[elm]] <- make_base_reg(island = "both",
                                                   outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                                   interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml"),
                                                   illegal = ILL,
                                                   offset = FALSE)
    names(res_data_list_interact)[elm] <- paste0("both_",SIZE,"_",ILL)
    elm <- elm + 1
  }
}

## Partial effects
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list_interact, FUN = make_APEs, rounding = 4)) %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Price elasticity heterogeneity across ownership and local market development") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3, 
                     "All" = 3),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("Price signal", 1, 2, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction with ", 3, 8, 
            italic = FALSE, bold = TRUE)  %>%
  pack_rows("Domestic private \n ownership", 3, 4, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Foreign ownership", 5, 6, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("# reachable mills", 7, 8, 
            italic = TRUE, bold = TRUE)  %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "8em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "8em",
              latex_valign = "b")

rm(res_data_list_interact)

### LUC DYNAMICS --------------------------------------------------------------------------------------

## Reressions
# infrastructure to store results
res_data_list_dyn <- list()
elm <- 1

dyn_list <- list("rapid", "slow")

# first estimate rapid over the whole period 
for(ILL in ill_status){
  res_data_list_dyn[[elm]] <- make_base_reg(island = "both", 
                                            outcome_variable = paste0("lucpfip_rapid_pixelcount"),
                                            illegal = ILL,
                                            offset = FALSE)
  names(res_data_list_dyn)[elm] <- paste0("both_rapid_",ILL)
  elm <- elm + 1
}

# then rapid and slow in comparable times: i.e. up to 2010
for(DYN in dyn_list){
  for(ILL in ill_status){
    res_data_list_dyn[[elm]] <- make_base_reg(island = "both", 
                                              outcome_variable = paste0("lucpfip_",DYN,"_pixelcount"),
                                              illegal = ILL,
                                              end_year = 2010,
                                              offset = FALSE)
    names(res_data_list_dyn)[elm] <- paste0("both_",DYN,"_",ILL)
    elm <- elm + 1
  }
}

## Partial effects
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list_dyn, FUN = make_APEs)) %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Price elasticity of immediate and transitional deforestation in Indonesian industrial plantations") %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "Immediate conversion" = 3,
                     "Immediate conversion" = 3,
                     "Transitional conversion" = 3),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "2002 - 2010" = 3,
                     "2002 - 2014" = 6),
                   align = "c",
                   strikeout = F) %>%
  # pack_rows("Interaction with \n # of reachable mills", 4, 6, # domestic private ownership  "Interaction with \n domestic private ownership"
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n foreign ownership", 7, 9,
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n # of reachable mills", 10, 12,
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.2em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "4em",
              latex_valign = "c")


rm(res_data_list_dyn)



### PRICE DYNAMICS ---------------------------------------------------------------------------------

## Regressions

res_data_list_prdyn <- list()
elm <- 1

size_list <- list("i", "sm", "a")
# # legality definition
# ill_def <- 2
# ill_status <- c(paste0("ill",ill_def), "all")#paste0("no_ill",ill_def), do not include legal, to lighten tables



## With SR only 
for(SIZE in size_list){
  #for(ILL in ill_status){
  res_data_list_prdyn[[elm]] <- make_base_reg(island = "both",
                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                              only_sr = TRUE,
                                              offset = FALSE)
  names(res_data_list_prdyn)[elm] <- paste0("both_",SIZE,"_onlySR")
  elm <- elm + 1
}
# With SR AND MDR
for(SIZE in size_list){
  res_data_list_prdyn[[elm]] <- make_base_reg(island = "both",
                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                              dynamics = TRUE,
                                              interact_regressors = TRUE,
                                              offset = FALSE)
  names(res_data_list_prdyn)[elm] <- paste0("both_",SIZE,"_prdyn")
  elm <- elm + 1
  
}
## Partial effects
rm(ape_mat)
ape_mat1 <- bind_cols(lapply(res_data_list_prdyn[1:length(size_list)], FUN = make_APEs, K = 1)) %>% as.matrix()
ape_mat2 <- bind_cols(lapply(res_data_list_prdyn[(length(size_list)+1):(2*length(size_list))], FUN = make_APEs, K = 2, rounding = 3)) %>% as.matrix()

# add 4 lines to ape_mat1 (equivalent to estimate and CI of MR and interaction)
ape_mat1 <- rbind(ape_mat1,matrix(ncol = ncol(ape_mat1), nrow = 4))
ape_mat1 <- ape_mat1[c(1,2,5,6,7,8,3,4),]
ape_mat1[is.na(ape_mat1)] <- ""

ape_mat <- cbind(ape_mat1, ape_mat2)
ape_mat <- ape_mat[,c(1,4,2,5,3,6)]

# tant qu'on en est à bricoler
ape_mat[6, "both_a_prdyn"] <- "[0.003; 0.05]" # to check: rerun make_APEs function with rounding = 3
ape_mat[6, "both_sm_prdyn"] <- "[0.003; 0.08]"
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 

ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Partial effects of short- and medium-run price signals on deforestation across Indonesian oil palm plantations") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 2,
                     "Smallholder plantations" = 2, 
                     "All plantations" = 2),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("Short-run price signal", 1, 2, # short run is indeed always before MR in the regressor vector construction in make_base_reg
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Medium-run price signal", 3, 4, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction", 5, 6, 
            italic = TRUE, bold = TRUE)  %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "13em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "5em",
              latex_valign = "c")



rm(res_data_list_prdyn)

### COMMODITY -----------------------------------------------------------------------------------------
# infrastructure to store results
res_data_list_commo <- list()
elm <- 1

size_list <- list("i", "sm", "a")

## With FFB only 
for(SIZE in size_list){
  #for(ILL in ill_status){
  res_data_list_commo[[elm]] <- make_base_reg(island = "both",
                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                              commo = "ffb",
                                              offset = FALSE)
  names(res_data_list_commo)[elm] <- paste0("both_",SIZE,"_ffb")
  elm <- elm + 1
}
## With FFB *and* CPO prices  
for(SIZE in size_list){
  res_data_list_commo[[elm]] <- make_base_reg(island = "both",
                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                              commo = c("ffb","cpo"), # important that ffb is indeed in first position of the vector. 
                                              interact_regressors = TRUE,
                                              offset = FALSE)
  names(res_data_list_commo)[elm] <- paste0("both_",SIZE,"_ffbcpo")
  elm <- elm + 1
}

## Partial effects
rm(ape_mat)
ape_mat1 <- bind_cols(lapply(res_data_list_commo[1:length(size_list)], FUN = make_APEs, K = 1)) %>% as.matrix()
ape_mat2 <- bind_cols(lapply(res_data_list_commo[(length(size_list)+1):(2*length(size_list))], FUN = make_APEs, K = 2)) %>% as.matrix()

# add 4 lines to ape_mat1 (equivalent to estimate and CI of MR and interaction)
ape_mat1 <- rbind(ape_mat1,matrix(ncol = ncol(ape_mat1), nrow = 4))
ape_mat1 <- ape_mat1[c(1,2,5,6,7,8,3,4),]
ape_mat1[is.na(ape_mat1)] <- ""

ape_mat <- cbind(ape_mat1, ape_mat2)
ape_mat <- ape_mat[,c(1,4,2,5,3,6)]

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 

ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Partial effects of FFB and CPO prices on deforestation across Indonesian oil palm plantations") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 2,
                     "Smallholder plantations" = 2, 
                     "All plantations" = 2),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("FFB price signal", 1, 2, # short run is indeed always before MR in the regressor vector construction in make_base_reg
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("CPO price signal", 3, 4, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction", 5, 6, 
            italic = TRUE, bold = TRUE)  %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "9em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "5em",
              latex_valign = "c")




rm(res_data_list_commo)









### PRICE VARIABILITY ----------------------------------------------------
# infrastructure to store results
res_data_list_variability <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")


for(SIZE in size_list){
  for(ILL in ill_status){
    res_data_list_variability[[elm]] <- make_base_reg(island = ISL,
                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                                      illegal = ILL,
                                                      price_variation = TRUE,
                                                      offset = FALSE)
    names(res_data_list_variability)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }
}

## PARTIAL EFFECTS
rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list_variability, FUN = make_APEs)) %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "r",
      caption = "Effects of price variability on deforestation across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3, 
                     "All" = 3),
                   align = "c",
                   strikeout = F) %>%
  # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") 

rm(ape_mat)
##### SPECIFICATION CHARTS #####


# These charts compare coefficients from a single dependent variable across different models/specifications
# Therefore we produce one such chart for CPO TOTAL, FFB TOTAL (and possibly for LRs and SRs too)


# The function is written in base R and exactly copied from https://github.com/ArielOrtizBobea/spec_chart/blob/master/spec_chart_function.R
schart <- function(data, labels=NA, highlight=NA, n=1, index.est=1, index.se=2, index.ci=NA,
                   order="asis", ci=.95, ylim=NA, axes=T, heights=c(1,1), leftmargin=11, offset=c(0,0), ylab="Coefficient", lwd.border=1,
                   lwd.est=4, pch.est=21, lwd.symbol=2, ref=0, lwd.ref=1, lty.ref=2, col.ref="black", band.ref=NA, col.band.ref=NA,length=0,
                   col.est=c("grey60", "red3"), col.est2=c("grey80","lightcoral"), bg.est=c("white", "white"),
                   col.dot=c("grey60","grey95","grey95","red3"),
                   bg.dot=c("grey60","grey95","grey95","white"),
                   pch.dot=c(22,22,22,22), fonts=c(2,1), adj=c(1,1),cex=c(1,1)) {
  
  # Authors: Ariel Ortiz-Bobea (ao332@cornell.edu).
  # Version: March 10, 2020
  # If you like this function and use it, please send me a note. It might motivate
  # me to imporove it or write new ones to share.
  
  # Description of arguments
  
  # Data:
  # data: data.frame with data, ideally with columns 1-2 with coef and SE, then logical variables.
  # labels: list of labels by group. Can also be a character vector if no groups. Default is rownames of data.
  # index.est: numeric indicating position of the coefficient column.
  # index.se: numeric indicating position of the SE column.
  # index.ci: numeric vector indicating position of low-high bars for SE. Can take up to 2 CI, so vector can be up to length 4
  
  # Arrangement and basic setup:
  # highlight: numeric indicating position(s) of models (row) to highlight in original dataframe.
  # n: size of model grouping. n=1 removes groupings. A vector yields arbitrary groupings.
  # order: whether models should be sorted or not. Options: "asis", "increasing", "decreasing"
  # ci: numeric indicating level(s) of confidence. 2 values can be indicated.
  # ylim: if one wants to set an arbitrary range for Y-axis
  
  # Figure layout:
  # heights: Ratio of top/bottom panel. Default is c(1,1) for 1 50/50 split
  # leftmargin: amount of space on the left margin
  # offset: vector of numeric with offset for the group and specific labels
  # ylab: Label on the y-axis of top panel. Default is "Coefficient"
  # lwd.border: width of border and other lines
  
  # Line and symbol styles and colors:
  # lwd.est: numeric indicating the width of lines in the top panel
  # ref: numeric vector indicating horizontal reference lines(s) Default is 0.
  # lty.ref: Style of reference lines. Default is dash line (lty=2).
  # lwd.ref. Width of reference lines. Default is 1.
  # col.ref: vector of colors of reference lines. Default is black.
  # band.ref: vector of 2 numerics indicating upper abdn lower height for a band
  # col.band.ref: color of this band
  # col.est: vector of 2 colors indicating for "other" and "highlighted" models
  # col.est2: same for outer confidence interval if more than 1 confidence interval
  # col.dot: vector of 4 colors indicating colors for borders of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # bg.dot : vector of 4 colors indicating colors for background of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # pch.dot: style of symbols in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # length: length of the upper notch on th vertical lines. default is 0.
  
  # Letter styles
  # fonts: numeric vector indicating font type for group (first) and other labels (second) (e.g. 1:normal, 2:bold, 3:italic)
  # adj: numeric vector indicating alignment adjustment for text label: 0 is left, .5 is center, 1 is right.
  # cex: numeric vector for size of fonts for top panel (first) and bottom panel (Second)
  
  # 1. Set up
  if (T) {
    # Arrange data
    d <- data
    rownames(d) <- 1:nrow(d)
    
    # Create ordering vector
    if (order=="asis")       o <- 1:length(d[,index.est])
    if (order=="increasing") o <- order(d[,index.est])
    if (order=="decreasing") o <- order(-d[,index.est])
    if (!is.numeric(d[,index.est])) {warning("index.est does not point to a numeric vector.") ; break}
    d <- d[o,]
    est <- d[,index.est] # Estimate
    if (length(index.ci)>1) {
      l1 <- d[,index.ci[1]]
      h1 <- d[,index.ci[2]]
      if (length(index.ci)>2) {
        l2 <- d[,index.ci[3]]
        h2 <- d[,index.ci[4]]
      }
    } else {
      if (!is.numeric(d[,index.se]))  {warning("index.se does not point to a numeric vector.") ; break}
      se  <- d[,index.se] # Std error
      ci <- sort(ci)
      a <- qnorm(1-(1-ci)/2)
      l1 <- est - a[1]*se
      h1 <- est + a[1]*se
      if (length(ci)>1) {
        l2 <- est - a[2]*se
        h2 <- est + a[2]*se
      }
    }
    
    # Table
    if (length(index.ci)>1) remove.index <- c(index.est,index.ci) else remove.index <- c(index.est,index.se)
    remove.index <- remove.index[!is.na(remove.index)]
    tab <- t(d[,-remove.index]) # get only the relevant info for bottom panel
    if (!is.list(labels) & !is.character(labels)) labels <- rownames(tab)
    
    # Double check we have enough labels
    if ( nrow(tab) != length(unlist(labels))) {
      print("Warning: number of labels don't match number of models.")
      labels <- rownames(tab)
    }
    
    # Plotting objects
    xs <- 1:nrow(d) # the Xs for bars and dots
    if (n[1]>1 & length(n)==1) xs <- xs + ceiling(seq_along(xs)/n) - 1 # group models by n
    if (length(n)>1) {
      if (sum(n) != nrow(d) ) {
        warning("Group sizes don't add up.")
      } else {
        idx <- unlist(lapply(1:length(n), function(i) rep(i,n[i])))
        xs <- xs + idx - 1
      }
    }
    h <- nrow(tab) + ifelse(is.list(labels),length(labels),0) # number of rows in table
    # Location of data and labels
    if (is.list(labels)) {
      index <- unlist(lapply(1:length(labels), function(i) rep(i, length(labels[[i]])) ))
      locs <- split(1:length(index),index)
      locs <- lapply(unique(index), function(i) {
        x <- locs[[i]]+i-1
        x <- c(x,max(x)+1)
      })
      yloc  <- unlist(lapply(locs, function(i) i[-1])) # rows where data points are located
      yloc2 <- sapply(locs, function(i) i[1]) # rows where group lables are located
    } else {
      yloc <- 1:length(labels)
    }
    
    # Range
    if (is.na(ylim[1]) | length(ylim)!=2) {
      if (length(index.ci)>2 | length(ci)>1) {
        ylim <- range(c(l2,h2,ref)) # range that includes reference lines
      } else {
        ylim <- range(c(l1,h1,ref))
      }
      ylim <- ylim + diff(ylim)/10*c(-1,1) # and a bit more
    }
    xlim <- range(xs) #+ c(1,-1)
  }
  
  # 2. Plot
  if (T) {
    #par(mfrow=c(2,1), mar=c(0,leftmargin,0,0), oma=oma, xpd=F, family=family)
    layout(t(t(2:1)), height=heights, widths=1)
    par(mar=c(0,leftmargin,0,0), xpd=F)
    
    # Bottom panel (plotted first)
    plot(est, xlab="", ylab="", axes=F, type="n", ylim=c(h,1), xlim=xlim)
    lapply(1:nrow(tab), function(i) {
      # Get colors and point type
      type <- ifelse(is.na(tab[i,]),3,ifelse(tab[i,]==TRUE,1, ifelse(tab[i,]==FALSE,2,NA)))
      type <- ifelse(names(type) %in% paste(highlight) & type==1,4,type) # replace colors for baseline model
      col <- col.dot[type]
      bg  <- bg.dot[type]
      pch <- as.numeric(pch.dot[type])
      sel <- is.na(pch)
      # Plot points
      points(xs, rep(yloc[i],length(xs)), col=col, bg=bg, pch=pch, lwd=lwd.symbol)
      points(xs[sel], rep(yloc[i],length(xs))[sel], col=col[sel], bg=bg[sel], pch=pch.dot[3]) # symbol for missing value
      
    })
    par(xpd=T)
    if (is.list(labels)) text(-offset[1], yloc2, labels=names(labels), adj=adj[1], font=fonts[1], cex=cex[2])
    # Does not accomodate subscripts
    text(-rev(offset)[1], yloc , labels=unlist(labels), adj=rev(adj)[1], font=fonts[2], cex=cex[2])
    # Accomodates subscripts at the end of each string
    if (F) {
      labels1 <- unlist(labels)
      lapply(1:length(labels1), function(i) {
        a  <- labels1[i]
        a1 <- strsplit(a,"\\[|\\]")[[1]][1]
        a2 <- rev(strsplit(a,"\\[|\\]")[[1]])[1]
        if (identical(a1,a2))  a2 <- NULL
        text(-rev(offset)[1], yloc[i], labels=bquote(.(a1)[.(a2)]), adj=adj[2], font=fonts[2], cex=cex[2])
      })
    }
    par(xpd=F)
    
    # Top panel (plotted second)
    colvec  <- ifelse(colnames(tab) %in% paste(highlight), col.est[2], col.est[1])
    bg.colvec  <- ifelse(colnames(tab) %in% paste(highlight), bg.est[2], bg.est[1])
    colvec2 <- ifelse(colnames(tab) %in% paste(highlight),col.est2[2], col.est2[1])
    plot(est, xlab="", ylab="", axes=F, type="n", ylim=ylim, xlim=xlim)
    # Band if present
    if (!is.na(band.ref[1])) {
      rect(min(xlim)-diff(xlim)/10, band.ref[1], max(xlim)+diff(xlim)/10, band.ref[2], col=col.band.ref, border=NA)
    }
    # Reference lines
    abline(h=ref, lty=lty.ref, lwd=lwd.ref, col=col.ref)
    # Vertical bars
    if (length(ci)>1 | length(index.ci)>2) arrows(x0=xs, y0=l2, x1=xs, y1=h2, length=length, code=3, lwd=rev(lwd.est)[1], col=colvec2, angle=90)
    arrows(x0=xs, y0=l1, x1=xs, y1=h1, length=length, code=3, lwd=lwd.est[1]     , col=colvec, angle=90)
    points(xs, est, pch=pch.est, lwd=lwd.symbol, col=colvec, bg=bg.colvec)
    # Axes
    if (axes) {
      axis(2, las=2, cex.axis=cex[1], lwd=lwd.border)
      axis(4, labels=NA, lwd=lwd.border)
    }
    mtext(ylab, side=2, line=3.5, cex=cex[1])
    box(lwd=lwd.border)
    
  }
  
} 


## We now produce the "data" argument for the function schart, i.e. we run one regression, and bind in a dataframe 
# the coefficients, the SE, and indicator variables for the correponding specifications. 
make_spec_chart_df <- function(island,
                               start_year = 2002, 
                               end_year = 2014,
                               outcome_variable = "lucpfip_pixelcount", # LHS. One of "lucfip_pixelcount", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount"p
                               catchment = "CR",  
                               alt_cr = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                               nearest_mill = FALSE, # whether the ibs variables should be attributed to parcels as from the nearest mill or inverse distance weighted average.  
                               margin = "both",
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
                               fe = "lonlat + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                               offset = FALSE, # Logical. Should the log of the remaining forest be added as an offset.  
                               lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                               controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml"), # , "wa_prex_cpo_imp1"character vectors of names of control variables (don't specify lags in their names)
                               remaining_forest = FALSE, # Logical. If TRUE, the remaining forest is added as a control
                               interaction_terms = NULL, # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1"), # may be one or several of the controls specified above. 
                               interacted = "regressors",
                               interact_regressors = TRUE, # if there are two regressors (e.g. ffb and cpo), should their interaction be included in the model? 
                               pya_ov = FALSE, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                               illegal = "all", # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                               min_forest_2000 = 0, 
                               min_coverage = 0, # fraction, from 0 to 1. Minimum share of reachable IBS over all reachable (UML), for an obs.to be included in sample. 
                               weights = FALSE, # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
                               # additional argument to this function: 
                               variable = "cpo",
                               cluster = "reachable" # "subdistrict"
){
  ### Get coefficient and SE for a given specification passed to make_base_reg.
  # call make_base_reg
  res_data_list <- make_base_reg(island = island,
                                 start_year = start_year, 
                                 end_year = end_year,
                                 outcome_variable = outcome_variable,
                                 catchment = catchment,
                                 alt_cr = alt_cr,  
                                 nearest_mill = nearest_mill, 
                                 margin = margin,
                                 restr_marg_def = restr_marg_def,
                                 commo = commo, 
                                 x_pya = x_pya, 
                                 dynamics = dynamics, 
                                 yoyg = yoyg, 
                                 only_sr = only_sr,
                                 log_prices = log_prices,
                                 short_run = short_run, 
                                 imp = imp, 
                                 distribution = distribution, 
                                 fe = fe, 
                                 offset = offset,
                                 lag_or_not = lag_or_not, 
                                 controls = controls,
                                 remaining_forest = remaining_forest, 
                                 interaction_terms = interaction_terms, # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1"), # may be one or several of the controls specified above. 
                                 interacted = interacted,
                                 interact_regressors = interact_regressors, # if there are two regressors (e.g. ffb and cpo), should their interaction be included in the model? 
                                 pya_ov = pya_ov, 
                                 illegal = illegal, # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                                 min_forest_2000 = min_forest_2000, 
                                 min_coverage = min_coverage,
                                 weights = weights)
  
  # compute SE the way specified by cluster
  ape_mat <- make_APEs_1regr(res_data_list, CLUSTER = cluster, stddev = FALSE) 
  # reg_summary <- summary(reg_res, se = "cluster", cluster = cluster)
  # get coeff and SE
  # coeff_and_se <- reg_summary$coeftable[match(variable, commo), 1:2]
  
  ### make indicator variables that will be used to label specifications. 
  # /!\ THE ORDER HERE MATTERS !
  ind_var <- data.frame(
    # outcome variable
    #"larger_forest_def" = FALSE,
    #"rapid" = FALSE,
    ## Data cleaning
    "imp1" = FALSE,
    "lag_or_not" = FALSE,
    ## Sampling
    "min_forest_2000" = FALSE,
    "min_coverage" = FALSE,
    ## Catchment model
    #"nearest_mill" = FALSE,
    #"intensive_only" = FALSE,
    "alt_cr" = FALSE,
    "CA" = FALSE,
    # Price signal past year average
    "only_sr" = FALSE,
    "ya_3" = FALSE,
    "ya_4" = FALSE,
    "ya_5" = FALSE,
    # distribution assumptions
    "poisson" = FALSE,
    "quasipoisson" = FALSE,
    "negbin" = FALSE,
    # controls
    # "two_commo" = FALSE,
    # "price_dynamics" = FALSE,
    "control_own" = FALSE,
    "n_reachable_uml_control" = FALSE, 
    "ov_lag1" = FALSE,
    "ov_lag4" = FALSE,
    "ngb_ov_lag4" = FALSE,
    "prex_cpo_control" = FALSE,
    "baseline_forest_trend" = FALSE,
    #"remaining_forest" = FALSE,
    # interactions
    # "own_interact" = FALSE,
    # "n_reachable_uml_interact" = FALSE,
    # "prex_cpo_interact" = FALSE,
    # "remaining_forest_interact"  = FALSE,
    # weights
    #"weights" = FALSE,
    # fixed effects
    "unit_fe" = FALSE,
    "tw_fe" = FALSE,
    "unit_provyear_fe" = FALSE,
    "unit_distryear_fe" = FALSE,
    "unit_subdistryear_fe" = FALSE,
    # "unit_villageyear_fe" = FALSE,
    # standard errors
    "unit_cluster" = FALSE, 
    "reachable_cluster" = FALSE,
    "village_cluster" = FALSE,
    "subdistrict_cluster" = FALSE,
    "district_cluster" = FALSE,
    "twoway_cluster" = FALSE
  )
  
  # change the OV indicator variable to TRUE 
  #if(grepl("lucf", outcome_variable)){ind_var[,"larger_forest_def"] <- TRUE}
  #if(grepl("rapid",outcome_variable)){ind_var[,"rapid"] <- TRUE}
  
  ## Data cleaning
  # set the indicator variable for the data imputation 
  if(imp == 1){ind_var[,"imp1"] <- TRUE}
  # set indicator variable for lagging or not
  if(lag_or_not == "_lag1"){ind_var[,"lag_or_not"] <- TRUE}
  
  ## sampling
  # minimums 
  if(min_forest_2000>0){ind_var[,"min_forest_2000"] <- TRUE}
  if(min_coverage>0){ind_var[,"min_coverage"] <- TRUE}
  
  ## Catchment model
  # nearest mill rather than inverse-distance weighted average
  if(nearest_mill){ind_var[,"nearest_mill"] <- TRUE}
  # margin
  #if(margin == "intensive"){ind_var[,"intensive_only"] <- TRUE}
  # Catchment 
  if(alt_cr){ind_var[,"alt_cr"] <- TRUE} 
  if(catchment == "CA"){ind_var[,"CA"] <- TRUE} 
  
  ## Past year average
  # only Short Run
  if(only_sr){ind_var[,"only_sr"] <- TRUE}
  # set # year average indicator variables (+1 because we look at the total effect)
  ind_var[,grepl(paste0("ya_", x_pya+1), colnames(ind_var))] <- TRUE
  
  ## set the indicator variables for the distribution assumptions
  ind_var[,distribution] <- TRUE
  
  # set indicator variables for fixed effects
  if(fe == "lonlat"){ind_var[,"unit_fe"] <- TRUE}
  if(fe == "lonlat + year"){ind_var[,"tw_fe"] <- TRUE}
  if(fe == "lonlat + province_year"){ind_var[,"unit_provyear_fe"] <- TRUE}
  if(fe == "lonlat + district_year"){ind_var[,"unit_distryear_fe"] <- TRUE}
  if(fe == "lonlat + subdistrict_year"){ind_var[,"unit_subdistryear_fe"] <- TRUE}
  if(fe == "lonlat + village_year"){ind_var[,"unit_villageyear_fe"] <- TRUE}
  
  ## set indicator variables for controls
  # second commo 
  # if(length(commo) == 2){ind_var[,"two_commo"] <- TRUE}
  # if(dynamics){ind_var[,"price_dynamics"] <- TRUE}
  # ownership
  if(any(grepl("pct_own", controls))){ind_var[,"control_own"] <- TRUE}
  # n_reachable_uml
  if(any(grepl("n_reachable_uml", controls))){ind_var[,"n_reachable_uml_control"] <- TRUE}
  # 1-year lagged outcome variable
  if(any(grepl(paste0(outcome_variable, "_lag1"), controls))){ind_var[,"ov_lag1"] <- TRUE}
  # 4-year lagged outcome variable
  if(any(grepl(paste0(outcome_variable, "_lag4"), controls))){ind_var[,"ov_lag4"] <- TRUE}
  # spatial lag deforestation
  if(any(grepl("ngb_ov_lag4", controls))){ind_var[,"ngb_ov_lag4"] <- TRUE}
  # prex_cpo
  if(any(grepl("prex_cpo", controls))){ind_var[,"prex_cpo_control"] <- TRUE}
  # baseline forest trend
  if("baseline_forest_trend" %in% controls){ind_var[,"baseline_forest_trend"] <- TRUE}
  # remaining forest
  if("remain_pf_pixelcount" %in% controls){ind_var[,"remaining_forest"] <- TRUE}
  
  ## interactions
  # if(any(grepl("own_nat_priv",interaction_terms)) | any(grepl("own_for",interaction_terms))){ind_var[,"own_interact"] <- TRUE}
  # if(any(grepl("n_reachable_uml",interaction_terms))){ind_var[,"n_reachable_uml_interact"] <- TRUE}
  # if(any(grepl("prex_cpo",interaction_terms))){ind_var[,"prex_cpo_interact"] <- TRUE}
  # if(any(grepl("remain_pf",interaction_terms))){ind_var[,"remaining_forest_interact"] <- TRUE}
  
  
  # weights indicator variable
  #if(weights){ind_var[,"weights"] <- TRUE}
  #if(distribution == "negbin"){ind_var[,"weights"] <- FALSE} # turn it back to FALSE if negbin distribution
  
  # clustering
  if(length(cluster) == 1 & cluster == "lonlat"){ind_var[,"unit_cluster"] <- TRUE}
  if(length(cluster) == 1 & cluster == "reachable"){ind_var[,"reachable_cluster"] <- TRUE}
  if(length(cluster) == 1 & cluster == "village"){ind_var[,"village_cluster"] <- TRUE}
  if(length(cluster) == 1 & cluster == "district"){ind_var[,"district_cluster"] <- TRUE}
  if(length(cluster) == 1 & cluster == "subdistrict"){ind_var[,"subdistrict_cluster"] <- TRUE}
  if(length(cluster) == 2){ind_var[,"twoway_cluster"] <- TRUE}
  
  ### Bind together coeff, SE, and specification labels
  spec_df <- cbind(ape_mat[1,],ape_mat[2,], ind_var)
  
  return(spec_df)
}


### COMPUTE THE DATASETS FOR EACH ISLAND, AND COMMODITY OF INTEREST

# we make the specification charts for 6 main estimates
SIZE <- "a"
ISL <- "both"
# for(ISL in c("Sumatra", "Kalimantan", "both")){
#   for(SIZE in c("i","sm")){


reg_stats_indvar_list <- list()
i <- 1

### Add to the list the specifications that are to be compared in the chart. 
# These are particular departures from the preferred specification (which arguments are already set by default)
# Specifying only arguments for which it's changing + those defining the estimate we are interested in plotting the spec chart. 

## the main specification
reg_stats_indvar_list[["main"]] <- make_spec_chart_df(island = ISL,
                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"))
i <- i+1

## Sampling 
# Minimum forest coverage in 2000
reg_stats_indvar_list[["min_forest_2000"]] <- make_spec_chart_df(island = ISL,
                                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                 min_forest_2000 = 0.5)
i <- i+1

# Minimum IBS/UML ratio 
reg_stats_indvar_list[["min_coverage"]] <- make_spec_chart_df(island = ISL,
                                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                              min_coverage = 0.5)
i <- i+1

## Data cleaning
# Stronger imputations (imp2)
reg_stats_indvar_list[["imp"]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                     imp = 2)
i <- i+1


# Not lagged controls
reg_stats_indvar_list[["lag_or_not"]] <- make_spec_chart_df(island = ISL,
                                                            outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                            lag_or_not = "")
i <- i+1

## For an alternative forest definition
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucf",SIZE,"p_pixelcount")) 
# i <- i+1

## Alternative past year average
reg_stats_indvar_list[["only_sr"]] <- make_spec_chart_df(island = ISL,
                                                         outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                         only_sr = TRUE)
i <- i+1

for(XPYA in c(2, 4)){
  reg_stats_indvar_list[[paste0(XPYA,"_pya")]] <- make_spec_chart_df(island = ISL,
                                                                     outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                     x_pya = XPYA)
  i <- i+1
}

## Catchment model 

# Nearest mill
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  nearest_mill = TRUE)
# i <- i+1

# For intensive margin only
# reg_stats_indvar_list[["margin"]] <- make_spec_chart_df(island = ISL,
#                                                         outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                         margin = "intensive") 
# i <- i+1

# For alternative catchment radius
reg_stats_indvar_list[["alt_cr"]] <- make_spec_chart_df(island = ISL,
                                                        outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                        alt_cr = TRUE) 
i <- i+1  

# For catchment area
reg_stats_indvar_list[["catchment"]] <- make_spec_chart_df(island = ISL,
                                                           outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                           cluster = "subdistrict", # variable reachable is not computed in the CA stream. 
                                                           catchment = "CA") 
i <- i+1


i # 12

## RAPID LUCFP 
if(SIZE == "i"){
  # loops over critical parameters (not repeating because the outcome is different)
  for(IMP in c(1, 2)){
    for(XPYA in c(2, 3, 4)){
      for(LAG in c("_lag1", "")){
        reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                         outcome_variable = paste0("lucpfip_rapid_pixelcount"),
                                                         imp = IMP,
                                                         x_pya = XPYA,
                                                         lag_or_not = LAG)
        i <- i+1
      }
    }
  }    
  i # 28
  
  # in CA
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   outcome_variable = paste0("lucpfip_rapid_pixelcount"),
                                                   catchment = "CA") 
  i <- i+1
  
  # With 2way clustering
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   outcome_variable = paste0("lucpfip_rapid_pixelcount"),
                                                   SE = "twoway")
  i <- i+1
  
  # negbin
  # reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
  #                                                  outcome_variable = paste0("lucpfip_rapid_pixelcount"),
  #                                                  distribution = "negbin")
  # i <- i+1
  
  # unit FE
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   outcome_variable = paste0("lucpfip_rapid_pixelcount"),
                                                   fe = "lonlat")
  i <- i+1
}


## For alternative distributional assumptions
for(DISTR in c("poisson", "negbin")){#
  reg_stats_indvar_list[[paste0(DISTR,"_distribution")]] <- make_spec_chart_df(island = ISL,
                                                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                               distribution = DISTR)
  i <- i+1
}

## For alternative fixed effects
for(FE in c("lonlat", "lonlat + year", 
            "lonlat + province_year", "lonlat + subdistrict_year")){#, "lonlat + village_year"
  reg_stats_indvar_list[[paste0(FE," fe")]] <- make_spec_chart_df(island = ISL,
                                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                  fe = FE)
  i <- i+1
}
i  
## For an alternative standard error computation (two-way clustering)
c <- 1
for(CLT in list("lonlat", "village", "subdistrict", "district", c("lonlat","district_year"))){#"reachable" ,  
  reg_stats_indvar_list[[paste0("cluster_", c)]] <- make_spec_chart_df(island = ISL,
                                                                       outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                       cluster = CLT)
  c <- c + 1
}

i # 22

### For alternative control sets, (with their interaction each time)
# Without controls
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c(), 
#                                                  interaction_terms = NULL)

ov_lag1 <- paste0("lucpf",SIZE,"p_pixelcount_lag1")
ov_lag4 <- paste0("lucpf",SIZE,"p_pixelcount_lag4")

control_sets_list <- list(c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"), # Ownership control only
                          c("n_reachable_uml"), # remoteness control only
                          c(ov_lag1), # lagged outcome variable
                          c("wa_prex_cpo_imp1"), # prex cpo control only
                          c("baseline_forest_trend"), # baseline forest trend only
                          # c("remain_pf_pixelcount"), # remaining forest only
                          
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
                            ov_lag1), # Ownership and lagged ov control
                          c("n_reachable_uml", 
                            ov_lag1), # remoteness and lagged ov control                          
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
                            "wa_prex_cpo_imp1"), # Ownership and wa_prex_cpo_imp1 control
                          c("n_reachable_uml", 
                            "wa_prex_cpo_imp1"), # remoteness and wa_prex_cpo_imp1 control
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
                            "baseline_forest_trend"), # Ownership and baseline_forest_trend control
                          c("n_reachable_uml", 
                            "baseline_forest_trend"), # remoteness and baseline_forest_trend control
                          # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
                          #   "remain_pf_pixelcount"), # Ownership and remain_pf_pixelcount control
                          # c("n_reachable_uml", 
                          #   "remain_pf_pixelcount"), # remoteness and remain_pf_pixelcount control
                          
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                            ov_lag1), # Ownership, remoteness, and lagged ov control   
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                            "wa_prex_cpo_imp1"), # Ownership, remoteness, and wa_prex_cpo_imp1 control
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                            "baseline_forest_trend"), # Ownership, remoteness, and baseline_forest_trend control
                          # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                          #   "remain_pf_pixelcount"), # Ownership, remoteness, and remain_pf_pixelcount control
                          
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                            ov_lag1, "wa_prex_cpo_imp1"), # two basic + lagged ov + wa_prex_cpo_imp1 
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                            ov_lag1, "baseline_forest_trend"), # two basic + lagged ov + baseline_forest_trend
                          c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                            "wa_prex_cpo_imp1", "baseline_forest_trend") # two basic + wa_prex_cpo_imp1 + baseline_forest_trend 
)

c <- 1
for(ctrl in control_sets_list){
  reg_stats_indvar_list[[paste0("ctrl_",c)]] <- make_spec_chart_df(island = ISL,
                                                                   outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                   controls = ctrl)
  c <- c + 1
}

## Add the control set with the 4-year lagged deforestation 
reg_stats_indvar_list[[paste0("ctrl_",c)]] <- make_spec_chart_df(island = ISL,
                                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                 controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                                                                              ov_lag4))
c <- c + 1

## Add the control set with the 4-lagged deforestation in the neighbor grid cells
reg_stats_indvar_list[[paste0("ctrl_",c)]] <- make_spec_chart_df(island = ISL,
                                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                 controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
                                                                              "ngb_ov_lag4"))
c <- c + 1



# # two basic + wa_prex_cpo_imp1 and remain_pf_pixelcount
# ctrl <- c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
#           "wa_prex_cpo_imp1", "remain_pf_pixelcount")
# reg_stats_indvar_list[[paste0(ctrl, sep = "_")]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = ctrl)
# i <- i+1

# # two basic + baseline_forest_trend and remain_pf_pixelcount
# ctrl <- c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
#           "baseline_forest_trend", "remain_pf_pixelcount")
# reg_stats_indvar_list[[paste0(ctrl, sep = "_")]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = ctrl)
# i <- i+1

# ## All controls
# ctrl <- c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
#           "wa_prex_cpo_imp1", "baseline_forest_trend", "remain_pf_pixelcount")
# reg_stats_indvar_list[[paste0(ctrl, sep = "_")]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = ctrl)
# i <- i+1



## Now with less interactions
# # With all controls, but no interaction
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = NULL)
# i <- i+1
#   
# # With all controls, but only ownership interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"))
# i <- i+1
# 
# # With all controls, but only remoteness interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("n_reachable_uml"))
# i <- i+1
# 
# # With all controls, but only wa_prex_cpo_imp1 interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_prex_cpo_imp1"))
# i <- i+1
# 
# # With all controls, but only ownership and remoteness interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml"))
# i <- i+1
# 
# # With all controls, but only ownership and wa_prex_cpo_imp1 interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "wa_prex_cpo_imp1"))
# i <- i+1
# 
# # With all controls, but only remoteness and wa_prex_cpo_imp1 interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("n_reachable_uml", "wa_prex_cpo_imp1"))
# i <- i+1


# convert to dataframe to be able to chart
reg_stats_indvar <- bind_rows(reg_stats_indvar_list)
#reg_stats_indvar <- reg_stats_indvar[,-c(ncol(reg_stats_indvar))]

# save it 
if(sum(duplicated(reg_stats_indvar))==0 ){ # i.e. 50 currently & nrow(reg_stats_indvar)+1 == i
  saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",SIZE,"_02072021")))
} else{print(paste0("SOMETHING WENT WRONG in spec_chart_df_",ISL,"_",SIZE))}

#}else{reg_stats_indvar <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",ISL)))}
#   }
# }

# one_25 <- reg_stats_indvar[,c(1:25)] 
# nearest_mill <- one_25[,"nearest_mill"]
# one_25 <- one_25[,-4]
# one_25[,"nearest_mill"] <- nearest_mill
# new <- cbind(one_25, reg_stats_indvar[,c(26:34)])
# reg_stats_indvar <- new

# reg_stats_indvar <- dplyr::filter(reg_stats_indvar, nearest_mill == FALSE)
# reg_stats_indvar <- dplyr::select(reg_stats_indvar, -nearest_mill)

### PLOTTING 
### GIVE HERE THE ISLAND, THE OUTCOME AND THE DATE FOR WHICH YOU WANT THE SPEC CHART TO BE PLOTTED
scdf <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_both_a_02072021")))

# some modifications for now because scdf run with some "mistakes"
# scdf <- dplyr::select(scdf, -weights)
# scdf <- dplyr::mutate(scdf, own_interact = (own_for_interact | own_nat_priv_interact))
# scdf <- dplyr::select(scdf, -own_for_interact, -own_nat_priv_interact)
# scdf <- scdf[,c(1:23,29,24:28)]
# 
# scdf$ngb_ov_lag4 <- FALSE
# scdf <- cbind(scdf[,1:19], scdf[,ncol(scdf)], scdf[,20:31])
# names(scdf)[20] <- "ngb_ov_lag4"
# 
# colnames(scdf) %>% all.equal(colnames(reg_stats_indvar))
# scdf <- rbind(scdf, reg_stats_indvar)

# make labels for the chart
schart_labels <- list(#"Dependent variable:" = c("Larger forest definition"),
  "Data cleaning:" = c("Stronger imputations", 
                       "Lagged IBS variables"),
  
  "Sampling:" = c("Minimum forest cover in 2000",
                  "Minimum IBS to UML mill ratio"),
  
  "Catchment model" = c(#"Nearest mill",
    #"Intensive margin expansion only",
    "Alternative catchment radius",
    "2-hour driving catchment area"), 
  
  "Price signal averaged over:" = c("1 year",
                                    "3 years", 
                                    "4 years", 
                                    "5 years"),        
  
  "Distribution assumption:" = c("Poisson",
                                 "Quasi-Poisson", 
                                 "Negative Binomial"),
  
  
  "Controls:" = c(#paste0(toupper(c("ffb", "cpo")[!grepl(VAR, c("ffb", "cpo"))]), " price signal"), 
    "Ownership",
    "# reachable mills", 
    "1-y lag deforestation",
    "4-y lag deforestation",
    "4-y neighbors' deforestation",
    "% CPO exported", 
    "Baseline forest trend"),#"Remaining forest""Neighbors' outcomes, 4-year lagged"
  # "Interaction with:" = c(#paste0(toupper(c("ffb", "cpo")[!grepl(VAR, c("ffb", "cpo"))]), " price signal"), 
  #   "Ownership",
  #   "# reachable UML mills", 
  #   "% CPO exported", 
  #   "Remaining forest"),
  #"Weights" = "", 
  
  "Fixed effects:" = c("Plantation", 
                       "Plantation and year", 
                       "Plantation and province-year",
                       "Plantation and district-year", 
                       "Plantation and subdistrict-year"),
  
  "Level of SE clustering:" = c("Plantation", 
                                "Set of reachable mills",
                                "Village", 
                                "Subdistrict", 
                                "District",
                                "Plantation and district-year")
) 

# # remove rows where SEs could not be computed
# problems <- reg_stats_indvar[is.na(reg_stats_indvar$`Std. Error`),]
# reg_stats_indvar <- reg_stats_indvar[!is.na(reg_stats_indvar$`Std. Error`),]
# rownames(reg_stats_indvar) <- seq(1, nrow(reg_stats_indvar))


# find position of model to highlight in original data frame
a <- scdf
a <- a[#a$larger_forest_def==FALSE & 
  #a$rapid == FALSE &
  a$imp1 &
    a$lag_or_not &
    
    a$min_forest_2000 == FALSE & 
    a$min_coverage == FALSE & 
    
    #a$nearest_mill == FALSE & 
    #a$intensive_only == FALSE & 
    a$alt_cr == FALSE &
    a$CA == FALSE &
    
    a$only_sr == FALSE &
    a$ya_4 &
    
    #a$CR_30km &  
    a$quasipoisson &
    
    a$unit_distryear_fe &
    #a$two_commo == FALSE & 
    a$control_own & 
    a$n_reachable_uml_control  &
    a$ov_lag1 == FALSE &
    a$ov_lag4 == FALSE &
    a$ngb_ov_lag4 == FALSE & 
    a$prex_cpo_control == FALSE &
    a$baseline_forest_trend == FALSE &
    # a$remaining_forest == FALSE & 
    # a$ngb_ov_lag4 == FALSE & 
    # a$own_interact & 
    # a$n_reachable_uml_interact & 
    # a$prex_cpo_interact & 
    # a$remaining_forest_interact == FALSE & 
    #a$weights == FALSE & 
    a$reachable_cluster, ] 

model_idx <- row.names(a)
# this is our baseline model(s)
a[model_idx,]

# are there duplicated estimates (this is not expected) 
scdf[duplicated(scdf[,c(1,2)]),] # Poisson may be exactly the same as quasipoisson at the level of rounding CIs

schart(scdf, 
       labels = schart_labels,
       order="increasing",
       heights=c(.6,1.5),
       pch.dot=c(20,20,20,20),
       ci=c(.95),
       cex=c(0.6,0.7),
       highlight=model_idx, 
       #col.est=c(rgb(0,0,0,0.1),rgb(0,0.2,0.6, 0.1)),
       col.est=c("grey70", "red3"),
       #col.est2=c(rgb(0,0,0,0.08),"lightblue"),
       #col.dot=c(rgb(0,0,0,0.12),"grey95","grey95",rgb(0,0.4,0.6,0.3)),
       col.dot=c("grey70","grey95","grey95","red3"),
       #bg.dot=c(rgb(0,0,0,0.12),"grey95","grey95",rgb(0,0.4,0.6,0.3)),
       bg.dot=c("grey60","grey95","grey95","white"),
       adj=c(1,1), 
       #offset=c(10.4,10.1),
       leftmargin = 12,
       ylab = "Elasticity estimates",
       lwd.est = 5.8,
       lwd.symbol = 1
       #pch.est=20
)

# Should be saved in .pdf landscape A4 for correct vizualization under Latex.  





#### COUNTERFACTUAL APPLICATION #### 

# Produce the approximate population sample

# same sample selection but with some enlargements: 
# those within CR (intensive margin) of UML mills, 
# with positive deforestation at least one year, 
# with NAs in any covariate, 
# with remaining positive forest extent (but that should not change the cross-section size) 
# outside RSPO. 
# Rather annual than full period. This means that from the above selection, we take the number of grid cells. 
# Or the annual average number of grid cells. The former is better in line with the idea of effect on Indonesian oil palm sector as it was at the end of the study period. 


# Grid cells within 50km of a known (UML) mill. 
lucpfip <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                    parcel_size/1000,"km_",
                                    "82km_UML_CR.rds")))
lucpfsmp <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfsmp_panel_",
                                     parcel_size/1000,"km_",
                                     "82km_UML_CR.rds")))
lucpfip_dyn <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_dynamics_",
                                        parcel_size/1000,"km_",
                                        "82km_UML_CR.rds")))

# keep only year before 2015 (after they mean nothing since we plantation data are from 2015)
# and after 2001 (because the d_clean final sample used for estimation runs from 2002 due 4-year average variables being only observed from 2002 on.)
lucpfip <- lucpfip[lucpfip$year < 2015 & lucpfip$year > 2001,] 
lucpfip_dyn <- lucpfip_dyn[lucpfip_dyn$year < 2015 & lucpfip_dyn$year > 2001,] 
lucpfsmp <- lucpfsmp[lucpfsmp$year < 2015 & lucpfsmp$year > 2001,] 

# remove coordinates, except in one data frame
lucpfsmp <- dplyr::select(lucpfsmp, -lat, -lon, -idncrs_lat, -idncrs_lon)
lucpfip_dyn <- dplyr::select(lucpfip_dyn, -lat, -lon, -idncrs_lat, -idncrs_lon)

# join them
d <- inner_join(lucpfip, lucpfip_dyn, by = c("lonlat", "year")) 
d <- inner_join(d, lucpfsmp, by = c("lonlat", "year")) 
rm(lucpfip, lucpfsmp, lucpfip_dyn)

## make a variable that counts rapid and slow lucpfp events (this is only computed for primary forest as of now)
# THIS IS GOING TO BE THE MAIN OUTCOME VARIABLE INSTEAD OF lucpfip_pixelcount_total
d$lucpfip_pixelcount <- d$lucpfip_rapid_pixelcount + d$lucpfip_slow_pixelcount

# make variable that counts lucpfp events on both small and medium sized plantations 
d$lucpfsmp_pixelcount <- d$lucpfsp_pixelcount_total + d$lucpfmp_pixelcount_total


## make a variable that counts lucpfp events on all types of plantations
d <- mutate(d, lucpfap_pixelcount = lucpfip_pixelcount + lucpfsmp_pixelcount) 


## MAKE THE lonlat 
uni_lonlat <- unique(d$lonlat)
d <- mutate(d, 
            lonlat = match(lonlat, uni_lonlat)) 
rm(uni_lonlat)

# some arrangements
d <- dplyr::arrange(d, lonlat, lonlat, year)
row.names(d) <- seq(1,nrow(d))


### RESTRICT THE GRID CELL SAMPLE 

## to islands of interest 
# d has no island variable so far. But it is the merger of prepare_ scripts run over Sumatra and Kalimantan, which are the islands of interest. 



## to within 50km a UML mill
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
uml <- uml[!is.na(uml$lat),]
uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
uml <- st_transform(uml, crs = indonesian_crs)
uml <- st_buffer(uml, dist = 5e4)
uml <- st_union(uml)
# uml <- st_as_sf(uml)
# uml$CR <- TRUE

d <- st_as_sf(d, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs)

d_cs <- d[!duplicated(d$lonlat), ]
sgbp <- st_within(d_cs, uml)
d_cs$CR <- rep(FALSE, nrow(d_cs))
d_cs$CR[lengths(sgbp) > 0] <- TRUE
sum(d_cs$CR)

d_cs <- st_drop_geometry(d_cs)
d <- left_join(d, d_cs[, c("lonlat", "CR")], by = "lonlat")
d_save <- d
#d <- filter(d, d$CR)


## With some deforestation at least once 
# If it's smallholders
d$smallholders <- (d$lonlat %in% d[-obs2remove(fml = as.formula(paste0("lucpfsmp_pixelcount ~ lonlat")),
                                               data = d, 
                                               family = "poisson"),]$lonlat )


d <- d[-obs2remove(fml = as.formula(paste0("lucpfap_pixelcount ~ lonlat")),
                   data = d, 
                   family = "poisson"),]




## With some variables possibly NAs - ... need not do anything 



## That are not in RSPO
rspo <- st_read("input_data/RSPO_supply_bases/RSPO-certified_oil_palm_supply_bases_in_Indonesia.shp")
rspo <- st_transform(rspo, crs = indonesian_crs)
# set up a year variable $
# we assume that the year of probable behavior change is the year after the date of the icletdate variable
# BUT THIS MAY CHANGE 
rspo$year <- sub(pattern = ".*, ", 
                 replacement = "", 
                 x = rspo$icletdate)
rspo$year[rspo$year=="no letter"] <- "2015"
rspo$year <- rspo$year %>% as.integer()
rspo$year <- rspo$year +1

d$rspo_cert <- rep(FALSE, nrow(d))

n_with_rspo <- nrow(d)
gc_with_rspo <- length(unique(d$lonlat))

for(y in min(rspo$year):max(d$year)){
  # select d from the given year
  d_cs <- d[d$year == y, c("lonlat", "year", "rspo_cert")]
  # select supply bases already certified this year
  rspo_cs <- rspo[rspo$year <= y,] %>% st_geometry()
  
  sgbp <- st_within(x = d_cs, y = rspo_cs)
  
  d$rspo_cert[d$year == y][lengths(sgbp) == 1] <- TRUE
}

# - are not in an intended RSPO-certified supply base
d <- d[d$rspo_cert==FALSE,]

# how many in RSPO ? It removes observations but no grid cell, because all RSPO grid cells have obs. 
# before RSPO started
n_with_rspo - nrow(d)
gc_with_rspo - length(unique(d$lonlat))



## And finally identify illegal deforestation in this population, but do not restrict the sample to it yet. 
# OIL PALM CONCESSIONS
cns <- st_read(file.path("input_data/oil_palm_concessions"))
cns <- st_transform(cns, crs = indonesian_crs)

# LEGAL LAND USE 
llu <- st_read(file.path("input_data/kawasan_hutan/Greenorb_Blog/final/KH-INDON-Final.shp"))
llu <- st_transform(llu, crs = indonesian_crs)
unique(llu$Fungsi)
names(llu)[names(llu) == "Fungsi"] <- "llu"

# restrict llu to provinces of interest
llu <- llu[llu$Province == "Sumatra Utara" |
             llu$Province == "Riau" |
             llu$Province == "Sumatra Selantan" |
             llu$Province == "Papua Barat" |
             llu$Province == "Kalimantan Timur" |
             llu$Province == "Kalimantan Selatan" |
             llu$Province == "Kalimantan Tengah" |
             llu$Province == "Kalimantan Barat" |
             llu$Province == "Bengkulu" |
             llu$Province == "Lampung" |
             llu$Province == "Jambi" |
             llu$Province == "Bangka Belitung" |
             llu$Province == "Kepuluan Riau" |
             llu$Province == "Sumatra Barat" |
             llu$Province == "Aceh", ]


# OIL PALM CONCESSIONS
# We do not observe whether a grid cell is within a concession annually. 
# Therefore we only proceed with a cross section
d_cs <- d[!duplicated(d$lonlat),]
sgbp <- st_within(d_cs, cns)
d_cs$concession <- rep(FALSE, nrow(d_cs))
d_cs$concession[lengths(sgbp) > 0] <- TRUE

d_cs <- st_drop_geometry(d_cs)

d <- left_join(d, d_cs[,c("lonlat", "concession")], by = "lonlat")

# note that some d fall within more than one concession record. There may be several reasons for concession overlaps 
# like renewal of concession, with our withour aggrandisement. For our purpose, it only matters that there is at least one 
# concession record. 


# LEGAL LAND USE 
# this is quite long (~5min)
d_cs <- d[!duplicated(d$lonlat),] # the point of this is to give back spatial class to d_cs
d_cs <- st_join(x = d_cs, 
                y = st_make_valid(llu[,"llu"]), # st_make_valid bc thrown error otherwise 
                join = st_within, 
                left = TRUE)

# some grid cells seem to fall within overlapping llu shapes though. 
# It's really marginal (12 instances). Just remove the duplicates it produces.
d_cs <- d_cs[!duplicated(d_cs$lonlat),]

# merge back with panel 
d <- st_drop_geometry(d)
d_cs <- st_drop_geometry(d_cs)

d <- left_join(d, d_cs[,c("lonlat", "llu")], by = "lonlat")

unique(d$llu)

# ILLEGAL LUCFP 
# one possible link to shed light on accronyms http://documents1.worldbank.org/curated/pt/561471468197386518/pdf/103486-WP-PUBLIC-DOC-107.pdf

d <- dplyr::mutate(d,
                   illegal1 = (!concession & (llu != "HPK" | llu == "<NA>")), # it's not in concession and not in a convertible forest zone
                   illegal2 = (!concession & (llu == "KSA/KPA" | # it's not in concession and it's in a permanent forest zone designation
                                                llu == "KSA" |
                                                llu == "KPA" |
                                                llu == "KSAL" |
                                                llu == "HP" |
                                                llu == "HPT" |
                                                llu == "HL")))

# yields many missing in illegal because many grid cells are within a mising land use legal classification
# parcels[!duplicated(parcels$lonlat) & !is.na(parcels$llu), c("lonlat", "concession", "llu", "illegal1", "illegal2")]

# restrict the sample to the pre-defined legal status 






### COUNT THE AVERAGE ANNUAL NUMBER OF GRID CELLS ###
aggr_factor_all <- ddply(d, "year", summarise,
                         annual_avg_n = length(unique(lonlat)))[,"annual_avg_n"] %>% mean() %>% round(0)

aggr_factor_sm <- ddply(d[d$smallholders, ], "year", summarise,
                        annual_avg_n = length(unique(lonlat)))[,"annual_avg_n"] %>% mean() %>% round(0)

aggr_factor_ill <- ddply(d[!is.na(d$illegal2) & d$illegal2 == TRUE, ], "year", summarise,
                         annual_avg_n = length(unique(lonlat)))[,"annual_avg_n"] %>% mean() %>% round(0)
# because of the many NAs in illegal variables, we do not know the legal status for a large proportion of griod cells. 
# therefore we cannot say that "of the total aggregated effect, a given amount comes from illegal deforestation". 

### COMPUTE THE AVERAGE PARTIAL EFFECT ###




res_data_all <- make_base_reg(island = "both",
                              outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                              output_full = FALSE)

res_data_sm <- make_base_reg(island = "both",
                             outcome_variable = paste0("lucpfsmp_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                             output_full = FALSE)

res_data_ill <- make_base_reg(island = "both",
                              outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                              ill = "ill2",
                              output_full = FALSE)


make_counterfactuals <- function(res_data, aggr_factor){
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it is necessary that d_clean is in the main environment so that vcov.fixest can retrieve cluster variable
  coeff <- reg_res$coefficients[1]
  names(coeff) <- NULL
  
  avg_fv_ha <- mean(reg_res$fitted.values)*pixel_area_ha #  11.81169 ha
  
  avg_ha <- mean(d$lucpfap_pixelcount)*pixel_area_ha # 11.75246 ha
  
  print(avg_fv_ha*aggr_factor)
  # avg_ha*aggr_factor
  
  # one std. dev. in price 
  # remove fixed effect variations from the regressor
  reg_sd <- fixest::feols(fml = as.formula(paste0(
    names(reg_res$coefficients)[1],
    " ~ 1 | ", 
    paste0(reg_res$fixef_vars, collapse = "+"))),
    data = d_clean)
  
  # and take the remaining standard deviation
  sd_price_change <- sd(reg_sd$residuals)
  print(paste0("std. dev. price change is 100* : ",sd_price_change))
  
  # Price decrease necessary to halve deforestation (in pct.)
  delta_defo <- -0.29 # Indonesian unconditional target in NDC = 29% below BAU by 2030 incl. LULUCF
  tax <- (((1 + delta_defo)^(1/coeff)) - 1)
  
  # store application results in table
  price_changes <- c(sd_price_change, -0.01, tax)
  app_mat <- matrix(nrow = 3, ncol = length(price_changes))
  row.names(app_mat) <- c("Relative change (%)", "Total change (ha)", "Potential CO2 revenues (M$)")
  colnames(app_mat) <- as.character(price_changes)
  
  for(rel_price_change in price_changes){
    effect <- ((1+rel_price_change)^(coeff) - 1)
    app_mat[1, match(rel_price_change, price_changes)] <- effect*100
    
    defo <- (effect)*avg_fv_ha*aggr_factor
    app_mat[2, match(rel_price_change, price_changes)] <- defo
    
    app_mat[3, match(rel_price_change, price_changes)] <- -defo*174*(44/12)*5*1e-6
    # 174 is the Carbon loss associated with deforestation, in ton (Guillaume et al. 2018)
    # 44/12 is the C to CO2 conversion factor 
    # 5 is the CO2 price paid by Norway in the first REDD+ payment in $/tCO2eq. https://www.regjeringen.no/en/aktuelt/noreg-betaler-530-millionar-for-redusert-avskoging-i-indonesia/id2722135/
    # We give the final value in million $
  }
  
  app_mat[1,] <- app_mat[1,] %>% round(2)
  app_mat[2,] <- app_mat[2,] %>% formatC(digits = 0, format = "f")
  app_mat[3,] <- as.numeric(app_mat[3,]) %>% formatC(digits = 1, format = "f")
  
  return(app_mat) 
}
cf_mat_all <- make_counterfactuals(res_data_all, aggr_factor = aggr_factor_all)
cf_mat_sm <- make_counterfactuals(res_data_sm, aggr_factor = aggr_factor_sm)
cf_mat_ill <- make_counterfactuals(res_data_ill, aggr_factor = aggr_factor_ill)

cf_mat <- cbind(cf_mat_sm, cf_mat_ill, cf_mat_all)
# 
# # opportunity cost for the oil palm sector associated with a 1% decrease in prices
# ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# annual_avg_cpo_price <- mean(ibs$cpo_price_imp1, na.rm= T)
# annual_avg_cpo_production <- 19963.57143 * 1000 # between 2001 and 2014, according to USDA (which gives it in 1000 MT): https://www.indexmundi.com/agriculture/?country=id&commodity=palm-oil&graph=production
# cpo_oc <- 0.01 * annual_avg_cpo_price * annual_avg_cpo_production


# ACTUALLY IT IS NOT RELEVANT TO DISPLAY ILLEGAL AND SM, 
# because aggr_factor_ill is underestimated because many missings on illegal, so the participation of illegal to the total change
# is not meaningfull. 
cf_mat <- cf_mat_all
colnames(cf_mat) <- NULL
options(knitr.table.format = "latex")
kable(cf_mat, booktabs = T, align = "r",
      caption = "Counterfactual annual effects of different CPO price changes on deforestation in Indonesia") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "+1 std. dev." = 1, # 
                     "-1%" = 1,
                     "-19%" = 1
                     # "+1 std. dev." = 1, # 
                     # "-19%" = 1,
                     # "+1 std. dev." = 1, # 
                     # "-19%" = 1
  ),
  bold = F,
  align = "c") %>%
  
  column_spec(column = 1,
              width = "9em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(cf_mat))),
              width = "5em",
              latex_valign = "b")

rm(cf_mat)
