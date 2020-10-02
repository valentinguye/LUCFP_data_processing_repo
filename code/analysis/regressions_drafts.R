
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table", 
                   "foreign", "readstata13", "readxl",
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


# # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
troublePackages <- c("leaflet", "leaflet.providers", "png")
# Attempt to load packages from user's default libraries.
lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing

### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/reg_results")


### INDONESIAN CRS 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE
parcel_size <- 3000

### SET NUMBER OF THREADS USED BY {FIXEST} TO ONE (TO REDUCE R SESSION CRASH)
getFixest_nthreads()

### Set dictionary for names of variables to display in regression tables 
setFixest_dict(c(parcel_id = "grid cell",
                 lucpfip_ha_total = "LUCPFIP (ha)", 
                 lucpfip_pixelcount_total = "Land use change from primary forest to industrial oil palm plantations (LUCPFIP)", 
                 lucfip_ha_30th = "LUCFIP (30 pct. canopy density, ha)",
                 lucfip_ha_60th = "LUCFIP (60 pct. canopy density, ha)",
                 lucfip_ha_90th = "LUCFIP (90 pct. canopy density, ha)",
                 lucfip_pixelcount_30th = "LUCFIP (30 pct. canopy density, pixels)",
                 lucfip_pixelcount_60th = "LUCFIP (60 pct. canopy density, pixels)",
                 lucfip_pixelcount_90th = "LUCFIP (90 pct. canopy density, pixels)",
                 # No time dynamics FFB variables 
                 wa_ffb_price_imp1_3ya = "FFB price signal, 3 year average",
                 wa_ffb_price_imp1_3ya_lag1 = "FFB price signal, 3 year average (lagged)",
                 wa_ffb_price_imp1_4ya = "FFB price signal, 4 year average",
                 wa_ffb_price_imp1_4ya_lag1 = "FFB price signal, 4 year average (lagged)",
                 wa_ffb_price_imp1_5ya = "FFB price signal, 5 year average",
                 wa_ffb_price_imp1_5ya_lag1 = "FFB price signal, 5 year average (lagged)",
                 wa_ffb_price_imp1_yoyg_3ya = "FFB price signal y-o-y growth rate, 3 year average",
                 wa_ffb_price_imp1_yoyg_3ya_lag1 = "FFB price signal y-o-y growth rate, 3 year average (lagged)",
                 wa_ffb_price_imp1_yoyg_4ya = "FFB price signal y-o-y growth rate, 4 year average",
                 wa_ffb_price_imp1_yoyg_4ya_lag1 = "FFB price signal y-o-y growth rate, 4 year average (lagged)",
                 wa_ffb_price_imp1_yoyg_5ya = "FFB price signal y-o-y growth rate, 5 year average",
                 wa_ffb_price_imp1_yoyg_5ya_lag1 = "FFB price signal y-o-y growth rate, 5 year average (lagged)",
                 # No time dynamics CPO variables 
                 wa_cpo_price_imp1_3ya = "CPO price signal, 3 year average",
                 wa_cpo_price_imp1_3ya_lag1 = "CPO price signal, 3 year average (lagged)",
                 wa_cpo_price_imp1_4ya = "CPO price signal, 4 year average",
                 wa_cpo_price_imp1_4ya_lag1 = "CPO price signal, 4 year average (lagged)",
                 wa_cpo_price_imp1_5ya = "CPO price signal, 5 year average",
                 wa_cpo_price_imp1_5ya_lag1 = "CPO price signal, 5 year average (lagged)",
                 wa_cpo_price_imp1_yoyg_3ya = "CPO price signal y-o-y growth rate, 3 year average",
                 wa_cpo_price_imp1_yoyg_3ya_lag1 = "CPO price signal y-o-y growth rate, 3 year average (lagged)",
                 wa_cpo_price_imp1_yoyg_4ya = "CPO price signal y-o-y growth rate, 4 year average",
                 wa_cpo_price_imp1_yoyg_4ya_lag1 = "CPO price signal y-o-y growth rate, 4 year average (lagged)",
                 wa_cpo_price_imp1_yoyg_5ya = "CPO price signal y-o-y growth rate, 5 year average",
                 wa_cpo_price_imp1_yoyg_5ya_lag1 = "CPO price signal y-o-y growth rate, 5 year average (lagged)",
                 # SR FFB variables
                 wa_ffb_price_imp1 = "FFB price signal",
                 wa_ffb_price_imp1_lag1 = "FFB price signal (lagged)",
                 wa_ffb_price_imp1_yoyg = "FFB price signal y-o-y growth rate",
                 wa_ffb_price_imp1_yoyg_lag1 = "FFB price signal y-o-y growth rate (lagged)",
                 wa_ffb_price_imp1_dev_2pya = "FFB price signal deviation from 2 past year average",
                 wa_ffb_price_imp1_dev_2pya_lag1 = "FFB price signal deviation from 2 past year average (lagged)",
                 wa_ffb_price_imp1_dev_3pya = "FFB price signal deviation from 3 past year average",
                 wa_ffb_price_imp1_dev_3pya_lag1 = "FFB price signal deviation from 3 past year average (lagged)",
                 wa_ffb_price_imp1_dev_4pya = "FFB price signal deviation from 4 past year average",
                 wa_ffb_price_imp1_dev_4pya_lag1 = "FFB price signal deviation from 4 past year average (lagged)",
                 # LR FFB variables
                 wa_ffb_price_imp1_2pya = "FFB price signal, 2 past year average",
                 wa_ffb_price_imp1_2pya_lag1 = "FFB price signal, 2 past year average (lagged)",
                 wa_ffb_price_imp1_yoyg_2pya = "FFB price signal y-o-y growth rate, 2 past year average",
                 wa_ffb_price_imp1_yoyg_2pya_lag1 = "FFB price signal y-o-y growth rate, 2 past year average (lagged)",
                 wa_ffb_price_imp1_3pya = "FFB price signal, 3 past year average",
                 wa_ffb_price_imp1_3pya_lag1 = "FFB price signal, 3 past year average (lagged)",
                 wa_ffb_price_imp1_yoyg_3pya = "FFB price signal y-o-y growth rate, 3 past year average",
                 wa_ffb_price_imp1_yoyg_3pya_lag1 = "FFB price signal y-o-y growth rate, 3 past year average (lagged)",
                 wa_ffb_price_imp1_4pya = "FFB price signal, 4 past year average",
                 wa_ffb_price_imp1_4pya_lag1 = "FFB price signal, 4 past year average (lagged)",
                 wa_ffb_price_imp1_yoyg_4pya = "FFB price signal y-o-y growth rate, 4 past year average",
                 wa_ffb_price_imp1_yoyg_4pya_lag1 = "FFB price signal y-o-y growth rate, 4 past year average (lagged)",
                 # SR CPO variables
                 wa_cpo_price_imp1 = "CPO price signal",
                 wa_cpo_price_imp1_lag1 = "CPO price signal (lagged)",
                 wa_cpo_price_imp1_yoyg = "CPO price signal y-o-y growth rate",
                 wa_cpo_price_imp1_yoyg_lag1 = "CPO price signal y-o-y growth rate (lagged)",
                 wa_cpo_price_imp1_dev_2pya = "CPO price signal deviation from 2 past year average",
                 wa_cpo_price_imp1_dev_2pya_lag1 = "CPO price signal deviation from 2 past year average (lagged)", 
                 wa_cpo_price_imp1_dev_3pya = "CPO price signal deviation from 3 past year average",
                 wa_cpo_price_imp1_dev_3pya_lag1 = "CPO price signal deviation from 3 past year average (lagged)", 
                 wa_cpo_price_imp1_dev_4pya = "CPO price signal deviation from 4 past year average",
                 wa_cpo_price_imp1_dev_4pya_lag1 = "CPO price signal deviation from 4 past year average (lagged)", 
                 # LR CPO variables
                 wa_cpo_price_imp1_2pya = "CPO price signal, 2 past year average",
                 wa_cpo_price_imp1_2pya_lag1 = "CPO price signal, 2 past year average (lagged)",
                 wa_cpo_price_imp1_yoyg_2pya = "CPO price signal y-o-y growth rate, 2 past year average",
                 wa_cpo_price_imp1_yoyg_2pya_lag1 = "CPO price signal y-o-y growth rate, 2 past year average (lagged)",
                 wa_cpo_price_imp1_3pya = "CPO price signal, 3 past year average",
                 wa_cpo_price_imp1_3pya_lag1 = "CPO price signal, 3 past year average (lagged)",
                 wa_cpo_price_imp1_yoyg_3pya = "CPO price signal y-o-y growth rate, 3 past year average",
                 wa_cpo_price_imp1_yoyg_3pya_lag1 = "CPO price signal y-o-y growth rate, 3 past year average (lagged)",
                 wa_cpo_price_imp1_4pya = "CPO price signal, 4 past year average",
                 wa_cpo_price_imp1_4pya_lag1 = "CPO price signal, 4 past year average (lagged)",
                 wa_cpo_price_imp1_yoyg_4pya = "CPO price signal y-o-y growth rate, 4 past year average",
                 wa_cpo_price_imp1_yoyg_4pya_lag1 = "CPO price signal y-o-y growth rate, 4 past year average (lagged)",
                 ## controls
                 lucpfip_pixelcount_total_lag1 = "LUCPFIP (pixels, lagged)",
                 lucpfip_pixelcount_total_2pya = "LUCPFIP (pixels, 2 past year average)",
                 lucpfip_pixelcount_total_3pya = "LUCPFIP (pixels, 3 past year average)",
                 lucpfip_pixelcount_total_4pya = "LUCPFIP (pixels, 4 past year average)",
                 n_reachable_uml = "# reachable UML mills",
                 n_reachable_uml_lag1 = "# reachable UML mills (lagged)",
                 wa_pct_own_cent_gov_imp = "Local government mill ownership (pct.)",
                 wa_pct_own_cent_gov_imp_lag1 = "Local government mill ownership (pct., lagged)",
                 wa_pct_own_loc_gov_imp = "Local government mill ownership (pct.)",
                 wa_pct_own_loc_gov_imp_lag1 = "Local government mill ownership (pct., lagged)",
                 wa_pct_own_nat_priv_imp = "Domestic private mill ownership (pct.)",
                 wa_pct_own_nat_priv_imp_lag1 = "Domestic private mill ownership (pct., lagged)",
                 wa_pct_own_for_imp = "Foreign mill ownership (pct.)",
                 wa_pct_own_for_imp_lag1 = "Foreign mill ownership (pct., lagged)", 
                 wa_prex_cpo_imp1 = "Percentage CPO exported",
                 wa_prex_cpo_imp1_lag1 = "Percentage CPO exported (lagged)",
                 wa_prex_cpo_imp2 = "Percentage CPO exported",
                 wa_prex_cpo_imp2_lag1 = "Percentage CPO exported (lagged)"
))



# OK SO NOW FIXED EFFECT VARIABLES ARE DEFINED IN add_parcel_variables.R AS geographic_year 
# the only change implied is that now district^year can also be specified as district_year.


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


catchment_radius <- 5e4
island <- c("Sumatra", "Kalimantan", "Papua")
outcome_variable <- "lucpfip_pixelcount_total"
commo <- c("ffb","cpo")
dynamics <- "both"

##### COMPARE DISTRIBUTIONS / ESTIMATORS ##### 

# # gauss and ols are equivalent
# gauss <- fixest::feglm(fml, data = d,
#               family = "gaussian")
# 
# ols <- fixest::feols(fml, data=d)
# 
# # gauss with log outcome variable and poisson/quasipoisson wth log ling are *comparable* but not 
# # necessarily equivalent. 
# est_pois <- fixest::feglm(fml, data = d,
#                           family = "poisson")
# est_gaus = update(est_pois, log(lucpfip_pixelcount_total+1) ~ ., family="gaussian")
# est_qpois <- fixest::feglm(fml, data = d,
#                            family = "quasipoisson")
# etable(est_pois, est_qpois, est_gaus, se = "White")
# 
# ### ### ### ###

compare_estimators <- function(catchment_radius, island, outcome_variable, commo, dynamics){
  
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))
  # names(d)
  # 
  # length(unique(d$parcel_id)) # length cross sections: 47648
  # unique(d$year) # years from 2001 to 2015
  
  ## Island subsampling
  d <- d[d$island %in% island,]
  
  ## Base LU subsampling
  # df30 <- d[d$any_fc2000_30th,]
  # length(unique(df30$parcel_id)) # 41844
  # dpf <- d[d$any_pfc2000_total,]
  # length(unique(dpf$parcel_id)) # 24959 
  if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
    d <- d[d$any_pfc2000_total,]
  }
  if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
     outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
     outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
    d <- d[d$any_fc2000_30th,]
  }
  

  ## Specifications
  if(outcome_variable == "lucpfip_pixelcount_total"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucpfip_pixelcount_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucpfip_pixelcount_total ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucpfip_pixelcount_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    if(commo == "both" & dynamics == "SR"){
      specification <- lucpfip_pixelcount_total ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
        
    if(commo == "both" & dynamics == "both"){
      specification <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  if(outcome_variable == "lucpfip_ha_total"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucpfip_ha_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucpfip_ha_total ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucpfip_ha_total ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "SR"){
      specification <- lucpfip_ha_total ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "both"){
      specification <- lucpfip_ha_total ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  if(outcome_variable == "lucfip_pixelcount_30th"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucfip_pixelcount_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucfip_pixelcount_30th ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucfip_pixelcount_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "SR"){
      specification <- lucfip_pixelcount_30th ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "both"){
      specification <- lucfip_pixelcount_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  if(outcome_variable == "lucfip_ha_30th"){
    if(commo == "ffb" & dynamics == "LR"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "SR"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "ffb" & dynamics == "both"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "LR"){
      specification <- lucfip_ha_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "SR"){
      specification <- lucfip_ha_30th ~ wa_cpo_price_imp1_dev_3pya_lag1 + 
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "cpo" & dynamics == "both"){
      specification <- lucfip_ha_30th ~ wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "SR"){
      specification <- lucfip_ha_30th ~  wa_ffb_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "LR"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1  + wa_cpo_price_imp1_yoyg_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
    
    if(commo == "both" & dynamics == "both"){
      specification <- lucfip_ha_30th ~ wa_ffb_price_imp1_yoyg_3pya_lag1 + wa_ffb_price_imp1_dev_3pya_lag1 +
        wa_cpo_price_imp1_yoyg_3pya_lag1 + wa_cpo_price_imp1_dev_3pya_lag1 +
        wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
        n_reachable_uml_lag1 | 
        parcel_id
    }
  }
  
  ## Regressions
  poiglm_est <- fixest::feglm(specification, data = d, family = "poisson", notes = FALSE)
  
  qpoiglm_est <- fixest::feglm(specification, data = d, family = "quasipoisson", notes = FALSE)
  
  poimlm_est <- fixest::femlm(specification, data = d, family = "poisson", notes = FALSE)
  
  nb_est <- fixest::fenegbin(specification, data = d, notes = FALSE)
  
  gauglm_est <- fixest::feglm(specification, data = d, family = "gaussian", notes = FALSE)
  
  # title
  if(length(island) == 1){
    table_title <- paste0("Estimator comparisons, ",island, 
                        " catchment radius of ", catchment_radius/1000,"km, ") 
  }else{
    table_title <- paste0("Estimator comparisons, all islands, catchment radius of ", 
                          catchment_radius/1000,"km, ") 
  }
  
  # # file 
  # if(length(island) == 1){
  #   table_file <- file.path(paste0("outputs/regressions/tables/est_comparisons_",
  #                                  island,"_",
  #                                  parcel_size/1000,"km_",
  #                                  catchment_radius/1000,"km_",
  #                                  outcome_variable, 
  #                                  ".tex")) 
  # }else{
  #   table_file <- file.path(paste0("outputs/regressions/tables/est_comparisons_all_",
  #                                  parcel_size/1000,"km_",
  #                                  catchment_radius/1000,"km_",
  #                                  outcome_variable, 
  #                                  ".tex")) 
  # }


  return(etable(poiglm_est, qpoiglm_est, poimlm_est, nb_est, gauglm_est, 
         cluster = ~ parcel_id, 
         tex = TRUE,
         # file = table_file, 
         # replace = TRUE,
         title = table_title,
         subtitles = c("Poisson GLM", "Quasi-Poisson GLM", "Poisson MLM", "Negative Binomial", "Gaussian GLM"),
         family = FALSE,
         coefstat = "confint",
         sdBelow = TRUE,
         dict = TRUE))  
  
 rm(d)
}


OV <- "lucpfip_pixelcount_total"
COMMO <- "ffb"
DYN <- "both"
CR <- 1e4

# "All" islands
ISL <- c("Sumatra", "Kalimantan", "Papua")
for(CR in c(1e4, 3e4, 5e4)){
    for(COMMO in c("ffb", "cpo", "both")){
      for(DYN in c("SR", "LR", "both")){
        for(OV in c("lucpfip_pixelcount_total", "lucpfip_ha_total", 
                    "lucfip_pixelcount_30th", "lucfip_ha_30th")){
        compare_estimators(catchment_radius = CR, 
                     island = ISL, 
                     outcome_variable = OV, 
                     commo = COMMO, 
                     dynamics = DYN)
      }
    }
  }
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### COMPARE SPECIFICATIONS ALONG FULL SET OF FIXED-EFFECTS ##### 
# (With QUASIPOISSON)


catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
outcome_variable <- "lucpfip_pixelcount_total"
dynamics <- FALSE
commo <- c("ffb", "cpo")
yoyg <- FALSE
short_run <- "unt level" # from c("unt_level", "dev", "yoyg")
imp <- 1
x_pya <- 3 # 2, 3 or 4
lag_or_not <- "_lag1" # from c("_lag1", "")
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
oneway_cluster <- ~parcel_id


compare_fe_all_islands <- function(catchment_radius, # c(1e4, 3e4, 5e4)
                                   island, # c("Sumatra", "Kalimantan", "Papua", c("Sumatra", "Kalimantan", "Papua"))
                                   outcome_variable, # c("lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucfip_pixelcount_intact", "lucfip_pixelcount_degraded", "lucfip_pixelcount_total",)
                                   dynamics, # TRUE or FALSE
                                   commo, # c("ffb", "cpo", c("ffb", "cpo"))
                                   yoyg, # TRUE or FALSE 
                                   short_run, # sub or full vector from c("unt level", "yoyg", "dev") 
                                   imp, # c(1,2)
                                   x_pya, # c(2, 3, 4)
                                   lag_or_not, # c("_lag1", "")
                                   controls, # character vectors of base names of controls (don't specify in their names)
                                   weights,
                                   oneway_cluster # formula if clusters are to be specified (see argument cluster in ?fixest::etable)
                                   ){
  
  
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
  
  # list to be filled with formulae of different fixed-effect models
  fe_model_list <- list() 
  
  # character vector of fixed effect specifications to be compared in a single regression table. 
  # Should be given in a fixest format.
  if(length(island) == 1){
    fixed_effects <- c("parcel_id", 
                       "parcel_id + year", 
                       "parcel_id + province^year", 
                       "parcel_id + district^year")
  }else{
    fixed_effects <- c("parcel_id", 
                       "parcel_id + year", 
                       "parcel_id + island^year", 
                       "parcel_id + province^year", 
                       "parcel_id + district^year") 
  } 
  
  # list model formulae with different fixed effects
  for(fe in fixed_effects){
    fe_model_list[[match(fe, fixed_effects)]] <- as.formula(paste0(outcome_variable,
                                                                " ~ ",
                                                                paste0(regressors, collapse = "+"),
                                                                " + ",
                                                                paste0(paste0(controls,lag_or_not), collapse = "+"),
                                                                " | ",
                                                                fe))
  }
  
# run quasipoisson regression for each fixed-effect model 
  if(weights == TRUE){
    var_weights <- d$sample_coverage_lag1/100 
    fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                          data = d, 
                          family = "quasipoisson", 
                          notes = FALSE, 
                          weights = var_weights)
  }else{
    fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                          data = d, 
                          family = "quasipoisson", 
                          notes = FALSE)
  }

  
  
  ### Return tables
  # title
  if(length(island) == 1){
    table_title <- paste0("Fixed-effect comparisons, ",island, 
                          " catchment radius of ", catchment_radius/1000,"km ") 
  }else{
    table_title <- paste0("Fixed-effect comparisons, all islands, catchment radius of ", 
                          catchment_radius/1000,"km ") 
  }
  
  # esttable(fe_reg_list,            
  #          se = "cluster", 
  #          drop = c("own", "reachable"))
  
  ## ONE WAY CLUSTERING
  if(length(island)>1){
    etable(fe_reg_list, 
           #cluster = oneway_cluster,
           se = "cluster",
           tex = TRUE,
           # file = table_file, 
           # replace = TRUE,
           title = table_title,
           # subtitles = c("FE: grid cell", "FE: year", "FE: grid cell + year", "FE: grid cell + island*year", "FE: grid cell + province*year", "FE: grid cell + district*year"),
           family = TRUE,
           drop = c("own", "reachable"),
           coefstat = "confint",
           sdBelow = FALSE,
           yesNoFixef = "X",
           fitstat = c("sq.cor"),
           dict = TRUE)
  }else{
    etable(fe_reg_list,
           #cluster = oneway_cluster,
           se = "cluster",
           tex = TRUE,
           # file = table_file, 
           # replace = TRUE,
           title = table_title,
           # subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
           # subtitles = c("FE: grid cell", "FE: year", "FE: grid cell + year", "FE: grid cell + island*year", "FE: grid cell + province*year", "FE: grid cell + district*year"),
           family = TRUE,
           drop = c("own", "reachable"),
           coefstat = "confint",
           sdBelow = FALSE,
           yesNoFixef = "X",
           fitstat = c("sq.cor"),
           dict = TRUE)
  }
  
  
  # ## TWO WAY CLUSTERING
  # if(length(island)>1){
  #   etable(qpoiglm_est_lucfip_ufe, 
  #          qpoiglm_est_lucfip_tfe, 
  #          qpoiglm_est_lucfip_twfe, 
  #          qpoiglmspec_lucfip_iyfe, 
  #          qpoiglmspec_lucfip_dyfe, 
  #          cluster = twoway_cluster,
  #          se = "twoway",
  #          tex = TRUE,
  #          # file = table_file, 
  #          # replace = TRUE,
  #          title = table_title,
  #          subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
  #          family = FALSE,
  #          drop = c("own", "reachable"),
  #          coefstat = "confint",
  #          sdBelow = TRUE,
  #          dict = TRUE)
  # }else{
  #   etable(qpoiglm_est_lucfip_ufe, 
  #          qpoiglm_est_lucfip_tfe, 
  #          qpoiglm_est_lucfip_twfe, 
  #          qpoiglmspec_lucfip_dyfe, 
  #          cluster = twoway_cluster,
  #          se = "twoway",
  #          tex = TRUE,
  #          # file = table_file, 
  #          # replace = TRUE,
  #          title = table_title,
  #          subtitles = c("Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM", "Quasi-Poisson GLM"),
  #          family = FALSE,
  #          drop = c("own", "reachable"),
  #          coefstat = "confint",
  #          sdBelow = TRUE,
  #          dict = TRUE)
  # }
  rm(d)
}



# Run on all islands
# "All" islands
ISL <- c("Sumatra", "Kalimantan", "Papua")
CR <- 3e4
OV <- "lucpfip_pixelcount_total"
YOYG <- FALSE
DYN <- TRUE
XPYA <- 2
ISL <- "Kalimantan"

#for(YOYG in c(0, 1)){ # put this first because these are not comaparable measures and hence coeff. 
#for(CR in c(3e4, 5e4)){
  #for(SR in c("unt level", "dev")){
for(ISL in c("Sumatra", "Kalimantan")){
for(DYN in c(0,1)){
  for(XPYA in c(2, 3, 4)){
    compare_fe_all_islands(catchment_radius = CR, 
                           island = ISL,
                           outcome_variable = OV,
                           dynamics = DYN,
                           commo = c("cpo"), 
                           yoyg = YOYG,
                           short_run = "unt level", # does not matter if dynamics == FALSE
                           imp = 1,
                           x_pya = XPYA,
                           lag_or_not = "_lag1", # bien vérifier ça !  
                           controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml"),
                           weights = TRUE,
                           oneway_cluster = ~parcel_id)
  }
}
} 
#}
#}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### SPECIFICATION CHARTS / ROBUSTNESS CHECKS #####


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


# make labels for the chart
schart_labels <- list("Outcome variable:" = c("LUCPFIP", 
                                              "LUCFIP 30%"),
                      "Sample:" = c("30km catchment radius",
                                    "50km catchment radius", 
                                    "Data cleaning stronger imputations"),
                      "Distribution assumption:" = c("Poisson",
                                                     "Quasi-Poisson", 
                                                     "Negative Binomial"),
                      "Price signal averaged over:" = c("3 years", 
                                                        "4 years", 
                                                        "5 years"), 
                      "One-year lagged regressors" = "", 
                      "Fixed effects:" = c("Grid cell", 
                                           "Grid cell and year", 
                                           "Grid cell and province-year",
                                           "Grid cell and district-year"),
                      "Controls:" = c(paste0(toupper(c("ffb", "cpo")[!grepl(VAR, c("ffb", "cpo"))]), " price signal"), 
                                      "Ownership",
                                      "# reachable UML mills", 
                                      "% CPO exported", 
                                      "PYA outcome", 
                                      "PYA outcome neighbors"),
                      "Weights" = "", 
                      "Standard errors:" = c("Grid cell cluster", 
                                             "Two-way cluster")
) 



# We now produce the "data" argument for the function schart, i.e. we run regressions. 
# catchment_radius <- 3e4
# island <- "Kalimantan"
# outcome_variable <- "lucpfip_pixelcount_total"
# dynamics <- FALSE
# commo <- c("ffb", "cpo")
# yoyg <- FALSE
# short_run <- "unt level"
# imp <- 1
# distribution <- "quasipoisson"
# fixed_effects <- "parcel_id"
# x_pya <- 3
# lag_or_not <- "_lag1"
# controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
#               "n_reachable_uml")
# pya_ov <- FALSE
# weights <- TRUE
# cluster <- "cluster"
# variable <- "ffb"

gen_reg_results_in_df <- function(variable, 
                                  catchment_radius, # c(1e4, 3e4, 5e4)
                                  island, # c("Sumatra", "Kalimantan", "Papua", c("Sumatra", "Kalimantan", "Papua"))
                                  outcome_variable, # c("lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucfip_pixelcount_intact", "lucfip_pixelcount_degraded", "lucfip_pixelcount_total",)
                                  dynamics, # TRUE or FALSE
                                  commo, # c("ffb", "cpo", c("ffb", "cpo"))
                                  yoyg, # TRUE or FALSE 
                                  short_run, # sub or full vector from c("unt level", "yoyg", "dev") 
                                  imp, # c(1,2)
                                  distribution, # either "poisson", "quasipoisson", or "negbin"
                                  fixed_effects,
                                  x_pya, # c(2, 3, 4)
                                  lag_or_not, # c("_lag1", "")
                                  controls, # character vectors of base names of controls (don't specify in their names)
                                  pya_ov, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                                  weights,
                                  cluster # one of etable's "se" argument (esp. "cluster" or "twoway". Passed to se argument of etable. 
){
  
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
  
  # list model formulae with different fixed effects
  fml <- as.formula(paste0(outcome_variable,
                           " ~ ",
                           paste0(regressors, collapse = "+"),
                           " + ",
                           paste0(paste0(controls,lag_or_not), collapse = "+"),
                           " | ",
                           fixed_effects))
  
  if(distribution != "negbin"){
    # run feglm regressions with the specified family distribution
    if(weights == TRUE){
      var_weights <- d$sample_coverage_lag1/100 
      one_reg_results <- fixest::feglm(fml, 
                                       data = d, 
                                       family = distribution, 
                                       notes = TRUE, 
                                       weights = var_weights)
    }else{
      one_reg_results <- fixest::feglm(fml, 
                                       data = d, 
                                       family = distribution, 
                                       notes = TRUE)
    }
  }else{
    # run negative binomial regression with fixest::fenegbin
    # (and we cannot add weights)
    one_reg_results <- fixest::fenegbin(fml, 
                                        data = d, 
                                        notes = TRUE)
  }
  
  rm(d)
  
  # make the vector of stats needed to display in spec chart.
  reg_summary <- summary(one_reg_results, se = cluster)
  
  # if we want to compute the confidence int. within the function
  # confint <- confint(object = one_reg_results, 
  #                    parm = match(variable, commo), 
  #                    level = 0.95, 
  #                    se = cluster)
  # reg_stats <- cbind(reg_summary$coeftable[match(variable, commo),1:2], confint)
  
  reg_coeff_se <- reg_summary$coeftable[match(variable, commo), 1:2]
  
  ### make indicator variables that will be used to label specifications. 
  ind_var <- data.frame(# outcome variable
    "lucpfip_pixelcount_total" = FALSE, 
    "lucfip_pixelcount_30th" = FALSE,
    # sample
    "CR_30km" = FALSE,
    "CR_50km" = FALSE,
    "imp1" = FALSE, 
    # distribution assumptions
    "poisson" = FALSE,
    "quasipoisson" = FALSE,
    "negbin" = FALSE,
    # dependent variable
    "pya_3" = FALSE,
    "pya_4" = FALSE,
    "pya_5" = FALSE,
    "lag_or_not" = FALSE,
    # fixed effects
    "unit_fe" = FALSE,
    "tw_fe" = FALSE,
    "unit_provyear_fe" = FALSE,
    "unit_distryear_fe" = FALSE,
    # controls
    "two_commo" = FALSE, 
    "control_own" = FALSE,
    "n_reachable_uml_control" = FALSE, 
    "prex_cpo_control" = FALSE, 
    "pya_outcome_control" = FALSE,
    "pya_outcome_spatial" = FALSE,
    # weights
    "weights" = FALSE,
    # standard errors
    "oneway_cluster" = FALSE, 
    "twoway_cluster" = FALSE
  )
  
  # change the OV indicator variable to TRUE 
  ind_var[,outcome_variable] <- TRUE
  
  # sample  
  # change the CR indicator variable to TRUE for the corresponding CR
  ind_var[,grepl(paste0("CR_", catchment_radius/1000,"km"), colnames(ind_var))] <- TRUE
  
  # set the indicator variable for the data imputation 
  if(imp == 1){ind_var[,"imp1"] <- TRUE}
  
  # set the indicator variables for the distribution assumptions
  ind_var[,distribution] <- TRUE
  
  # dependent variable
  # set # year average indicator variables (+1 because we look at the total effect)
  ind_var[,grepl(paste0("pya_", x_pya+1), colnames(ind_var))] <- TRUE
  # set indicator variable for lagging or not
  if(lag_or_not == "_lag1"){ind_var[,"lag_or_not"] <- TRUE}
  
  # set indicator variables for fixed effects
  if(fixed_effects == "parcel_id"){ind_var[,"unit_fe"] <- TRUE}
  if(fixed_effects == "parcel_id + year"){ind_var[,"tw_fe"] <- TRUE}
  if(fixed_effects == "parcel_id + province_year"){ind_var[,"unit_provyear_fe"] <- TRUE}
  if(fixed_effects == "parcel_id + district_year"){ind_var[,"unit_distryear_fe"] <- TRUE}
  
  # set indicator variables for controls included
  # second commo 
  if(length(commo) == 2){ind_var[,"two_commo"] <- TRUE}
  # ownership
  if(any(grepl("pct_own", controls))){ind_var[,"control_own"] <- TRUE}
  # n_reachable_uml
  if(any(grepl("n_reachable_uml", controls))){ind_var[,"n_reachable_uml_control"] <- TRUE}
  # prex_cpo
  if(any(grepl("prex_cpo", controls))){ind_var[,"prex_cpo_control"] <- TRUE}
  # pya outcome
  if(pya_ov){ind_var[,"pya_outcome_control"] <- TRUE}
  # pya outcome spatial
  
  # weights indicator variable
  if(weights){ind_var[,"weights"] <- TRUE}
  if(distribution == "negbin"){ind_var[,"weights"] <- FALSE} # turn it back to FALSE if negbin distribution
  
  # clustering
  if(cluster == "cluster"){ind_var[,"oneway_cluster"] <- TRUE}
  if(cluster == "twoway"){ind_var[,"twoway_cluster"] <- TRUE}
  
  
  spec_df <- cbind(reg_coeff_se, ind_var)
  return(spec_df)
}

### GIVE HERE THE ISLAND AND COMMODITY FOR WHICH YOU WANT THE SPEC CHART TO BE EDITED
ISL <- c("Sumatra", "Kalimantan", "Papua")
ISL <- "Kalimantan"
VAR <- "cpo"

i <- 1
reg_stats_indvar_list <- list()

### RUN THESE FIRST LOOPS OF REGRESSIONS
for(OV in c("lucpfip_pixelcount_total", "lucfip_pixelcount_30th")){
  for(CR in c(3e4, 5e4)){
    for(IMP in c(1, 2)){
      #for(DISTR in c("quasipoisson", "negbin")){#
      for(XPYA in c(2, 3, 4)){
        for(LAG in c("_lag1", "")){
          for(FE in c("parcel_id", "parcel_id + district_year")){
            #for(WGH in c(TRUE, FALSE)){
            #for(CLST in c("cluster", "twoway")){
            reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                                island = ISL,
                                                                outcome_variable = OV,
                                                                catchment_radius = CR,
                                                                dynamics = FALSE,
                                                                commo = c("ffb", "cpo"),
                                                                yoyg = FALSE,
                                                                short_run = "unt level", # does not matter if dynamics == FALSE
                                                                imp = IMP,
                                                                distribution = "quasipoisson",
                                                                fixed_effects = FE,
                                                                x_pya = XPYA,
                                                                lag_or_not = LAG, # bien vérifier ça !
                                                                controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                             "n_reachable_uml"),
                                                                pya_ov = FALSE,
                                                                weights = TRUE,
                                                                cluster = "cluster"
            )
            i <- i+1
          }
        }
      }
    }
  }
}
#}
#}
#}

### THEN ADD RESULTS OF DEPARTURES FROM THE PREFERED SPECIFICATION

## For a different catchment radii
# (we do not do 10km because it yields too different results that make the chart hard to read)
# 
# for(FE in c("parcel_id", "parcel_id + district_year")){
# reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
#                                                     island = ISL,
#                                                     outcome_variable = "lucpfip_pixelcount_total",
#                                                     catchment_radius = 5e4,
#                                                     dynamics = FALSE,
#                                                     commo = c("ffb", "cpo"),
#                                                     yoyg = FALSE,
#                                                     short_run = "unt level", # does not matter if dynamics == FALSE
#                                                     imp = 1,
#                                                     distribution = "quasipoisson",
#                                                     fixed_effects = FE,
#                                                     x_pya = 3,
#                                                     lag_or_not = "_lag1", # bien vérifier ça !
#                                                     controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
#                                                                  "n_reachable_uml"),
#                                                     pya_ov = FALSE,
#                                                     weights = TRUE,
#                                                     cluster = "cluster"
# )
# i <- i+1
# }


# ## For different past year average
# for(XPYA in c(2, 4)){
#   for(FE in c("parcel_id", "parcel_id + district_year")){
#   reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
#                                                       island = ISL,
#                                                       outcome_variable = "lucpfip_pixelcount_total",
#                                                       catchment_radius = 3e4,
#                                                       dynamics = FALSE,
#                                                       commo = c("ffb", "cpo"),
#                                                       yoyg = FALSE,
#                                                       short_run = "unt level", # does not matter if dynamics == FALSE
#                                                       imp = 1,
#                                                       distribution = "quasipoisson",
#                                                       fixed_effects = FE,
#                                                       x_pya = XPYA,
#                                                       lag_or_not = "_lag1", # bien vérifier ça !
#                                                       controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
#                                                                    "n_reachable_uml"),
#                                                       pya_ov = FALSE,
#                                                       weights = TRUE,
#                                                       cluster = "cluster"
#   )
#   i <- i+1
# }
# }

## For poisson distribution
for(FE in c("parcel_id", "parcel_id + district_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = c("ffb", "cpo"),
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "poisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                   "n_reachable_uml"),
                                                      pya_ov = FALSE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}

## For different fixed effects
for(FE in c("parcel_id + year", "parcel_id + province_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = c("ffb", "cpo"),
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "quasipoisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                   "n_reachable_uml"),
                                                      pya_ov = FALSE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}

## CONTROLS
## Without the other commodity price signal
for(FE in c("parcel_id", "parcel_id + district_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = VAR,
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "quasipoisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                   "n_reachable_uml"),
                                                      pya_ov = FALSE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}

## Without ownership control
for(FE in c("parcel_id", "parcel_id + district_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = c("ffb", "cpo"),
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "quasipoisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("n_reachable_uml"),
                                                      pya_ov = FALSE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}

## Without n_reachable_uml
for(FE in c("parcel_id", "parcel_id + district_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = c("ffb", "cpo"),
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "quasipoisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),
                                                      pya_ov = FALSE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}


## With percentage CPO exported
for(FE in c("parcel_id", "parcel_id + district_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = c("ffb", "cpo"),
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "quasipoisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                   "n_reachable_uml",
                                                                   "wa_prex_cpo_imp1"),
                                                      pya_ov = FALSE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}

## With past year average outcome 
for(FE in c("parcel_id", "parcel_id + district_year")){
  reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                      island = ISL,
                                                      outcome_variable = "lucpfip_pixelcount_total",
                                                      catchment_radius = 3e4,
                                                      dynamics = FALSE,
                                                      commo = c("ffb", "cpo"),
                                                      yoyg = FALSE,
                                                      short_run = "unt level", # does not matter if dynamics == FALSE
                                                      imp = 1,
                                                      distribution = "quasipoisson",
                                                      fixed_effects = FE,
                                                      x_pya = 3,
                                                      lag_or_not = "_lag1", # bien vérifier ça !
                                                      controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                   "n_reachable_uml"),
                                                      pya_ov = TRUE,
                                                      weights = TRUE,
                                                      cluster = "cluster"
  )
  i <- i+1
}

## Without IBS mills sample coverage weights
for(CR in c(3e4, 5e4)){
  for(FE in c("parcel_id", "parcel_id + district_year")){
    reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                        island = ISL,
                                                        outcome_variable = "lucpfip_pixelcount_total",
                                                        catchment_radius = CR,
                                                        dynamics = FALSE,
                                                        commo = c("ffb", "cpo"),
                                                        yoyg = FALSE,
                                                        short_run = "unt level", # does not matter if dynamics == FALSE
                                                        imp = 1,
                                                        distribution = "quasipoisson",
                                                        fixed_effects = FE,
                                                        x_pya = 3,
                                                        lag_or_not = "_lag1", # bien vérifier ça !
                                                        controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                     "n_reachable_uml"),
                                                        pya_ov = FALSE,
                                                        weights = FALSE,
                                                        cluster = "cluster"
    )
    i <- i+1
  }
}

## With two way clustering
reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                    island = ISL,
                                                    outcome_variable = "lucpfip_pixelcount_total",
                                                    catchment_radius = 3e4,
                                                    dynamics = FALSE,
                                                    commo = c("ffb", "cpo"),
                                                    yoyg = FALSE,
                                                    short_run = "unt level", # does not matter if dynamics == FALSE
                                                    imp = 1,
                                                    distribution = "quasipoisson",
                                                    fixed_effects = "parcel_id + district_year",
                                                    x_pya = 3,
                                                    lag_or_not = "_lag1", # bien vérifier ça !
                                                    controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                 "n_reachable_uml"),
                                                    pya_ov = FALSE,
                                                    weights = TRUE,
                                                    cluster = "twoway"
)
i <- i+1

## Run some specifications with negbin 
# save it before negbin regressions which tend to make the R session crash.
# if(length(ISL) ==1){
# saveRDS(reg_stats_indvar_list, file.path(paste0("temp_data/reg_results/spec_chart_df_wo_negbin_",ISL,"_",VAR)))
# reg_stats_indvar_list <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_wo_negbin_",ISL,"_",VAR)))
# }else{
#   saveRDS(reg_stats_indvar_list, file.path(paste0("temp_data/reg_results/spec_chart_df_wo_negbin_all_",VAR)))
#   reg_stats_indvar_list <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_wo_negbin_all_",VAR)))
# }
# 
# i <- length(reg_stats_indvar_list)+1

for(CR in c(3e4, 5e4)){
  #for(IMP in c(1, 2)){
  #for(DISTR in c("quasipoisson", "negbin")){#
  for(XPYA in c(2, 3, 4)){
    #for(LAG in c("_lag1", "")){
    for(FE in c("parcel_id", "parcel_id + district_year")){
      reg_stats_indvar_list[[i]] <- gen_reg_results_in_df(variable = VAR, # this determines the variable the coefficient of which we are interested in.
                                                          island = ISL,
                                                          outcome_variable = "lucpfip_pixelcount_total",
                                                          catchment_radius = CR,
                                                          dynamics = FALSE,
                                                          commo = c("ffb", "cpo"),
                                                          yoyg = FALSE,
                                                          short_run = "unt level", # does not matter if dynamics == FALSE
                                                          imp = 1,
                                                          distribution = "negbin",
                                                          fixed_effects = FE,
                                                          x_pya = XPYA,
                                                          lag_or_not = "_lag1", # bien vérifier ça !
                                                          controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                                                       "n_reachable_uml"),
                                                          pya_ov = FALSE,
                                                          weights = FALSE,
                                                          cluster = "cluster"
      )
      i <- i+1
    }
  }
}
# }

# convert to dataframe to be able to chart
reg_stats_indvar <- bind_rows(reg_stats_indvar_list)

if(sum(duplicated(reg_stats_indvar))==0 & nrow(reg_stats_indvar)+1 == i){
  # save it 
  if(length(ISL) ==1){
    saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",VAR)))
  }else{
    saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_all_",VAR)))
  }
}else{print("SOMETHING WENT WRONG")}



# find position of model to highlight in original data frame
a <- reg_stats_indvar
model_idx <- a[a$lucpfip_pixelcount_total & 
                 a$CR_30km & 
                 a$imp1 &
                 a$quasipoisson &
                 a$pya_4 &
                 a$lag_or_not &
                 (a$unit_fe | a$unit_distryear_fe) &
                 a$two_commo & 
                 a$control_own & 
                 a$n_reachable_uml_control &
                 a$prex_cpo_control == FALSE &
                 a$pya_outcome_control == FALSE &
                 a$weights &
                 a$oneway_cluster, ] %>% rownames()
model_idx

schart(reg_stats_indvar, 
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
       ylab = paste0("Semi-elasticity to ", toupper(VAR), " price signal"),
       lwd.est = 5.8,
       lwd.symbol = 1
       #pch.est=20
)


### REMOVE SOME SPECIFICATIONS IF DESIRED 
ISL <- c("Sumatra", "Kalimantan", "Papua")
ISL <- "Kalimantan"
VAR <- "cpo"
if(length(ISL)==1){
  a <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",VAR)))
}else{
  a <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_all_",VAR)))
}
colnames(a)
nrow(a)
# remove regressions with: 
a <- dplyr::filter(a, !(CR_50km & negbin) & # negbin with 50km CR
                     !((pya_5 | pya_3) & negbin) & # negbin with other pya lengths than 4
                     !lucfip_pixelcount_30th & # lucfip outcome measure
                     lag_or_not & # not lagged regressors
                     imp1) # weaker imputations in data cleaning (imp2) 

# change indicator variables and labels to lighten the chart
funt <- function(x){any(x) & !all(x)} # returns TRUE iff there is any TRUE in the columns, but not only  
reg_stats_indvar_rstr <- cbind(a[,c(1,2)], select_if(a[,-c(1,2)], .predicate=funt))
colnames(reg_stats_indvar_rstr)

schart_labels_rstr <- list("Sample:" = c("30km catchment radius",
                                         "50km catchment radius"),
                           "Distribution assumption:" = c("Poisson",
                                                          "Quasi-Poisson", 
                                                          "Negative Binomial"),
                           "Price signal averaged over:" = c("3 years", 
                                                             "4 years", 
                                                             "5 years"), 
                           "Fixed effects:" = c("Grid cell", 
                                                "Grid cell and year", 
                                                "Grid cell and province-year",
                                                "Grid cell and district-year"),
                           "Controls:" = c(paste0(toupper(c("ffb", "cpo")[!grepl(VAR, c("ffb", "cpo"))]), " price signal"), 
                                           "Ownership",
                                           "# reachable UML mills", 
                                           "% CPO exported", 
                                           "PYA outcome"),
                           "Weights" = "", 
                           "Standard errors:" = c("Grid cell cluster", 
                                                  "Two-way cluster")
) 

# find position of model to highlight in original data frame
a <- reg_stats_indvar_rstr
model_idx_rstr <- a[a$CR_30km & 
                      a$quasipoisson &
                      a$pya_4 &
                      (a$unit_fe | a$unit_distryear_fe) &
                      a$two_commo & 
                      a$control_own & 
                      a$n_reachable_uml_control &
                      a$prex_cpo_control == FALSE &
                      a$pya_outcome_control == FALSE &
                      a$weights &
                      a$oneway_cluster, ] %>% rownames()
model_idx_rstr

schart(reg_stats_indvar_rstr, 
       labels = schart_labels_rstr,
       order="increasing",
       heights=c(.6,1.5),
       pch.dot=c(20,20,20,20),
       ci=c(.95),
       cex=c(0.6,0.7),
       highlight=model_idx_rstr, 
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
       ylab = paste0("Semi-elasticity to ", toupper(VAR), " price signal"),
       lwd.est = 5.8,
       lwd.symbol = 1
       #pch.est=20
)











##### ZERO INFLATED MODELLING ##### 
catchment_radius <- 3e4
outcome_variable <- "lucpfip_pixelcount_total"
x_pya <- 3
lag_or_not <- "_lag1"
fixed_effects <- "parcel_id + district^year"
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
weights <- TRUE
dynamics <- TRUE
commo <- c("ffb", "cpo")
yoyg <- FALSE
short_run <- "unt level"
imp <- 1
i <- 1

island_list <- list("Sumatra", "Kalimantan", c("Sumatra", "Kalimantan", "Papua"))
for(i in 1:length(island_list)){
  island <- island_list[[i]]
  
  # DATA FOR REGRESSIONS
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
  
  sum(d[,outcome_variable]==0)/length(d[,outcome_variable])   
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
  
  
  
  # Build the second dimension FE dummies
  d$prov_year <- paste0(d$province, "_", d$year)
  length(unique(d$year))*length(unique(d$province)) == length(unique(d$prov_year))
  
  d$distr_year <- paste0(d$district, "_", d$year)
  length(unique(d$year))*length(unique(d$district)) == length(unique(d$distr_year))
  
  # demean regressor variables by first dimensional FE : parcel_id.  
  for(var in c(outcome_variable, regressors, paste0(controls,lag_or_not))){
    unit_means <- ddply(d, "parcel_id", summarise,
                        !!as.symbol(paste0(var,"_um")) := mean(!!as.symbol(var), na.rm = TRUE))
    d <- merge(d, unit_means, by = "parcel_id")
    
    d[,paste0(var,"_dm")] <- d[,var] - d[,paste0(var,"_um")]
    #summary(d[d$parcel_id == 1525,paste0(var,"_dm")])
  }
  
  summary(d[d$parcel_id == 1525,paste0(var,"_dm")])
  View(d[!is.na(d[,paste0(var,"_um")]), c("parcel_id",var, paste0(var,"_um"), paste0(var,"_dm"))])
  
  manual_fml <- as.formula(paste0(paste0(outcome_variable,"_dm"), 
                                " ~ ", 
                                paste0(paste0(regressors, "_dm"), collapse = "+"), 
                                " + ", 
                                paste0(paste0(controls, lag_or_not, "_dm"), collapse = "+")))
  
  # does not work because the demeaned outcome variable now has negative values which is not compatible with Poisson. 
  manual_results <- fixest::feglm(manual_fml, 
                                  data = d, 
                                  family = "quasipoisson",
                                  notes = TRUE)
  
  # try with mixed model 
  mm_fml_fixed <- as.formula(paste0(outcome_variable,
                              " ~ ",
                              paste0(regressors, collapse = "+"),
                              " + ",
                              paste0(paste0(controls,lag_or_not), collapse = "+")))
# ça ne marche pas non plus.. 
  mm_results <- mixed_model(fixed = mm_fml_fixed, 
              random = ~1 | distr_year, 
              data = d,
              family = zi.poisson(), 
              zi_fixed = mm_fml_fixed)
  
  # compare with automatic fixed effects 
  
  auto_fml <- as.formula(paste0(outcome_variable,
                         " ~ ",
                         paste0(regressors, collapse = "+"),
                         " + ",
                         paste0(paste0(controls,lag_or_not), collapse = "+"),
                         " | parcel_id"))
  
  auto_results <- fixest::feglm(auto_fml, 
                                data = d, 
                                family = "quasipoisson",
                                notes = TRUE)
  
  # run quasipoisson regression for each fixed-effect model 
  if(weights == TRUE){
    var_weights <- d$sample_coverage_lag1/100 
    fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                          data = d, 
                          family = "quasipoisson", 
                          notes = FALSE, 
                          weights = var_weights)
  }else{
    fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                          data = d, 
                          family = "quasipoisson", 
                          notes = FALSE)
  }
  
  
  
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




##### INSTRUMENTAL VARIABLE ANALYSIS #####

# We use a control function approach 
# main idea: usual first stage
# extract the residuals: what the instruments + exo vars don't explain: 
# control for this endogenous variation in the 2nd stage, i.e. control for the 1st stage residuals
# as these residuals are estimated regressors, SEs are not correct. 
# thus, we estimate them with a bootstrap. 

### The instrument was constructed in add_variables.do and add_parcel_variables.R scripts.


catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
island <- "Sumatra"
outcome_variable <- "lucpfip_pixelcount_total"
dynamics <- FALSE
commo <- c("cpo")#"ffb", 
yoyg <- FALSE
short_run <- "unt level" # from c("unt_level", "dev", "yoyg")
imp <- 1
x_pya <- 3 # 2, 3 or 4
lag_or_not <- "_lag1" # from c("_lag1", "")
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml",
              paste0("wa_ffb_price_imp", imp, "_",x_pya+1,"ya"))
pya_ov <- FALSE
cluster_var <- "parcel_id"
fe <- "parcel_id + year"# needs to include year (because of the structural derivation of the first stage.
# so either "parcel_id + year" or "parcel_id + district^year"


iv <- 2 # 1,2,3,4 possible


### VARIABLES INVOLVED 

# we need a different kind of specifications as elsewhere, because we cannot have LR measures (X ya variables)
# as instrumented variables. 
endo_var <- paste0("wa_",commo,"_price_imp",imp,lag_or_not)

# add pya outcome variable 
if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}

# lag controls or not
if(lag_or_not=="_lag1"){controls <- paste0(controls,lag_or_not)}

instrument <- paste0("iv",iv,"_imp",imp,lag_or_not)


### FORMULAE

fml_1st <- as.formula(paste0(endo_var,
                             " ~ ", instrument,
                             " + ",
                             paste0(controls, collapse = "+"),
                             " | ", 
                             fe))

fml_2nd <- as.formula(paste0(outcome_variable,
                             " ~ ", endo_var,
                             " + est_res_1st +", 
                             paste0(controls, collapse = "+"),
                             " | ",
                             fe))


### DATA
d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                              parcel_size/1000,"km_",
                              catchment_radius/1000,"CR.rds")))


## Keep observations that: 

# - are in island or group of islands of interest
d <- d[d$island %in% island,]

# - were covered with some of the studied forest in 2000
if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
  d <- d[d$any_pfc2000_total,]
}
if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
   outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
   outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
  d <- d[d$any_fc2000_30th,]
}


# - are kept by fixest, i.e. that have no NA on any of the variables used (in either stage)
filter_vec <- base::rowSums(!is.na(d[,c("parcel_id", "year", outcome_variable, endo_var, controls, instrument)]))
filter_vec <- filter_vec == length(c("parcel_id", "year", outcome_variable, endo_var, controls, instrument))
d_nona <- d[filter_vec, c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument)]
anyNA(d_nona)

# - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
d_clean <- d_nona[-obs2remove(fml_2nd, d_nona, family = "poisson"),]

#rm(d, d_nona)


### COEFFICIENTS 

# # 1st stage 
est_1st <- fixest::feols(fml_1st,
                         data = d_clean)
# save estiamted residuals
d_clean$est_res_1st <- est_1st$residuals
# 2nd stage
est_2nd <- fixest::feglm(fml_2nd,
                         data = d_clean,
                         family = "quasipoisson")

# compare with standard 2sls
d_clean$est_fit_1st <- est_1st$fitted.values
std_fml_2nd <- as.formula(paste0(outcome_variable,
                                 " ~ est_fit_1st +",
                                 paste0(controls, collapse = "+"),
                                 " | ",
                                 fe))
# 2nd stage
std_est_2nd <- fixest::feglm(std_fml_2nd,
                         data = d_clean,
                         family = "quasipoisson")

# and with non-IV and same RHS
direct_fml <- as.formula(paste0(outcome_variable,
                                 " ~ ", endo_var, 
                                 " + ",
                                 paste0(controls, collapse = "+"),
                                 " | ",
                                 fe))
dir_est_2nd <- fixest::feglm(direct_fml,
                             data = d_clean,
                             family = "quasipoisson")




# library(pglm)
# fml_2nd_std <- as.formula(paste0(outcome_variable,
#                              " ~ ", endo_var,
#                              " + est_res_1st +", 
#                              paste0(controls, collapse = "+")))
# pglm_est_2nd <- pglm(fml_2nd_std, 
#                   data = d_clean_boot,
#                   index = c("parcel_id", "year"),
#                   family = "poisson", 
#                   effect = "twoways", 
#                   model = "within")
# 
# 
# bs_pglm <- sandwich::vcovBS(pglm_est_2nd, )


### CONTROL FUNCTION APPROACH WITH BOOTSTRAP

## Design the bootstrap sampling function 

# get the different cluster sizeS. This is necessary to cluster bootstraping with clusters of different sizes. 
# methodology comes from:
# https://stats.stackexchange.com/questions/202916/cluster-boostrap-with-unequally-sized-clusters/202924#202924
# this is a vector of same length as cluster_names. Values are the number of observations (the size) per cluster. 
# sizes <- sapply(cluster_names, FUN = function(name){
#                                       sum(as.character(d_clean[,cluster_var]) == name)
#                                      })
# sizes2 <- ddply(d_clean, cluster_var, summarise,
#               sizes = sum(duplicated(as.symbol(cluster_var))))
sizes <- table(d_clean[,cluster_var])
u_sizes <- sort(unique(sizes))

cl_names <- list()
n_clusters <- list()
for(s in u_sizes){
  # names and numbers of clusters of size s
  cl_names[[s]] <- names(sizes[sizes == s])
  n_clusters[[s]] <- length(cl_names[[s]])
}

sum(unlist(n_clusters)) == length(unique(d_clean$parcel_id))

# function that will tell boot::boot how to sample data at each replicate of statistic
ran.gen_cluster <- function(original_data, arg_list){

  
  
  cl_boot_dat <- NULL
  
  for(s in arg_list[["unique_sizes"]]){
    # resample cluster names 
    sample_cl_s <- sample(arg_list[["cluster_names"]][[s]], 
                          arg_list[["number_clusters"]][[s]], 
                          replace = TRUE) 
    sample_cl_s_tab <- table(sample_cl_s)
    
    for(n in 1:max(sample_cl_s_tab)){
      # vector to select obs. that are within the sampled clusters. 
      sel <- as.character(original_data[,arg_list[["cluster_variable"]]]) %in% names(sample_cl_s_tab[sample_cl_s_tab %in% n])
      
      # replicate the selected data n times
      clda <- original_data[sel,][rep(seq_len(nrow(original_data[sel,])), n), ]
      
      # we need to give a new cluster identifier during the resampling, otherwise a cluster sampled more than once 
      # will be "incorrectly treated as one large cluster rather than two distinct cluster" (by the fixed effects) (Cameron 2015)
      
      #identify row names without periods, and add ".0" 
      row.names(clda)[grep("\\.", row.names(clda), invert = TRUE)] <- paste0(grep("\\.", row.names(clda), invert = TRUE, value = TRUE),".0")
      # add the suffix due to the repetion after the existing cluster identifier. 
      clda[,arg_list[["cluster_variable"]]] <- paste0(clda[,arg_list[["cluster_variable"]]], 
                                                      sub(".*\\.","_",row.names(clda)))
      
      cl_boot_dat <- cl_boot_dat %>% rbind(clda) 
    }
  }  
  
  #test that the returned data are the same dimension as input
  # dim(cl_boot_dat)
  # dim(original_data)
  # dim(cl_boot_dat) == dim(d_clean)
  # 
  # # test new clusters are not duplicated (correct if anyDuplicated returns 0)
  # base::anyDuplicated(cl_boot_dat[,c("parcel_id","year")])

  return(cl_boot_dat)
}

# test the function
# test_boot_d <- ran.gen_cluster(original_data = d_clean, 
#                                arg_list = list(cluster_variable = cluster_var, 
#                                                unique_sizes = u_sizes, 
#                                                cluster_names = cl_names, 
#                                                number_clusters = n_clusters))
# dim(test_boot_d)
# dim(d_clean)

ctrl_fun_endo <- function(myfun_data, fsf, ssf){
  
  # 1st stage 
  est_1st <- fixest::feols(fsf, 
                           data = myfun_data)
  # save estiamted residuals
  myfun_data$est_res_1st <- est_1st$residuals
  # 2nd stage
  est_2nd <- fixest::feglm(ssf, 
                           data = myfun_data, 
                           family = "poisson")
  # statistics we want to evaluate the variance of:
  return(est_2nd$coefficients)
}

# Now see if we can replicate the boostrap SE. 
bootstraped <- boot(data = d_clean, 
                    statistic = ctrl_fun_endo, # 2 first arguments do not need to be called.
                    # the first one, arbitrarily called "myfun_data" is passed the previous "data" argument 
                    fsf = fml_1st,
                    ssf = fml_2nd,
                    ran.gen = ran.gen_cluster,
                    mle = list(cluster_variable = cluster_var, 
                               unique_sizes = u_sizes, 
                               cluster_names = cl_names, 
                               number_clusters = n_clusters),
                    sim = "parametric",
                    R = 200)

norm.ci(bootstraped, index = 1) # equivalent to : 
norm.ci(t0 = bootstraped$t0[1],
        t = bootstraped$t[,1])

# replicating norm.ci (bias-corrected)
i <- match(endo_var, names(bootstraped$t0))
2*bootstraped$t0[i]-mean(bootstraped$t[,i])-qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,i]))
2*bootstraped$t0[i]-mean(bootstraped$t[,i])+qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,i]))

sandwich::vcovCL(bootstraped, cluster = ~parcel_id)

# non bias corrected CI
bootstraped$t0[1] - qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,1]))
bootstraped$t0[1] + qnorm((1+0.95)/2)*sqrt(var(bootstraped$t[,1]))

# replicating Stata (non bias-corrected): 
bootstraped$t0[1] - qnorm((1+0.95)/2)*0.0281059 
bootstraped$t0[1] + qnorm((1+0.95)/2)*0.0281059 























##### REGRESSION TABLES ##### 
# We display grid cell FE and grid cell + district_year FE 
# We separate dynamics = FALSE and = TRUE in two tables, bc difficult to read otherwise.
# Group three island groups in one table ? (6 columns) --> addresses the objective of inter-groups comparisons. 
# Group weighted vs. not weighted: there are some != in Kalimantan. But it's not worth taking so much space.

# Later we might want to add "spatial model" and IV model columns. 
# Or group negbin (which is also not weighted) --> not too different, not enough justified. 

make_reg_tables <- function(island,
                            catchment_radius, # c(1e4, 3e4, 5e4)
                                  outcome_variable, # c("lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucfip_pixelcount_intact", "lucfip_pixelcount_degraded", "lucfip_pixelcount_total",)
                                  dynamics, # TRUE or FALSE
                                  commo, # c("ffb", "cpo", c("ffb", "cpo"))
                                  yoyg, # TRUE or FALSE 
                                  short_run, # sub or full vector from c("unt level", "yoyg", "dev") 
                                  imp, # c(1,2)
                                  distribution, # either "poisson", "quasipoisson", or "negbin"
                                  fe_set,
                                  x_pya, # c(2, 3, 4)
                                  lag_or_not, # c("_lag1", "")
                                  controls, # character vectors of base names of controls (don't specify in their names)
                                  pya_ov, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                                  weights,
                                  cluster # one of etable's "se" argument (esp. "cluster" or "twoway". Passed to se argument of etable. 
){
  
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
  
  # list to be filled with formulae of different fixed-effect models
  fe_model_list <- list() 
  fixed_effects <- fe_set 

  # list model formulae with different fixed effects
  for(fe in fixed_effects){
    fe_model_list[[match(fe, fixed_effects)]] <- as.formula(paste0(outcome_variable,
                                                                   " ~ ",
                                                                   paste0(regressors, collapse = "+"),
                                                                   " + ",
                                                                   paste0(paste0(controls,lag_or_not), collapse = "+"),
                                                                   " | ",
                                                                   fe))
  }
  
  # run regression for each fixed-effect model 
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    if(weights == TRUE){
      var_weights <- d$sample_coverage_lag1/100 
      fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                            data = d, 
                            family = distribution, 
                            notes = TRUE, 
                            weights = var_weights)
    }else{
      fe_reg_list <- lapply(fe_model_list, fixest::feglm,
                            data = d, 
                            family = distribution, 
                            notes = TRUE)
    }
  }else{ # no weights allowed in negative binomial
    fe_reg_list <- lapply(fe_model_list, fixest::fenegbin,
                          data = d, 
                          family = distribution, 
                          notes = TRUE)
  }
  
  rm(d)
  return(fe_reg_list)
}
# f_test <- function(x){
#   list_test <- list(paste0(x,1), paste0(x,2))
#   return(list_test)
# }
# lapply(island_list, f_test) %>% str() %>%  length()
island_list <- list(c("Sumatra", "Kalimantan", "Papua"), "Sumatra", "Kalimantan")
## WITHOUT dynamics
result_list_nodyn <- lapply(island_list, make_reg_tables, 
                            outcome_variable = "lucpfip_pixelcount_total",
                            catchment_radius = 3e4, 
                            dynamics = FALSE,
                            commo = c("ffb", "cpo"),
                            yoyg = FALSE, 
                            short_run = "unt level", 
                            imp = 1,
                            distribution = "quasipoisson",
                            fe_set = c("parcel_id", "parcel_id + district_year"),
                            x_pya = 4,
                            lag_or_not = "_lag1",
                            controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                         "n_reachable_uml"),
                            pya_ov = FALSE, 
                            weights = TRUE,
                            cluster = "cluster") %>% unlist(recursive = FALSE)
## WITH dynamics
result_list_dyn <- lapply(island_list, make_reg_tables, 
                            outcome_variable = "lucpfip_pixelcount_total",
                            catchment_radius = 3e4, 
                            dynamics = TRUE,
                            commo = c("ffb", "cpo"),
                            yoyg = FALSE, 
                            short_run = "unt level", 
                            imp = 1,
                            distribution = "quasipoisson",
                            fe_set = c("parcel_id", "parcel_id + district_year"),
                            x_pya = 4,
                            lag_or_not = "_lag1",
                            controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                                         "n_reachable_uml"),
                            pya_ov = FALSE, 
                            weights = TRUE,
                            cluster = "cluster") %>% unlist(recursive = FALSE)


# for(ISL in island_list){
#   results[[match(ISL, island_list)]] <- make_reg_tables()
# }

### Return tables
## WITHOUT time distinction between short and long term price signals 
# title
  table_title_nodyn <- paste0("LUCPFIP semi-elasticities to price signals") 
# LateX table
  etable(result_list_nodyn, 
         #cluster = oneway_cluster,
         se = "cluster",
         tex = TRUE,
         # file = table_file, 
         # replace = TRUE,
         title = table_title_nodyn,
         subtitles = c("All islands", "All islands", "Sumatra", "Sumatra", "Kalimantan", "Kalimantan"),
         family = TRUE,
         drop = c("own", "reachable"),
         coefstat = "confint",
         sdBelow = FALSE,
         yesNoFixef = "X",
         fitstat = c("sq.cor"),
         dict = TRUE, 
         powerBelow = -7)
  
## WITH the distinction between short and long term price signals
  # title
  table_title_dyn <- paste0("LUCPFIP semi-elasticities to short and long term price signals") 
  # LateX table
  etable(result_list_dyn, 
         #cluster = oneway_cluster,
         se = "cluster",
         tex = TRUE,
         # file = table_file, 
         # replace = TRUE,
         title = table_title_dyn,
         subtitles = c("All islands", "All islands", "Sumatra", "Sumatra", "Kalimantan", "Kalimantan"),
         family = TRUE,
         drop = c("own", "reachable"),
         coefstat = "confint",
         sdBelow = FALSE,
         yesNoFixef = "X",
         fitstat = c("sq.cor"),
         dict = TRUE, 
         powerBelow = -7)












##### COUNTERFACTUALS ##### 
  
  ### ISLAND & DISTRICT SHAPES
  island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
  names(island_sf)[names(island_sf)=="island"] <- "shape_des"
  
  island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)
  
  district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
  district_names <- read.dta13(file.path("temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta"))
  district_names <- district_names[!duplicated(district_names$bps_),]
  district_sf$d__2000 <- district_sf$d__2000 %>% as.character()
  district_names$bps_ <- district_names$bps_ %>% as.character()
  district_sf <- left_join(x = district_sf, y = district_names[,c("name_", "bps_")], 
                           by = c("d__2000" = "bps_"),
                           all = FALSE, all.x = FALSE, all.y = FALSE)
  district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)  

outcome_variable = "lucpfip_pixelcount_total"
dynamics = FALSE
commo = c("ffb", "cpo")#
yoyg = FALSE
short_run = "unt level"
imp = 1
distribution = "quasipoisson"
fe = c("parcel_id + district_year")
x_pya = 3
lag_or_not = "_lag1"

pya_ov = FALSE 
weights = FALSE
cluster = "cluster"

#i <- 1

island_list <- list("Sumatra", "Kalimantan")#, c("Sumatra", "Kalimantan", "Papua")

elmts <- rep(list(NA), length(island_list))
names(elmts) <- c("Sumatra", "Kalimantan")#, "All"

attributes <- c("est_results",
                "distr_lvl", "distr_cf")

for(i in 1:length(island_list)){
  
  controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
             "n_reachable_uml")
  
  island <- island_list[[i]]
  
  elmts[[i]] <- list()
  length(elmts[[i]]) <- length(attributes)
  names(elmts[[i]]) <- attributes
  
  ### VARIABLES INVOLVED 
  
  # variable of interests (called regressors pour l'instant)
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
  
  # if, on the other hand, we want to identify the short run effects, controlling for the long run's. 
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
  
  # lag controls or not
  if(lag_or_not=="_lag1"){controls <- paste0(controls,lag_or_not)}
  # add rspo_cert which is not to be lagged. 
  # controls <- c(controls, "rspo_cert")
  
  ### DATA FOR REGRESSIONS
  if(island == "Sumatra"){
    catchment_radius <- 3e4
  }else{
    catchment_radius <- 5e4}

  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.rds")))
  
  ## Keep observations that: 

  # - are in island or group of islands of interest
  d <- d[d$island %in% island,]
  
  # - were covered with some of the studied forest in 2000
  if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
    d <- d[d$any_pfc2000_total,]
  }
  if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
     outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
     outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
    d <- d[d$any_fc2000_30th,]
  }

  # - are not in an intended RSPO-certified supply base
  d <- d[d$rspo_cert==FALSE,]
  
  # - are kept by fixest, i.e. that have no NA on any of the variables used (in either stage)
  used_vars <- c("parcel_id", "year", "district", "island", 
                 "n_reachable_ibsuml_lag1", "sample_coverage_lag1", 
                 outcome_variable, regressors, controls)
  
  filter_vec <- base::rowSums(!is.na(d[,used_vars]))
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d[filter_vec, used_vars]
  if(anyNA(d_nona)){stop()}
  
  # - sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
  # probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
  # just remove them if any, so that there is no bug. 
  
  if(weights == TRUE){
    d_nona <- d_nona[d_nona$sample_coverage_lag1!=0,]
  }
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_nona$district_year <- paste0(d_nona$district,"_",d_nona$year)
  
  d_clean <- d_nona[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                                      d_nona, 
                                      family = "poisson"),]
    
  # note that obs2remove has to be done at the end, otherwise some units (along fe dimension)
  # can become "only-zero-outcome" units after other data removal.  

  ### FORMULA
  # list to be filled with formulae of different fixed-effect models
  
  fe_model <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                paste0(regressors, collapse = "+"),
                                " + ",
                                paste0(controls, collapse = "+"),
                                " | ",
                                fe))
  ### REGRESSIONS
  
  # run regression for each fixed-effect model 
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    if(weights == TRUE){
      var_weights <- d_clean$sample_coverage_lag1/100 
      fe_reg <- fixest::feglm(fe_model,
                            data = d_clean, 
                            family = distribution, 
                            notes = TRUE, 
                            weights = var_weights)
      
    }else{
      fe_reg <- fixest::feglm(fe_model,
                            data = d_clean, 
                            family = distribution, 
                            notes = TRUE)
    }
  }else{ # no weights allowed in negative binomial
    fe_reg <- fixest::fenegbin(fe_model,
                          data = d_clean, 
                          family = distribution, 
                          notes = TRUE)
  }

  
  # store it 
  elmts[[i]][["est_results"]] <- fe_reg
  
  # save coefficients 
  coeff_ffb <- fe_reg$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
  coeff_cpo <- fe_reg$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
  
  # save fitted values and linear predictors (which include FEs)
  d_clean$fitted.values <- fe_reg$fitted.values
  d_clean$linear.predictors <- fe_reg$linear.predictors
  
  ### PASSING ON TO DISTRICTS
  # dataframe of districts present in estimating sample d_clean ONLY
  distr_lvl <- data.frame(d_clean$district) %>% unique()
  colnames(distr_lvl) <- "district"
  
  
  ### AGGREGATION FACTOR
  # now we need to compute each district's number of obs.  
  # we do not only want those used in the estimation, as many were not included because of missing regressors
  # either because away from sample mills, or because at least one regressor was NA.
  # So we want to count how many grid cells were within the catchment radius of a UML mill each year, 
  # and sum over years.  
  d_aggr <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_reachable_geovars_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"km_UML_CR.rds")))
  
  ## Keep observations that: 
  # - are between 2001 and 2015
  d_aggr <- d_aggr[d_aggr$year<2016,]
  
  # - are in island or group of islands of interest
  d_aggr <- d_aggr[d_aggr$island %in% island,]
  # here it's 1172820 with the three islands
  
  # - are within the catchment radius of an UML mill this year
  # recall that this catchment radius is 30km if island is Sumatra and 50km if Kalimantan
  d_aggr <- d_aggr[d_aggr$n_reachable_uml > 0,] #933843
  # now should not be balanced anymore
  #length(unique(d_aggr$parcel_id))*length(unique(d_aggr$year))!=nrow(d_aggr)
  
  # - were covered with some of the studied forest in 2000
  # this we do not do as for regular data frame, but all these obs. are removed with 
  # obs2remove that removes those with only 0 across years. 
  d_aggr <- d_aggr[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                               data = d_aggr, 
                               family = "poisson"),] # 126082
  
  # count the number of observations per district
  # there are obs. from d_aggr that are in districts that are not in d_clean, therefore 
  # it's important that x = distr_n_parcels_within_umlcr and all.x = TRUE 
  # so that we also estimate lucfp in these districts.
  distr_n_parcels_within_umlcr <- ddply(d_aggr, "district", summarise,
                                        n_parcels_within_umlcr = length(parcel_id))
  distr_lvl <- merge(distr_n_parcels_within_umlcr, distr_lvl, 
                     by = "district",
                     all = TRUE)

  # these guys don't have an estimate of district level fitted value. 
  # give them the island level one. 
  # we are in an aggregating exercise anyways, and if we assume that our sample was randomly selected
  # among all these obs., then it makes sense to impute them the sample average. 
  
  # so put this n_parcels counting before the fitted values. Compute the average after each ddply above.
  # and replace NAs with it. 
  
  sum(distr_n_parcels_within_umlcr$n_parcels_within_umlcr)
  sum(distr_lvl$n_parcels_within_umlcr)
  
  ### FITTED VALUES 

  ## FITTED VALUE AT THE AVERAGE 
  fv_at_distr_avg <- ddply(d_clean, "district", summarise, 
                        fv_at_distr_avg = mean(linear.predictors))
  fv_at_distr_avg$fv_at_distr_avg <- exp(fv_at_distr_avg$fv_at_distr_avg)
  
  distr_lvl <- merge(distr_lvl, fv_at_distr_avg, by = "district", all = TRUE)
  
  # replace NAs with island average
  distr_lvl$fv_at_distr_avg <- replace(x = distr_lvl$fv_at_distr_avg, 
                                       list = is.na(distr_lvl$fv_at_distr_avg), 
                                       values = mean(fv_at_distr_avg$fv_at_distr_avg))
  
  distr_lvl$fv_at_distr_avg_ha <- distr_lvl$fv_at_distr_avg*(27.8*27.6)/(1e4)
  
  # One could also do it this way, equivalent, less dependent on package fixest, but longer. 
  #  covariates <- c(regressors, controls)
  # covariate_distr_avg <- distr_lvl
  # 
  # for(i in 1:length(covariates)){
  #   covariate_distr_avg_i <- ddply(d_clean, "district", summarise,
  #                                  !!as.symbol(covariates[i]) := mean(!!as.symbol(covariates[i])))   
  #   
  #   covariate_distr_avg <- merge(covariate_distr_avg, covariate_distr_avg_i, by = "district")
  # }
  # 
  # #mean(d_clean[d_clean$district == "Kota Pekan Baru", covariates[6]])
  # covariate_distr_avg <- as.matrix(dplyr::select(covariate_distr_avg, -district))
  # 
  # # check that covariates and coefficients are in the same order 
  # if(!all.equal(colnames(covariate_distr_avg), names(fe_reg_list[[1]]$coefficients)))
  #   {stop()}
  # 
  # Xb_at_distr_avg <- covariate_distr_avg %*% fe_reg_list[[1]]$coefficients
  # distr_lvl <- cbind(distr_lvl, Xb_at_distr_avg)
  # 
  # # add FEs
  # distr_avg_fe <- ddply(d_clean, "district", summarise, 
  #                      sumFE_davg = mean(sumFE))
  # distr_lvl <- merge(distr_lvl, distr_avg_fe, by = "district")
  # 
  # # take the exp of the sum of linear predictors and FE, and convert from pixels to hectares. 
  # distr_lvl$fv_atavg <- exp(distr_lvl$Xb_at_distr_avg + distr_lvl$sumFE_davg)
  # distr_lvl$fv_atavg_ha <- distr_lvl$fv_atavg*(27.8*27.6)/(1e4)
  
  ## FITTED VALUE AT THE MEDIAN 
  fv_at_distr_med <- ddply(d_clean, "district", summarise, 
                        fv_at_distr_med = median(linear.predictors))
  fv_at_distr_med$fv_at_distr_med <- exp(fv_at_distr_med$fv_at_distr_med)
  
  distr_lvl <- merge(distr_lvl, fv_at_distr_med, by = "district", all = TRUE)
  
  # replace NAs with island median
  distr_lvl$fv_at_distr_med <- replace(x = distr_lvl$fv_at_distr_med, 
                                       list = is.na(distr_lvl$fv_at_distr_med), 
                                       values = median(fv_at_distr_med$fv_at_distr_med))
  distr_lvl$fv_at_distr_med_ha <- distr_lvl$fv_at_distr_med*(27.8*27.6)/(1e4)
  
  
  ## AVERAGE FITTED VALUE 
  distr_avg_fv <- ddply(d_clean, "district", summarise, 
                        distr_avg_fv = mean(fitted.values))
  
  distr_lvl <- merge(distr_lvl, distr_avg_fv, by = "district", all = TRUE)
  
  # replace NAs with island average
  distr_lvl$distr_avg_fv <- replace(x = distr_lvl$distr_avg_fv, 
                                       list = is.na(distr_lvl$distr_avg_fv), 
                                       values = mean(distr_avg_fv$distr_avg_fv))
  distr_lvl$distr_avg_fv_ha <- distr_lvl$distr_avg_fv*(27.8*27.6)/(1e4)
  
  ## MEDIAN FITTED VALUE 
  distr_med_fv <- ddply(d_clean, "district", summarise, 
                        distr_med_fv = median(fitted.values))
  
  distr_lvl <- merge(distr_lvl, distr_med_fv, by = "district", all = TRUE)
  
  # replace NAs with island average
  distr_lvl$distr_med_fv <- replace(x = distr_lvl$distr_med_fv, 
                                    list = is.na(distr_lvl$distr_med_fv), 
                                    values = median(distr_med_fv$distr_med_fv))
  distr_lvl$distr_med_fv_ha <- distr_lvl$distr_med_fv*(27.8*27.6)/(1e4)
  
  
  ## AVERAGE TRUE LUCPFIP 
  distr_avg_y <- ddply(d_clean, "district", summarise, 
                        lucpfip_pixelcount_total_davg = mean(lucpfip_pixelcount_total))
  
  distr_lvl <- merge(distr_lvl, distr_avg_y, by = "district", all = TRUE)
  # replace NAs with island average
  distr_lvl$lucpfip_pixelcount_total_davg <- replace(x = distr_lvl$lucpfip_pixelcount_total_davg, 
                                                     list = is.na(distr_lvl$lucpfip_pixelcount_total_davg), 
                                                     values = mean(distr_avg_y$lucpfip_pixelcount_total_davg))
  distr_lvl[,"lucpfip_total_davg_ha"] <- distr_lvl[,"lucpfip_pixelcount_total_davg"]*(27.8*27.6)/(1e4)
  
  ## MEDIAN TRUE OUTCOME
  # much lower, given the skewed distribution of outcome. 
  distr_med_y <- ddply(d_clean, "district", summarise, 
                       lucpfip_pixelcount_total_dmed = median(lucpfip_pixelcount_total))
  
  distr_lvl <- merge(distr_lvl, distr_med_y, by = "district", all = TRUE)
  # replace NAs with island average
  distr_lvl$lucpfip_pixelcount_total_dmed <- replace(x = distr_lvl$lucpfip_pixelcount_total_dmed, 
                                                     list = is.na(distr_lvl$lucpfip_pixelcount_total_dmed), 
                                                     values = median(distr_med_y$lucpfip_pixelcount_total_dmed))
  
  distr_lvl[,"lucpfip_total_dmed_ha"] <- distr_lvl[,"lucpfip_pixelcount_total_dmed"]*(27.8*27.6)/(1e4)
  
  # the best "fit" seems to be the average fitted value. 
  
  # store it 
  elmts[[i]][["distr_lvl"]] <- distr_lvl
  
 
  ### COUNTERFACTUALS
  ## compute counterfactuals at the district level
  distr_cf <- distr_lvl[,c("district", 
                           "lucpfip_total_davg_ha", 
                           "distr_avg_fv_ha", "distr_med_fv_ha",
                           "fv_at_distr_avg_ha", "fv_at_distr_med_ha",
                           "n_parcels_within_umlcr")]
  
  # totals are counted in KHA, hence the /1e3

  # without tax
  
  distr_cf$total_baseline_true_kha <- (1/1e3)*distr_cf[,"n_parcels_within_umlcr"]*distr_cf[,"lucpfip_total_davg_ha"]
  distr_cf$total_baseline_avg_kha <- (1/1e3)*distr_cf[,"n_parcels_within_umlcr"]*distr_cf[,"distr_avg_fv_ha"]
  distr_cf$total_baseline_atavg_kha <- (1/1e3)*distr_cf[,"n_parcels_within_umlcr"]*distr_cf[,"fv_at_distr_avg_ha"]
  distr_cf$total_baseline_atmed_kha <- (1/1e3)*distr_cf[,"n_parcels_within_umlcr"]*distr_cf[,"fv_at_distr_med_ha"]
  
  #sum(distr_cf$total_baseline_avg_kha) #500 kha, in 30km in Sumatra, so rather close to what's in table 2. 
  
  # with counterfactual market instruments (taxes)
  # ffb
  for(tax in c(10, 50, 100)){
    distr_cf[,paste0("avg_cf_ffb_tax",tax,"_ha")] <- distr_cf$distr_avg_fv_ha*exp(coeff_ffb*(-tax))
    distr_cf[,paste0("total_cf_ffb_tax",tax,"_kha")] <- (1/1e3)*distr_cf[,"n_parcels_within_umlcr"]*distr_cf[,paste0("avg_cf_ffb_tax",tax,"_ha")]
  }
  # cpo
  for(tax in c(50, 100, 200)){
    distr_cf[,paste0("avg_cf_cpo_tax",tax,"_ha")] <- distr_cf$distr_avg_fv_ha*exp(coeff_cpo*(-tax))
    distr_cf[,paste0("total_cf_cpo_tax",tax,"_kha")] <- (1/1e3)*distr_cf[,"n_parcels_within_umlcr"]*distr_cf[,paste0("avg_cf_cpo_tax",tax,"_ha")]
  }

  # compute differences
  # ffb
  for(tax in c(10, 50, 100)){
    distr_cf[,paste0("diff_ffb_tax",tax,"_kha")] <- distr_cf$total_baseline_avg_kha - distr_cf[,paste0("total_cf_ffb_tax",tax,"_kha")]
  }
  # cpo
  for(tax in c(50, 100, 200)){
    distr_cf[,paste0("diff_cpo_tax",tax,"_kha")] <- distr_cf$total_baseline_avg_kha - distr_cf[,paste0("total_cf_cpo_tax",tax,"_kha")]
  }
  
  #How much should be compensated (lhs) for a y avoided lucfip
  
  #distr_cf$lucpfip_total_davg_ha = tax*coeff*fv
  distr_cf$tax_ZD = distr_cf$lucpfip_total_davg_ha/(coeff_cpo*distr_cf$distr_avg_fv_ha)
  all.equal(distr_cf$lucpfip_total_davg_ha, distr_cf$distr_avg_fv_ha)
  all.equal(mean(d_clean$lucpfip_pixelcount_total), mean(d_clean$fitted.values))

  # this number can be interpreted (at the island level/ or other level where the coeff is estimated)
  # as the premium that should have been transfered to all palm growers, through the price of cpo, 
  # to compensate them for an interdiction of deforestation. 
  
  # ça fait presque un million d'hectares, c'est pas encore assez par rapport au total accumulé sur la
  # période et sur les 3 îles... trouver pourquoi 
  # --> BECAUSE we do not add the estimate for all the districts that are not represented in our 
  # estimating sample d_clean. 

  
  # add the shapes to these districts
  distr_cf <- merge(distr_cf, district_sf[,"name_"], 
                                    by.x = "district",
                                    by.y = "name_", 
                                    all = FALSE) %>% st_as_sf()
  
  # store it 
  elmts[[i]][["distr_cf"]] <- distr_cf
  
  rm(d, d_nona, d_clean, d_aggr)
}


# check total areas: 
# ~514 kha, in 30km in Sumatra, so rather close to what's in table 2. 
sum(elmts[["Sumatra"]][["distr_cf"]]$total_baseline_true_kha) 
sum(elmts[["Sumatra"]][["distr_cf"]]$total_baseline_avg_kha) 
sum(elmts[["Sumatra"]][["distr_cf"]]$total_baseline_atavg_kha) 
sum(elmts[["Sumatra"]][["distr_cf"]]$total_baseline_atmed_kha) 
# ~908 kha, in 50km in Kalimantan, so rather not what's in table 2. 
sum(elmts[["Kalimantan"]][["distr_cf"]]$total_baseline_avg_kha) 
elmts[["Sumatra"]][["est_results"]]
elmts[["Kalimantan"]][["est_results"]]

# number of districts involved
elmts[["Kalimantan"]][["distr_lvl"]]$district %>% length() # 10 in d_clean but 28 in d_aggr
elmts[["Sumatra"]][["distr_lvl"]]$district %>% length() # 34 in d_clean but 43 in d_aggr




#### MAP ####

# bind df to plot
distr_cf <- rbind(elmts[["Sumatra"]][["distr_cf"]],
                  elmts[["Kalimantan"]][["distr_cf"]])

# settings for legend 
summary(distr_cf$total_baseline_avg_kha)

### MAP
# if we are to plot the doistricts' accumulated lucpfip 
max_val <- max(distr_cf$total_baseline_avg_kha)+1
# bins <- seq(from = 40, to = 180, by = 20)
cb <- colorNumeric("BuPu", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                   domain = 0:max_val,
                   na.color = "red", 
                   reverse = FALSE)
distr_cf%>% 
  leaflet() %>% 
  addTiles()%>%
  #addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  #addProviderTiles(providers$Stamen.TonerBackground) %>% 
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>% 
  addPolygons(fillColor = ~cb(distr_cf$total_baseline_avg_kha),
              label = ~total_baseline_avg_kha,
              #color = "#444444",
              stroke = FALSE,
              weight = 1,
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 1) %>% 
  addLegend(pal = cb,  
            values = distr_cf$total_baseline_avg_kha, 
            bins = 3, opacity = 0.5,
            labFormat = labelFormat(suffix = " kha"),
            title = "Predicted deforestation <br/> 2001-2015, <br/> tax = 0 (baseline)",
            position = "bottomright")


# do not do a new colorNumeric palette, we want all maps to be plotting along the same
# color-values 
distr_cf%>% 
  leaflet() %>% 
  addTiles()%>%
  #addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  #addProviderTiles(providers$Stamen.TonerBackground) %>% 
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>% 
  addPolygons(fillColor = ~cb(distr_cf$total_cf_cpo_tax200_kha),
              label = ~total_cf_cpo_tax100_kha,
              stroke = FALSE,
              weight = 1,
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 1) %>% 
  addLegend(pal = cb,  
            values = ~total_cf_cpo_tax200_kha, 
            bins = 4, opacity = 0.5,
            labFormat = labelFormat(suffix = " kha"),
            title = "Predicted deforestation <br/> 2001-2015, <br/> tax = $100/ton CPO",
            position = "bottomright")


# do not do a new colorNumeric palette, we want all maps to be plotting along the same
# color-values 
distr_cf%>% 
  leaflet() %>% 
  addTiles()%>%
  #addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  #addProviderTiles(providers$Stamen.TonerBackground) %>% 
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>% 
  addPolygons(fillColor = ~cb(total_cf_ffb_tax50_kha),
              label = ~district,
              stroke = FALSE,
              weight = 1,
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 1) %>% 
  addLegend(pal = cb,  
            values = ~total_cf_ffb_tax50_kha, 
            bins = 3, opacity = 0.5,
            labFormat = labelFormat(suffix = " kha"),
            title = "Predicted deforestation <br/> 2001-2015, <br/> tax = $50/ton FFB",
            position = "bottomright")

### MAP DIFFERENCE BASELINE - COUTERFACTUAL
# if we are to plot the differences between baseline and tax scenarii
### MAP
# if we are to plot the doistricts' accumulated lucpfip 
summary(distr_cf$diff_cpo_tax200_kha)
max_val <- max(distr_cf$diff_cpo_tax50_kha)+1
# bins <- seq(from = 40, to = 180, by = 20)
cb <- colorNumeric("BrBG", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                   domain = 0:max_val,
                   na.color = "red", 
                   reverse = FALSE)
distr_cf%>% 
  leaflet() %>% 
  addTiles()%>%
  #addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  #addProviderTiles(providers$Stamen.TonerBackground) %>% 
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>% 
  addPolygons(fillColor = ~cb(distr_cf$diff_cpo_tax50_kha),
              label = ~diff_cpo_tax200_kha,
              #color = "#444444",
              stroke = FALSE,
              weight = 1,
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 1) %>% 
  addLegend(pal = cb,  
            values = distr_cf$diff_cpo_tax50_kha, 
            bins = 4, opacity = 0.5,
            labFormat = labelFormat(suffix = " kha"),
            title = "Avoided conversion of<br/> 
            primary forest to oil palm<br/> 
            plantations 2001-2015",
            position = "bottomright") %>% 
  addLegend(title = "Counterfactual tax = $50/ton crude palm oil",
            position = "topright", 
            colors = NULL)

district_sf[district_sf$d__2000 == "_new",] %>% 
  leaflet() %>% 
  addTiles()%>%
  #addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  #addProviderTiles(providers$Stamen.TonerBackground) %>% 
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>% 
  addPolygons()


distr_cf[distr_cf$district=="Kab. Kotawaringin Timur","total_baseline_avg_kha"]

distr_cf$district %>% unique()

# bind df to plot
distr_cf <- rbind(elmts[["Sumatra"]][["distr_cf"]],
                  elmts[["Kalimantan"]][["distr_cf"]])
plot(distr_cf[,"total_baseline_avg_kha"])
plot(st_geometry(island_sf_prj), add = T)




# then tmap ? 
nrow(distr_cf_sf)
plot(distr_cf_sf[,"cf_ffb_tax10_ha"])

plot(st_geometry(distr_cf_sf))

### Generate demand and invert demand functions ### 
## FFB DEMAND FUNCTIONS
ffb_demand_avg <- function(tax, fitted_value_avg, aggr_factor, coeff_ffb){
  fitted_value_avg*aggr_factor/exp(coeff_ffb*tax)
}
ffb_demand_atavg <- function(tax, fitted_value_atavg, aggr_factor, coeff_ffb){
  fitted_value_atavg*aggr_factor/exp(coeff_ffb*tax)
}
ffb_demand_atmed <- function(tax, fitted_value_atmed, aggr_factor, coeff_ffb){
  fitted_value_atmed*aggr_factor/exp(coeff_ffb*tax)
}

# with aggr_factor = 1, we have the marginal effect. 
ffb_demand_avg(tax = 0, 
               fitted_value_avg = elmts[["Sumatra"]]["fitted_value_avg"],
               aggr_factor = 1,
               coeff_ffb = elmts[["Sumatra"]]["coeff_ffb"])
ffb_demand_atavg(tax = 0, 
                 fitted_value_atavg = elmts[["Sumatra"]]["fitted_value_atavg"],
                 aggr_factor = 1,
                 coeff_ffb = elmts[["Sumatra"]]["coeff_ffb"])
ffb_demand_atmed(tax = 0, 
                 fitted_value_atmed = elmts[["Sumatra"]]["fitted_value_atmed"],
                 aggr_factor = 1,
                 coeff_ffb = elmts[["Sumatra"]]["coeff_ffb"])


ffb_inv_demand_avg <- function(CF, fitted_value_avg, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(fitted_value_avg*aggr_factor) - log(CF))
}

ffb_inv_demand_atavg <- function(CF, fitted_value_atavg, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(fitted_value_atavg*aggr_factor) - log(CF))
}

ffb_inv_demand_atmed <- function(CF, fitted_value_atmed, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(fitted_value_atmed*aggr_factor) - log(CF))
}



### CPO DEMAND FUNCTIONS 
cpo_demand_avg <- function(tax, fitted_value_avg, aggr_factor, coeff_cpo){
  fitted_value_avg*aggr_factor/exp(coeff_cpo*tax)
}
cpo_demand_atavg <- function(tax, fitted_value_atavg, aggr_factor, coeff_cpo){
  fitted_value_atavg*aggr_factor/exp(coeff_cpo*tax)
}
cpo_demand_atmed <- function(tax, fitted_value_atmed, aggr_factor, coeff_cpo){
  fitted_value_atmed*aggr_factor/exp(coeff_cpo*tax)
}

cpo_demand_avg(tax = 0, 
               fitted_value_avg = elmts[["Sumatra"]]["fitted_value_avg"],
               aggr_factor = 1,
               coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])
cpo_demand_atavg(tax = 0, 
                 fitted_value_atavg = elmts[["Sumatra"]]["fitted_value_atavg"],
                 aggr_factor = 1,
                 coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])
cpo_demand_atmed(tax = 0, 
                 fitted_value_atmed = elmts[["Sumatra"]]["fitted_value_atmed"],
                 aggr_factor = 1,
                 coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])



cpo_inv_demand_avg <- function(CF, fitted_value_avg, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(fitted_value_avg*aggr_factor) - log(CF))
}

cpo_inv_demand_atavg <- function(CF, fitted_value_atavg, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(fitted_value_atavg*aggr_factor) - log(CF))
}

cpo_inv_demand_atmed <- function(CF, fitted_value_atmed, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(fitted_value_atmed*aggr_factor) - log(CF))
}























##### DEMAND FOR LUCFP #####
catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
island <- "Sumatra"
outcome_variable <- "lucfip_pixelcount_30th"
outcome_variable <- "lucpfip_pixelcount_total"
x_pya <- 3
lag_or_not <- "_lag1"
fixed_effects <- "parcel_id + district_year"
controls <- c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
weights <- TRUE

### Prepare the data ### 

elmts <- rep(list(NA), 3)
names(elmts) <- c("Sumatra", "Kalimantan", "All")

island_list <- list("Sumatra", "Kalimantan", c("Sumatra", "Kalimantan", "Papua"))
for(i in 1:length(island_list)){
    island <- island_list[[i]]
    elmts[[i]] <- rep(NA, 7)
    names(elmts[[i]]) <- c("aggr_factor", "coeff_ffb", "coeff_cpo", 
                           "fitted_value_atavg",
                           "fitted_value_atmed",
                           "fitted_value_avg",
                           "avg_lucpfip")
   
    # DATA FOR REGRESSIONS
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
    
    # DATA FOR AGGREGATION
    # the data are always at 50km (the largest CR)
    d_aggr <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                  parcel_size/1000,"km_50CR.rds")))
    
    # subsampling
    d_aggr <- d_aggr[d_aggr$island %in% island,]
    
    if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
      d_aggr <- d_aggr[d_aggr$any_pfc2000_total,]
    }
    if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
       outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
       outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
      d_aggr <- d_aggr[d_aggr$any_fc2000_30th,]
    }
    
    
    ### AGGREGATION FACTOR 
    length(unique(d_aggr$parcel_id))*length(unique(d_aggr$year)) == nrow(d_aggr)
    
    years <- unique(d_aggr$year)
    n_parcels_within_cr <- c() 
    for(t in years){
    n_parcels_within_cr[[match(t, years)]] <- d_aggr[d_aggr$year == t & d_aggr$n_reachable_uml > 0,] %>% nrow()
    }
    # we store this number in the first element of the vector of elements of each island 
    elmts[[i]]["aggr_factor"] <- sum(n_parcels_within_cr)

    
    # this is the total number of cells within CR and island, not the the number of cells used for estimation. 
    # this is matter of external validity.
    # therefore we read in the dataframes of grid cells within particular catchment radii, 
    # we adjust for the same restrictions (island and positive baseline forest extent) 
    # and count the number of parcels
    
    # n_cells <- length(unique(d$parcel_id))*length(unique(d$year))
    #names(n_cells) <- "n_cells"
    # then, multiply by 15, the number of years, if we want to compare with observed accumulated amouts
    # or if we want to say about the whole period. 
    # But keep the cross-sectional length (n_cells) if we want to say something like "each year", or "annually". 
    
    
    ### REGRESSION
    
    # there is no dynamics (distinction btw SR and LR effects) here as we aim at 
    # estimating the effect of a tax that remains in time. 
    regressors <- c(paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not),
                    paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)) 
    
    
    fml <- as.formula(paste0(outcome_variable,
                             " ~ ",
                             paste0(regressors, collapse = "+"),
                             " + ",
                             paste0(paste0(controls,lag_or_not), collapse = "+"),
                             " | ",
                             fixed_effects))
    
    if(weights == TRUE){
      var_weights <- d$sample_coverage_lag1/100 
      results <- fixest::feglm(fml,
                               data = d, 
                               family = "quasipoisson", 
                               notes = TRUE, 
                               weights = var_weights)
    }else{
      results <- fixest::feglm(fml,
                               data = d, 
                               family = "quasipoisson", 
                               notes = TRUE)
    }
    
    ## save coefficients 
    elmts[[i]]["coeff_ffb"] <- results$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
    elmts[[i]]["coeff_cpo"] <- results$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
    # names(mult_effect_ffb) <- NULL
    # names(mult_effect_cpo) <- NULL
    
    
    ## FITTED VALUE AT THE AVERAGE 
    covariate_means <- c(sapply(regressors, FUN = function(var){mean(d_aggr[,var], na.rm = TRUE)}),
                         sapply(paste0(controls, lag_or_not), FUN = function(var){mean(d_aggr[,var], na.rm = TRUE)}))
    
    linear_predictors_atavg <- c()
    for(covar in 1:length(covariate_means)){
      linear_predictors_atavg[covar] <- covariate_means[covar]*results$coefficients[covar]
    }
    # the issue with this is that covariate means are for the whole sample and not just the obs. 
    # used in the regression. (but is it very different, since those not used in the regression
    # are not precisely because they are missing and hence don't count in the mean)
    fitted_value_atavg <- exp(sum(linear_predictors_atavg) + mean(results$sumFE))
    elmts[[i]]["fitted_value_atavg"] <- fitted_value_atavg*(27.8*27.6)/(1e7)
  
    ## FITTED VALUE AT THE MEDIAN 
    covariate_medians <- c(sapply(regressors, FUN = function(var){median(d_aggr[,var], na.rm = TRUE)}),
                         sapply(paste0(controls, lag_or_not), FUN = function(var){median(d_aggr[,var], na.rm = TRUE)}))
    
    linear_predictors_atmed <- c()
    for(covar in 1:length(covariate_medians)){
      linear_predictors_atmed[covar] <- covariate_medians[covar]*results$coefficients[covar]
    }
  
    fitted_value_atmed <- exp(sum(linear_predictors_atmed) + median(results$sumFE))
    elmts[[i]]["fitted_value_atmed"] <- fitted_value_atmed*(27.8*27.6)/(1e7)
    
    ## AVERAGE FITTED VALUE 
    fitted_value_avg <- mean(results$fitted.values) # this is equal to mean(exp(results$linear.predictors))
    # it is to say that linear.predictors already encompass the fixed effects.
    # It's in pixelcount, convert it to thousand hectares
    elmts[[i]]["fitted_value_avg"] <- fitted_value_avg*(27.8*27.6)/(1e7)
    
    ## AVERAGE ACTUAL OUTCOME
    ## if we use actual lucpfip and not fitted values - this is just for informative purpose. 
    # If the sample over which the mean is computed is restricted to n_reachable_uml > 0, 
    # avg_lucpfip*aggr_factor gives the accumulated LUCPFIP in table 2. 
    # If the sample is not restricted, avg_lucpfip*nrow(d) gives the accumulated LUCFP area computed initially
    # in table 2 (over all parcels that are within a mill's CR at least one year)
    
    avg_lucpfip <- mean(d_aggr[d_aggr$n_reachable_uml > 0,outcome_variable]) # d$n_reachable_uml > 0
    elmts[[i]]["avg_lucpfip"] <- avg_lucpfip*(27.8*27.6)/(1e7)
    
    ## MARGINAL EFFECTS 
    # for a premium / tax of USD10/ton FFB (~1/3 of a std.dev) (and rescale to hectares by *1e3)
    # ffb_marginal_effect_avg <- coeff_ffb*fitted_value_avg*1e3*10 # average partial effect
    # ffb_marginal_effect_atavg <- coeff_ffb*fitted_value_atavg*1e3*10 # partial effect at average 
    # ffb_marginal_effect_atmed <- coeff_ffb*fitted_value_atmed*1e3*10 # partial effect at average 
    # ffb_marginal_effect_actual <- coeff_ffb*avg_lucpfip*1e3*10
    # cpo_marginal_effect_avg <- coeff_cpo*fitted_value_avg*1e3*10 # average partial effect
    # cpo_marginal_effect_atavg <- coeff_cpo*fitted_value_atavg*1e3*10 # partial effect at average 
    # cpo_marginal_effect_atmed <- coeff_cpo*fitted_value_atmed*1e3*10 # partial effect at average 
    # cpo_marginal_effect_actual <- coeff_cpo*avg_lucpfip*1e3*10
    # 
    # ffb_marginal_effect_avg 
    # ffb_marginal_effect_atavg 
    # ffb_marginal_effect_atmed 
    # ffb_marginal_effect_actual 
    # cpo_marginal_effect_avg 
    # cpo_marginal_effect_atavg 
    # cpo_marginal_effect_atmed 
    # cpo_marginal_effect_actual 
    
  
    ######### test manual computation of partial effects at average with mfx package ###################
    # library(mfx)
    # fml2 <- lucpfip_pixelcount_total ~ wa_ffb_price_imp1_4ya_lag1 + wa_cpo_price_imp1_4ya_lag1 + 
    #   wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + 
    #   wa_pct_own_for_imp_lag1 + n_reachable_uml_lag1
    # 
    # mfx <- poissonmfx(formula = fml2, 
    #            data = d,
    #            atmean = TRUE)
    # 
    # results <- fixest::feglm(fml2,
    #                          data = d, 
    #                          family = "poisson", 
    #                          notes = TRUE)
    # 
    # # now reproduce exactly the same, manually: 
    # 
    # 
    # coeff_ffb_mfx <- mfx$fit$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
    # coeff_cpo_mfx <- mfx$fit$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
    # coeff_ffb_fixest <- results$coefficients[paste0("wa_ffb_price_imp1_",x_pya+1,"ya",lag_or_not)] 
    # coeff_cpo_fixest <- results$coefficients[paste0("wa_cpo_price_imp1_",x_pya+1,"ya",lag_or_not)]
    # 
    # # they are not exactly equal 
    # coeff_ffb_mfx == coeff_ffb_fixest
    # round(coeff_ffb_mfx, 6) == round(coeff_ffb_fixest,6)
    # 
    # ## at average, manually
    # covariate_means <- c(sapply(regressors, FUN = function(var){mean(d[,var], na.rm = TRUE)}),
    #                      sapply(paste0(controls, lag_or_not), FUN = function(var){mean(d[,var], na.rm = TRUE)}))
    # 
    # # based on mfx (glm)
    # linear_predictors_mfx <- c()
    # for(i in 1:length(covariate_means)){
    #   linear_predictors_mfx[i] <- covariate_means[i]*mfx$fit$coefficients[i+1] # +1 bc of the intercept
    # }
    # fitted_value_atavg_mfx <- exp(sum(linear_predictors_mfx) + mfx$fit$coefficients[1]) # linear predictors don't integrate the intercept.
    # ffb_marginal_effect_atavg_mfx <- coeff_ffb_mfx*fitted_value_atavg_mfx
    # cpo_marginal_effect_atavg_mfx <- coeff_cpo_mfx*fitted_value_atavg_mfx
    # 
    # # based on fixest coeffs
    # linear_predictors_fixest <- c()
    # for(i in 1:length(covariate_means)){
    #   linear_predictors_fixest[i] <- covariate_means[i]*results$coefficients[i+1]
    # }
    # 
    # fitted_value_atavg_fixest <- exp(sum(linear_predictors_fixest) + results$coefficients[1])
    # ffb_marginal_effect_atavg_fixest <- coeff_ffb_fixest*fitted_value_atavg_fixest
    # cpo_marginal_effect_atavg_fixest <- coeff_cpo_fixest*fitted_value_atavg_fixest
    # 
    # # because coeffs are not exactly equal, linear predictors are not either, nor fitted values and marginal effects
    # linear_predictors_mfx == linear_predictors_fixest
    # round(linear_predictors_mfx, 6) == round(linear_predictors_fixest, 6)
    # fitted_value_atavg_mfx == fitted_value_atavg_fixest
    # ffb_marginal_effect_atavg_mfx == ffb_marginal_effect_atavg_fixest
    # round(ffb_marginal_effect_atavg_mfx, 6) == round(ffb_marginal_effect_atavg_fixest, 6)
    # 
    # # but this does not explain the gap with the marginal effect directly computed by mfx
    # mfx$mfxest[1,1] == ffb_marginal_effect_atavg_mfx
    # round(mfx$mfxest[1,1], 6) == round(ffb_marginal_effect_atavg_mfx,6) # is FALSE
    # 
    # mfx <- poissonmfx(formula = fml2, 
    #            data = d,
    #            atmean = FALSE)
    # 
    # # recover the fitted value at the average from mfx marginal effect
    # # bc we know the marginal effect at average is the product of the coeff and the fitted value at the average
    # mfx$mfxest[1,1]/coeff_ffb_mfx == fitted_value_atavg_mfx
    # # is false, therefore it's really the fitted value at the average that is not computed exactly the same in 
    # # mfx and manually 
    # 
    # ## AVERAGE FITTED VALUE 
    # fitted_value_avg <- mean(results$fitted.values)
    # ffb_marginal_effect_avg <- coeff_ffb*fitted_value_avg # average partial effect
    # cpo_marginal_effect_avg <- coeff_cpo*fitted_value_avg # average partial effect
    # 
    # poissonmfx(formula = fml2, 
    #                   data = d,
    #                   atmean = FALSE)
    # 
    # # so for average partial effect, results are equal. 
    # # NOT for partial effect at the average... NOW WITH THE INTERCEPT ADDED IT'S VERY SIMILAR, BUT STILL NOT EQUAL !!!
    # # --> try to compute manually on a balanced panel (remove records that have at least one rhs missing)
    # # or build a basic balanced df. 
    # # it's not coeff that differ since for average partial effect it works. 
    # # voir ce qu'est glm$effects (mfx$fit$effects)
    # 
    # 
    # ### and now with {margins}
    # library(margins)
    #### end of test ####
}
    
  


### Generate demand and invert demand functions ### 
## FFB DEMAND FUNCTIONS
ffb_demand_avg <- function(tax, fitted_value_avg, aggr_factor, coeff_ffb){
  fitted_value_avg*aggr_factor/exp(coeff_ffb*tax)
}
ffb_demand_atavg <- function(tax, fitted_value_atavg, aggr_factor, coeff_ffb){
  fitted_value_atavg*aggr_factor/exp(coeff_ffb*tax)
}
ffb_demand_atmed <- function(tax, fitted_value_atmed, aggr_factor, coeff_ffb){
  fitted_value_atmed*aggr_factor/exp(coeff_ffb*tax)
}
ffb_demand_actual <- function(tax, avg_lucpfip, aggr_factor, coeff_ffb){
  avg_lucpfip*aggr_factor/exp(coeff_ffb*tax)
}

ffb_demand_avg(tax = 0, 
               fitted_value_avg = elmts[["Kalimantan"]]["fitted_value_avg"],
               aggr_factor = elmts[["Kalimantan"]]["aggr_factor"],
               coeff_ffb = elmts[["Kalimantan"]]["coeff_ffb"])
ffb_demand_atavg(tax = 0, 
               fitted_value_atavg = elmts[["Kalimantan"]]["fitted_value_atavg"],
               aggr_factor = elmts[["Kalimantan"]]["aggr_factor"],
               coeff_ffb = elmts[["Kalimantan"]]["coeff_ffb"])
ffb_demand_atmed(tax = 0, 
                 fitted_value_atmed = elmts[["Kalimantan"]]["fitted_value_atmed"],
                 aggr_factor = elmts[["Kalimantan"]]["aggr_factor"],
                 coeff_ffb = elmts[["Kalimantan"]]["coeff_ffb"])
ffb_demand_actual(tax = 0, 
                  avg_lucpfip = elmts[["Kalimantan"]]["avg_lucpfip"],
                  aggr_factor = elmts[["Kalimantan"]]["aggr_factor"],
                  coeff_ffb = elmts[["Kalimantan"]]["coeff_ffb"]) 

#avg_lucpfip*nrow(d) 
# we get exactly the same amount as the actual accumulated LUCPFIP in table 2 (in its first 
# version where accumulation was over all parcels that were within a mill's CR at least one year. 

# Invert demands  
ffb_inv_demand_avg <- function(CF, fitted_value_avg, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(fitted_value_avg*aggr_factor) - log(CF))
}

ffb_inv_demand_atavg <- function(CF, fitted_value_atavg, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(fitted_value_atavg*aggr_factor) - log(CF))
}

ffb_inv_demand_atmed <- function(CF, fitted_value_atmed, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(fitted_value_atmed*aggr_factor) - log(CF))
}

ffb_inv_demand_actual <- function(CF, avg_lucpfip, aggr_factor, coeff_ffb){
  (1/coeff_ffb)*(log(avg_lucpfip*aggr_factor) - log(CF))
}


#### CPO DEMAND FUNCTIONS ####
cpo_demand_avg <- function(tax, fitted_value_avg, aggr_factor, coeff_cpo){
  fitted_value_avg*aggr_factor/exp(coeff_cpo*tax)
}
cpo_demand_atavg <- function(tax, fitted_value_atavg, aggr_factor, coeff_cpo){
  fitted_value_atavg*aggr_factor/exp(coeff_cpo*tax)
}
cpo_demand_atmed <- function(tax, fitted_value_atmed, aggr_factor, coeff_cpo){
  fitted_value_atmed*aggr_factor/exp(coeff_cpo*tax)
}
cpo_demand_actual <- function(tax, avg_lucpfip, aggr_factor, coeff_cpo){
  avg_lucpfip*aggr_factor/exp(coeff_cpo*tax)
}

cpo_demand_avg(tax = 0, 
               fitted_value_avg = elmts[["Sumatra"]]["fitted_value_avg"],
               aggr_factor = elmts[["Sumatra"]]["aggr_factor"],
               coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])
cpo_demand_atavg(tax = 0, 
                 fitted_value_atavg = elmts[["Sumatra"]]["fitted_value_atavg"],
                 aggr_factor = elmts[["Sumatra"]]["aggr_factor"],
                 coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])
cpo_demand_atmed(tax = 0, 
                 fitted_value_atmed = elmts[["Sumatra"]]["fitted_value_atmed"],
                 aggr_factor = elmts[["Sumatra"]]["aggr_factor"],
                 coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])
cpo_demand_actual(tax = 0, 
                  avg_lucpfip = elmts[["Sumatra"]]["avg_lucpfip"],
                  aggr_factor = elmts[["Sumatra"]]["aggr_factor"],
                  coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"]) 

#avg_lucpfip*nrow(d) 
# we get exactly the same amount as the actual accumulated LUCPFIP in table 2 (in its first 
# version where accumulation was over all parcels that were within a mill's CR at least one year. 

# Invert demands  
cpo_inv_demand_avg <- function(CF, fitted_value_avg, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(fitted_value_avg*aggr_factor) - log(CF))
}

cpo_inv_demand_atavg <- function(CF, fitted_value_atavg, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(fitted_value_atavg*aggr_factor) - log(CF))
}

cpo_inv_demand_atmed <- function(CF, fitted_value_atmed, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(fitted_value_atmed*aggr_factor) - log(CF))
}

cpo_inv_demand_actual <- function(CF, avg_lucpfip, aggr_factor, coeff_cpo){
  (1/coeff_cpo)*(log(avg_lucpfip*aggr_factor) - log(CF))
}


### PLOT DEMAND FOR DEFORESTATION
library(ggplot2)

ffb_axis_right_bound <- ffb_demand_atavg(tax = 0, 
                                         fitted_value_atavg = elmts[["All"]]["fitted_value_atavg"],
                                         aggr_factor = elmts[["All"]]["aggr_factor"],
                                         coeff_ffb = elmts[["All"]]["coeff_ffb"])

ggplot(data = data.frame(CF=c(0, ffb_axis_right_bound)), 
       aes(x=CF)) + 
  stat_function(fun=ffb_inv_demand_atavg, 
                args = list(fitted_value_atavg = elmts[["All"]]["fitted_value_atavg"],
                                         aggr_factor = elmts[["All"]]["aggr_factor"],
                                         coeff_ffb = elmts[["All"]]["coeff_ffb"])) + 
  stat_function(fun=ffb_inv_demand_atavg, 
                args = list(fitted_value_atavg = elmts[["Sumatra"]]["fitted_value_atavg"],
                            aggr_factor = elmts[["Sumatra"]]["aggr_factor"],
                            coeff_ffb = elmts[["Sumatra"]]["coeff_ffb"]), 
                colour = "blue") +
  stat_function(fun=ffb_inv_demand_atavg, 
                args = list(fitted_value_atavg = elmts[["Kalimantan"]]["fitted_value_atavg"],
                            aggr_factor = elmts[["Kalimantan"]]["aggr_factor"],
                            coeff_ffb = elmts[["Kalimantan"]]["coeff_ffb"]),
                colour = "red") + 
  labs(x = "Land use change from primary forest to industrial plantations (kha)", y = "tax on LUCFP (USD/ton FFB)", 
       title = "Demand for deforestation, 2001 - 2015")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Bon on s'arrête là pour l'instant. Le problème avec le fait de plotter toutes les îles 
# ensemble c'est qu'elles ne sont pas sur la même scale. On pourrait faire en sorte que la partie 
# tax < 0 disparaisse pour Kalimantan, mais bon même là c'est pas hyper intéressant de les 
# avoir toutes les 3 sur le même graph. Peut être qu'il sera plus pertinent in fine d'avoir 
# des graphs par île avec des courbes pour les différents groupes (industrial, medium, small plantations)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# plot provisoire pour mettre dans le latex au 17/08/2020, juste pour illustrer la partie 
# demand for deforestation. Mais ce n'est clairement pas le plot définitif

ffb_axis_right_bound <- ffb_demand_atmed(tax = 0, 
                                         fitted_value_atmed = elmts[["All"]]["fitted_value_atmed"],
                                         aggr_factor = elmts[["All"]]["aggr_factor"],
                                         coeff_ffb = elmts[["All"]]["coeff_ffb"])
ggplot(data = data.frame(CF=c(0, ffb_axis_right_bound)), 
       aes(x=CF)) + 
  stat_function(fun=ffb_inv_demand_atmed, 
                args = list(fitted_value_atmed = elmts[["All"]]["fitted_value_atmed"],
                            aggr_factor = elmts[["All"]]["aggr_factor"],
                            coeff_ffb = elmts[["All"]]["coeff_ffb"])) + 
  labs(x = "Land use change from primary forest to industrial plantations (kha)", y = "tax on LUCFP (USD/ton FFB)")

# effect of a premium of 25$/ton FFB
ffb_demand_atavg(0) - ffb_demand_atavg(25)
# effect of a premium of 100$/ton CPO
cpo_demand_atavg(0) - cpo_demand_atavg(100)



 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##### TESTING ZONE ##### 

## test whether the FE are in linear.predictors 
# linear.predictors, sumFE and fitted.values all have the same length: the number of obs. used in the reg
length(results$linear.predictors)==length(results$sumFE)
length(results$sumFE)==length(results$fitted.values)

# the fitted.values are equal to the exponential of the linear.predictors (here)
results$linear.predictors[1] 
results$sumFE[1]
results$fitted.values[4] == exp(results$linear.predictors[4])
results$linear.predictors[1] + results$sumFE[1]

fml_nofe <- as.formula(paste0(outcome_variable,
                              " ~ ",
                              paste0(regressors, collapse = "+"),
                              " + ",
                              paste0(paste0(controls,lag_or_not), collapse = "+")))

results_nofe <- fixest::feglm(fml_nofe,
                              data = d, 
                              family = "quasipoisson", 
                              notes = TRUE)
length(results_nofe$linear.predictors) - length(results$linear.predictors)

any(results$linear.predictors %in% results_nofe$linear.predictors)

any(results_nofe$linear.predictors %in% results$linear.predictors)

results_nofe$fitted.values[4] == exp(results_nofe$linear.predictors[4])

# so linear predictors are different when there are FE and where there is no. 
# which means that linear predictors INCLUDE FE. 

if(commo == "ffb" & dynamics == "both-yoyg" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- as.symbol(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "ffb" & dynamics == "both-3pya" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-yoyg" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1  +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-3pya" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year 
}
if(commo == "both" & dynamics == "both-yoyg" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
      wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_yoyg_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_yoyg_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "both" & dynamics == "both-3pya" & lagged == TRUE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
      wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
      wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
      n_reachable_uml_lag1 | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya_lag1 + wa_ffb_price_imp1_3pya_lag1 +
    wa_cpo_price_imp1_dev_3pya_lag1 + wa_cpo_price_imp1_3pya_lag1 +
    wa_pct_own_loc_gov_imp_lag1 + wa_pct_own_nat_priv_imp_lag1 + wa_pct_own_for_imp_lag1 +
    n_reachable_uml_lag1 | 
    parcel_id + district^year
}
if(commo == "ffb" & dynamics == "both-yoyg" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "ffb" & dynamics == "both-3pya" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-yoyg" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya  +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "cpo" & dynamics == "both-3pya" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year 
}
if(commo == "both" & dynamics == "both-yoyg" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
      wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_yoyg_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_yoyg_3pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}
if(commo == "both" & dynamics == "both-3pya" & lagged == FALSE){
  
  # unit FE
  spec_lucfip_ufe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_3pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id
  
  # time FE
  spec_lucfip_tfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    year
  
  # unit FE + year FE i.e. TWFE
  spec_lucfip_twfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + year
  
  # unit FE + island*year FE
  if(length(island) > 1){
    spec_lucfip_iyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
      wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
      wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
      n_reachable_uml | 
      parcel_id + island^year
  }
  
  # unit FE + province*year FE
  spec_lucfip_pyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + province^year
  
  # unit FE + district*year FE
  spec_lucfip_dyfe <- get(outcome_variable) ~ wa_ffb_price_imp1_dev_3pya + wa_ffb_price_imp1_4pya +
    wa_cpo_price_imp1_dev_3pya + wa_cpo_price_imp1_4pya +
    wa_pct_own_loc_gov_imp + wa_pct_own_nat_priv_imp + wa_pct_own_for_imp +
    n_reachable_uml | 
    parcel_id + district^year
}


### Regressions
qpoiglm_est_lucfip_ufe <- fixest::feglm(spec_lucfip_ufe, data = d, family = "quasipoisson", notes = FALSE)
qpoiglm_est_lucfip_tfe <- fixest::feglm(spec_lucfip_tfe, data = d, family = "quasipoisson", notes = FALSE)
qpoiglm_est_lucfip_twfe <- fixest::feglm(spec_lucfip_twfe, data = d, family = "quasipoisson", notes = FALSE)
if(length(island) > 1){
  qpoiglmspec_lucfip_iyfe <- fixest::feglm(spec_lucfip_iyfe, data = d, family = "quasipoisson", notes = FALSE)
}
qpoiglmspec_lucfip_dyfe <- fixest::feglm(spec_lucfip_dyfe, data = d, family = "quasipoisson", notes = FALSE)
qpoiglmspec_lucfip_pyfe <- fixest::feglm(spec_lucfip_pyfe, data = d, family = "quasipoisson", notes = FALSE)

# # give it a panel class for lag analyses
# pd <- panel(d, panel.id = c("parcel_id", "year"))

# on the sample of parcels covered with a positive extent of forest of 30% tree canopy density in 2000.
# using the area in hectare of LUCFIP on forest of 30% tree canopy density (outside industrial plantations in 2000)
est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1), data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+district, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+district+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island^year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_ffb_price_imp1)|parcel_id+island+district^year, data = df30_01)
etable(est_f30)



### CPO
# on the sample of parcels covered with a positive extent of forest of 30% tree canopy density in 2000.
# using the area in hectare of LUCFIP on forest of 30% tree canopy density (outside industrial plantations in 2000)
est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1), data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+district, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+district+year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island^year, data = df30_01)
etable(est_f30)

est_f30 <- feglm(lucfip_ha_30th~log(wa_cpo_price_imp1)|parcel_id+island+district^year, data = df30_01)
etable(est_f30)







# PREMIER SET 126 -31
# 2EME SET 41.85 - 12
# 3EME SET 219.15 - 66
36+123.2+39.6+15.4+4.95




