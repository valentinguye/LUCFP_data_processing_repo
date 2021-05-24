# regressions_CR.R from 31012021 commit 

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

### Set dictionary for names of variables to display in regression tables 
setFixest_dict(c(parcel_id = "grid cell",
                 lucpfip_pixelcount_total = "Land use change from primary forest to industrial oil palm plantations", 
                 lucpfip_pixelcount = "Industrial plantations", 
                 lucpfip_rapid_pixelcount = "Rapid conversion to industrial plantations",
                 lucpfip_slow_pixelcount = "Slow conversion to industrial plantations",
                 lucpfsmp_pixelcount = "Smallholder plantations", 
                 lucfip_pixelcount = "Land use change from 30% tree cover forest to industrial oil palm plantations", 
                 lucfsmp_pixelcount = "Land use change from 30% tree cover forest to small or medium-sized oil palm plantations", 
                 lucpfap_pixelcount = "Overall", 
                 lucfap_pixelcount = "Overall", 
                 
                 lucfip_ha_30th = "LUCFIP (30 pct. canopy density, ha)",
                 lucfip_ha_60th = "LUCFIP (60 pct. canopy density, ha)",
                 lucfip_ha_90th = "LUCFIP (90 pct. canopy density, ha)",
                 lucfip_pixelcount_60th = "LUCFIP (60 pct. canopy density, pixels)",
                 lucfip_pixelcount_90th = "LUCFIP (90 pct. canopy density, pixels)",
                 # No time dynamics FFB variables 
                 wa_ffb_price_imp1_3ya = "FFB price, 3-year average",
                 wa_ffb_price_imp1_3ya_lag1 = "FFB price, 3-year average",#(lagged)
                 wa_ffb_price_imp1_4ya = "FFB price, 4-year average",
                 wa_ffb_price_imp1_4ya_lag1 = "FFB price, 4-year average",#(lagged)
                 wa_ffb_price_imp1_5ya = "FFB price, 5-year average",
                 wa_ffb_price_imp1_5ya_lag1 = "FFB price, 5-year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_3ya = "FFB price y-o-y growth rate, 3-year average",
                 wa_ffb_price_imp1_yoyg_3ya_lag1 = "FFB price y-o-y growth rate, 3-year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_4ya = "FFB price y-o-y growth rate, 4-year average",
                 wa_ffb_price_imp1_yoyg_4ya_lag1 = "FFB price y-o-y growth rate, 4-year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_5ya = "FFB price y-o-y growth rate, 5-year average",
                 wa_ffb_price_imp1_yoyg_5ya_lag1 = "FFB price y-o-y growth rate, 5-year average",#(lagged)
                 ln_wa_ffb_price_imp1_3ya = "FFB price, 3-year average",
                 ln_wa_ffb_price_imp1_3ya_lag1 = "FFB price, 3-year average",#(lagged)
                 ln_wa_ffb_price_imp1_4ya = "FFB price, 4-year average",
                 ln_wa_ffb_price_imp1_4ya_lag1 = "FFB price, 4-year average",#(lagged)
                 ln_wa_ffb_price_imp1_5ya = "FFB price, 5-year average",
                 ln_wa_ffb_price_imp1_5ya_lag1 = "FFB price, 5-year average",#(lagged)
                 # No time dynamics CPO variables 
                 wa_cpo_price_imp1_3ya = "CPO price, 3-year average",
                 wa_cpo_price_imp1_3ya_lag1 = "CPO price, 3-year average",#(lagged)
                 wa_cpo_price_imp1_4ya = "CPO price, 4-year average",
                 wa_cpo_price_imp1_4ya_lag1 = "CPO price, 4-year average",#(lagged)
                 wa_cpo_price_imp1_5ya = "CPO price, 5-year average",
                 wa_cpo_price_imp1_5ya_lag1 = "CPO price, 5-year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_3ya = "CPO price y-o-y growth rate, 3-year average",
                 wa_cpo_price_imp1_yoyg_3ya_lag1 = "CPO price y-o-y growth rate, 3-year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_4ya = "CPO price y-o-y growth rate, 4-year average",
                 wa_cpo_price_imp1_yoyg_4ya_lag1 = "CPO price y-o-y growth rate, 4-year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_5ya = "CPO price y-o-y growth rate, 5-year average",
                 wa_cpo_price_imp1_yoyg_5ya_lag1 = "CPO price y-o-y growth rate, 5-year average",#(lagged)
                 ln_wa_cpo_price_imp1_3ya = "CPO price, 3-year average",
                 ln_wa_cpo_price_imp1_3ya_lag1 = "CPO price, 3-year average",#(lagged)
                 ln_wa_cpo_price_imp1_4ya = "CPO price, 4-year average",
                 ln_wa_cpo_price_imp1_4ya_lag1 = "CPO price, 4-year average",#(lagged)
                 ln_wa_cpo_price_imp1_5ya = "CPO price, 5-year average",
                 ln_wa_cpo_price_imp1_5ya_lag1 = "CPO price, 5-year average",#(lagged)
                 # SR FFB variables
                 wa_ffb_price_imp1 = "FFB price",
                 wa_ffb_price_imp1_lag1 = "FFB price",#(lagged)
                 wa_ffb_price_imp1_yoyg = "FFB price y-o-y growth rate",
                 wa_ffb_price_imp1_yoyg_lag1 = "FFB price y-o-y growth rate",#(lagged)
                 wa_ffb_price_imp1_dev_2pya = "FFB price deviation from 2 past year average",
                 wa_ffb_price_imp1_dev_2pya_lag1 = "FFB price deviation from 2 past year average",#(lagged)
                 wa_ffb_price_imp1_dev_3pya = "FFB price deviation from 3 past year average",
                 wa_ffb_price_imp1_dev_3pya_lag1 = "FFB price deviation from 3 past year average",#(lagged)
                 wa_ffb_price_imp1_dev_4pya = "FFB price deviation from 4 past year average",
                 wa_ffb_price_imp1_dev_4pya_lag1 = "FFB price deviation from 4 past year average",#(lagged)
                 ln_wa_ffb_price_imp1 = "FFB price",
                 ln_wa_ffb_price_imp1_lag1 = "FFB price",#(lagged)
                 # LR FFB variables
                 wa_ffb_price_imp1_2pya = "FFB price, 2 past year average",
                 wa_ffb_price_imp1_2pya_lag1 = "FFB price, 2 past year average",#(lagged)
                 wa_ffb_price_imp1_3pya = "FFB price, 3 past year average",
                 wa_ffb_price_imp1_3pya_lag1 = "FFB price, 3 past year average",#(lagged)
                 wa_ffb_price_imp1_4pya = "FFB price, 4 past year average",
                 wa_ffb_price_imp1_4pya_lag1 = "FFB price, 4 past year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_2pya = "FFB price y-o-y growth rate, 2 past year average",
                 wa_ffb_price_imp1_yoyg_2pya_lag1 = "FFB price y-o-y growth rate, 2 past year average",#(lagged)                
                 wa_ffb_price_imp1_yoyg_3pya = "FFB price y-o-y growth rate, 3 past year average",
                 wa_ffb_price_imp1_yoyg_3pya_lag1 = "FFB price y-o-y growth rate, 3 past year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_4pya = "FFB price y-o-y growth rate, 4 past year average",
                 wa_ffb_price_imp1_yoyg_4pya_lag1 = "FFB price y-o-y growth rate, 4 past year average",#(lagged)
                 ln_wa_ffb_price_imp1_2pya = "FFB price, 2 past year average",
                 ln_wa_ffb_price_imp1_2pya_lag1 = "FFB price, 2 past year average",#(lagged)
                 ln_wa_ffb_price_imp1_3pya = "FFB price, 3 past year average",
                 ln_wa_ffb_price_imp1_3pya_lag1 = "FFB price, 3 past year average",#(lagged)
                 ln_wa_ffb_price_imp1_4pya = "FFB price, 4 past year average",
                 ln_wa_ffb_price_imp1_4pya_lag1 = "FFB price, 4 past year average",#(lagged)                 
                 # SR CPO variables
                 wa_cpo_price_imp1 = "CPO price",
                 wa_cpo_price_imp1_lag1 = "CPO price",#(lagged)
                 wa_cpo_price_imp1_yoyg = "CPO price y-o-y growth rate",
                 wa_cpo_price_imp1_yoyg_lag1 = "CPO price y-o-y growth rate",#(lagged)
                 wa_cpo_price_imp1_dev_2pya = "CPO price deviation from 2 past year average",
                 wa_cpo_price_imp1_dev_2pya_lag1 = "CPO price deviation from 2 past year average",#(lagged) 
                 wa_cpo_price_imp1_dev_3pya = "CPO price deviation from 3 past year average",
                 wa_cpo_price_imp1_dev_3pya_lag1 = "CPO price deviation from 3 past year average",#(lagged) 
                 wa_cpo_price_imp1_dev_4pya = "CPO price deviation from 4 past year average",
                 wa_cpo_price_imp1_dev_4pya_lag1 = "CPO price deviation from 4 past year average",#(lagged) 
                 ln_wa_cpo_price_imp1 = "CPO price",
                 ln_wa_cpo_price_imp1_lag1 = "CPO price",#(lagged)
                 # LR CPO variables
                 wa_cpo_price_imp1_2pya = "CPO price, 2 past year average",
                 wa_cpo_price_imp1_2pya_lag1 = "CPO price, 2 past year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_2pya = "CPO price y-o-y growth rate, 2 past year average",
                 wa_cpo_price_imp1_yoyg_2pya_lag1 = "CPO price y-o-y growth rate, 2 past year average",#(lagged)
                 wa_cpo_price_imp1_3pya = "CPO price, 3 past year average",
                 wa_cpo_price_imp1_3pya_lag1 = "CPO price, 3 past year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_3pya = "CPO price y-o-y growth rate, 3 past year average",
                 wa_cpo_price_imp1_yoyg_3pya_lag1 = "CPO price y-o-y growth rate, 3 past year average",#(lagged)
                 wa_cpo_price_imp1_4pya = "CPO price, 4 past year average",
                 wa_cpo_price_imp1_4pya_lag1 = "CPO price, 4 past year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_4pya = "CPO price y-o-y growth rate, 4 past year average",
                 wa_cpo_price_imp1_yoyg_4pya_lag1 = "CPO price y-o-y growth rate, 4 past year average",#(lagged)
                 ln_wa_cpo_price_imp1_2pya = "CPO price, 2 past year average",
                 ln_wa_cpo_price_imp1_2pya_lag1 = "CPO price, 2 past year average",#(lagged)
                 ln_wa_cpo_price_imp1_3pya = "CPO price, 3 past year average",
                 ln_wa_cpo_price_imp1_3pya_lag1 = "CPO price, 3 past year average",#(lagged)
                 ln_wa_cpo_price_imp1_4pya = "CPO price, 4 past year average",
                 ln_wa_cpo_price_imp1_4pya_lag1 = "CPO price, 4 past year average",#(lagged)
                 ## interactions 
                 n_reachable_uml_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "# reachable mills X FFB price",
                 n_reachable_uml_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "# reachable mills X CPO price",
                 wa_pct_own_loc_gov_imp_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "local gvt mill ownership X FFB price",
                 wa_pct_own_loc_gov_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "local gvt mill ownership X CPO price",
                 wa_pct_own_nat_priv_imp_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "private mill ownership X FFB price",
                 wa_pct_own_nat_priv_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "private mill ownership X CPO price",
                 wa_pct_own_for_imp_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "foreign mill ownership X FFB price",
                 wa_pct_own_for_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "foreign mill ownership X CPO price",
                 wa_prex_cpo_imp1_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "pct. CPO exported X FFB price",
                 wa_prex_cpo_imp1_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "pct. CPO exported X CPO price",
                 ## controls
                 lucpfip_pixelcount_total_lag1 = "LUCPFIP (pixels, lagged)",
                 lucpfip_pixelcount_total_2pya = "LUCPFIP (pixels, 2 past year average)",
                 lucpfip_pixelcount_total_3pya = "LUCPFIP (pixels, 3 past year average)",
                 lucpfip_pixelcount_total_4pya = "LUCPFIP (pixels, 4 past year average)",
                 lucpfsmp_pixelcount_total_lag1 = "LUCPFSMP (pixels, lagged)",
                 lucpfsmp_pixelcount_total_2pya = "LUCPFSMP (pixels, 2 past year average)",
                 lucpfsmp_pixelcount_total_3pya = "LUCPFSMP (pixels, 3 past year average)",
                 lucpfsmp_pixelcount_total_4pya = "LUCPFSMP (pixels, 4 past year average)",
                 n_reachable_uml = "# reachable UML mills",
                 n_reachable_uml_lag1 = "# reachable UML mills",#(lagged)
                 wa_pct_own_cent_gov_imp = "Local government mill ownership",
                 wa_pct_own_cent_gov_imp_lag1 = "Local government mill ownership",
                 wa_pct_own_loc_gov_imp = "Local government mill ownership",
                 wa_pct_own_loc_gov_imp_lag1 = "Local government mill ownership",
                 wa_pct_own_nat_priv_imp = "Domestic private mill ownership",
                 wa_pct_own_nat_priv_imp_lag1 = "Domestic private mill ownership",
                 wa_pct_own_for_imp = "Foreign mill ownership",
                 wa_pct_own_for_imp_lag1 = "Foreign mill ownership", 
                 wa_prex_cpo_imp1 = "Percentage CPO exported",
                 wa_prex_cpo_imp1_lag1 = "Percentage CPO exported",#(lagged)
                 wa_prex_cpo_imp2 = "Percentage CPO exported",
                 wa_prex_cpo_imp2_lag1 = "Percentage CPO exported"#(lagged)
))



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
# catchment = "CR"
# outcome_variable = "lucpfap_pixelcount"
# island = "both"
# alt_cr = FALSE
# nearest_mill = FALSE
# margin = "both"
# restr_marg_def = TRUE
# commo = c("cpo")
# x_pya = 3
# dynamics = FALSE
# log_prices = TRUE
# yoyg = FALSE
# only_sr = FALSE
# short_run = "full"
# imp = 1
# distribution = "quasipoisson"
# fe = "parcel_id + district_year"#
# offset = FALSE
# lag_or_not = "_lag1"
# controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml")#, "wa_prex_cpo_imp1""wa_pct_own_loc_gov_imp",
# remaining_forest = FALSE
# interaction_terms = NULL # "illegal2"  #c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1")
# interact_regressors = TRUE
# interacted = "regressors"
# spatial_control = FALSE
# pya_ov = FALSE
# illegal = "all"# "ill2" #
# weights = FALSE
# min_forest_2000 = 0
# min_coverage = 0
# output_full = FALSE
# 
# rm(catchment,outcome_variable,island,alt_cr,commo,x_pya,dynamics,log_prices,yoyg,short_run,imp,distribution,fe,remaining_forest,offset,lag_or_not,controls,interaction_terms ,interacted,spatial_control,pya_ov,illegal,
#    nearest_mill, weights)

make_base_reg <- function(island,
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
                          spatial_control = FALSE, # logical, if TRUE, adds ~30min computation. Should the average of neighbors' outcome variable be added in the RHS. 
                          pya_ov = FALSE, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
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
  # add pya outcome variable 
  if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}
  
  # lag controls that are from IBS 
  select_ibs_controls <- grepl("pct_own", controls) | grepl("prex_", controls)
  if(lag_or_not=="_lag1" & length(controls)>0 & any(select_ibs_controls)){
    # find them
    controls[select_ibs_controls] <- sapply(controls[select_ibs_controls], FUN = paste0, lag_or_not)
  }
  
  # add spatial control or not
  if(spatial_control){controls <- c(controls, "ngb_ov_lag4")}
  
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
  
  ## Keep observations that: 
  
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
  if(grepl("_slow_",outcome_variable)){
    d <- dplyr::filter(d, year < 2011)
    d <- d[d$year<2011,]
  }
  
  # remove year 2015 as it is only available for industrial plantations
  d <- dplyr::filter(d, year < 2015)
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
  
  
  # should be applied to d_clean rather no ? or at least d_nona...
  # to speed up computation, and also because we do not want to control for values in neighboring cells that are not entering the regressions
  # these values are not bringing any bias by definition as they are excluded
  if(spatial_control){
    
    # restrict the data set to years where the paste outcome is available
    d <- d[d$year > 2004,]
    
    # d is a balanced panel so far so we can take any first record of a parcel to 
    # keep the most general cross section
    nrow(d) == length(unique(d$parcel_id))*length(unique(d$year))
    d_cs <- d[!duplicated(d$parcel_id),c("parcel_id", "year", "lat", "lon", "island",outcome_variable)]
    
    # spatial
    d_cs <- st_as_sf(d_cs, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
    
    # identify neighbors
    # this definition of neighbors includes the 8 closest, surounding, grid cells. 
    d_buf <- st_buffer(d_cs, dist = parcel_size - 10)
    row.names(d_buf) <- d_buf$parcel_id
    sgbp <- st_intersects(d_buf)
    
    
    # iterate on parcel ~30 minutes
    d[,"ngb_ov_lag4"] <- rep(NA, nrow(d))
    for(p_i in 1:length(sgbp)){
      # remove own index
      p_i_neighbors <- sgbp[[p_i]][sgbp[[p_i]] != p_i]
      # get the parcel_id corresponding to the sgbp index
      p_i_neighbors <- as.numeric(attr(sgbp,"region.id")[p_i_neighbors])
      # compute mean of outcome variable lags in the PANEL dataset
      for(y in unique(d$year)){
        d[d$parcel_id == as.numeric(attr(sgbp,"region.id")[p_i]) & d$year == y, "ngb_ov_lag4"]  <- mean(d[d$parcel_id %in% p_i_neighbors & d$year == y, paste0(outcome_variable,"_lag4")],na.rm = TRUE)
      }
    } 
    
    # empty neighbor sets (those that have no neighbors but themselves) and  end up with NaN as mean(integer(0)) 
    # turn them into NA
    d[is.nan(d$ngb_ov_lag3),"ngb_ov_lag3"] <- NA
    d[is.nan(d$ngb_ov_lag4),"ngb_ov_lag4"] <- NA
    
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


#### DESCRIPTIVE STATISTICS #### 
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
  sample_2 <- sample_2[-obs2remove(fml = as.formula(paste0(dep_var_2, " ~ parcel_id + district_year")),
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
  
  sample_1 <- left_join(sample_1[,c("parcel_id","year")], d_all[,c("parcel_id","year",variables_1)], by = c("parcel_id","year"))#
  sample_2 <- left_join(sample_2[,c("parcel_id","year")], d_all[,c("parcel_id","year",variables_2)], by = c("parcel_id","year"))
  
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
  
  length(unique(sample_1$parcel_id)) %>% paste0(" number of grid cells in sample") %>%  print()
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
  
  length(unique(sample_2$parcel_id)) %>% paste0(" total number of grid cells") %>%  print()
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
                     "# grid cells = 8368 \n # grid cell-year = 91620" = 3, 
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
                   strikeout = F) %>% 
  # pack_rows("Ownership signals (%)", 3, 5, 
  #           italic = TRUE)  %>%
  # column_spec(column = c(2,3,5,6,8,9),
  #             width = "3em") %>%
  # column_spec(column = c(4,7,10),
  #             width = "9em") %>%
  
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

#### Table of deforestation in different catchment radius / sample ####

### TOTAL DEFORESTATION 
rows <- c("Sumatra", "Kalimantan", "both")
cols <- c("sample", "30km", "50km", "total")

accu_lucpfp <- matrix(ncol = length(cols), nrow = length(rows)) # 4 cols for sample, 30km, 50km and total. 3 rows for Sumatra, Kalimantan, and both. 
row.names(accu_lucpfp) <- rows
colnames(accu_lucpfp) <- cols

## Prepare polygons of three Indonesian islands of interest 
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"
island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)


#island <- "Sumatra"
for(island in c("Sumatra", "Kalimantan")){
  # LUCPFIP
  defo_dyna <- list()
  for(DYNA in c("rapid", "slow")){
    brick_lucpfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",DYNA,"_",island,"_",parcel_size/1000,"km_total.tif")))
    
    # select 2001-2015 layers 
    layer_names <- paste0("parcel_lucpfip_",DYNA,"_",island,"_",parcel_size/1000,"km_total.",c(1:14))
    
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
    
    # select 2001-2015 layers 
    layer_names <- paste0("parcel_lucpf",size,"p_",island,"_",parcel_size/1000,"km_total.",c(1:14))
    
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
accu_lucpfp["Sumatra", "30km"] <- dplyr::filter(d_30_suma, n_reachable_ibs > 0 & year > 2000 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Sumatra", "50km"] <- dplyr::filter(d_50_suma, n_reachable_ibs > 0 & year > 2000 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Kalimantan", "30km"] <- dplyr::filter(d_30_kali, n_reachable_ibs > 0 & year > 2000 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Kalimantan", "50km"] <- dplyr::filter(d_50_kali, n_reachable_ibs > 0 & year > 2000 & year < 2015)$lucpfap_pixelcount %>% sum()

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

# Finally, add figures from Austin et al. 2017 (SI) - no because not comparable enough because we had smallholders without breaking down here. 
# austin <- matrix(nrow = 2, ncol = ncol(accu_lucpfp))
# austin["Sumatra", 4] <- "(448)"
# austin["Kalimantan", 4] <- "(1021)"

#### Print the LateX table code ####
options(knitr.table.format = "latex") 

colnames(accu_lucpfp) <- NULL

kable(accu_lucpfp, booktabs = T, align = "c", 
      caption = "Deforestation accumulated over 2001-2014, in kha.") %>% 
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
              width = "5em") %>% 
  footnote(general = c(""),
           threeparttable = TRUE, 
           escape = TRUE) 


### INDUSTRIAL VS SMALLHOLDERS
res_data_list_des <- list()
for(SIZE in c("i", "sm")){
  res_data_list_des[[paste0("both_",SIZE)]] <- make_base_reg(island = "both",
                                                             outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                                             output_full = FALSE)
}
# # export for analysis on other plateform
# d_clean_sf <- st_as_sf(d_clean, coords = c("lon","lat"), crs = 4326)
# st_write(d_clean_sf, file.path("temp_data/actual_analysis_data_Kalimantan_a"), driver = "ESRI Shapefile", delete_dsn = TRUE)



des_table <- make_desstats(sample_1 = res_data_list_des[["both_i"]][[2]],
                           sample_2 = res_data_list_des[["both_sm"]][[2]])
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
                     "# grid cells = 8368 \n # grid cell-year = 91620" = 3, 
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
                   strikeout = F) %>% 
  # pack_rows("Ownership signals (%)", 3, 5, 
  #           italic = TRUE)  %>%
  # column_spec(column = c(2,3,5,6,8,9),
  #             width = "3em") %>%
  # column_spec(column = c(4,7,10),
  #             width = "9em") %>%
  
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 


#### DESCRIPTIVE MAP #####
res_data <- make_base_reg(island = "both",
                          outcome_variable = paste0("lucpfap_pixelcount"))

d_clean <- res_data[[2]]
d_clean_cs <- d_clean[!duplicated(d_clean$parcel_id),]

d_cs <- st_as_sf(d_clean_cs, coords = c("lon", "lat"), crs = 4326)
d_cs <- st_transform(d_cs, crs = indonesian_crs)

d_cs <- st_buffer(d_cs, dist = 1500)
st_geometry(d_cs) <- sapply(st_geometry(d_cs), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
d_geo <- st_union(st_geometry(d_cs))
d_geo <- st_transform(d_geo, crs = 4326)


# d_cs <- ddply(d_clean, "parcel_id", summarise, 
#               accu_lucfp = sum(lucpfap_pixelcount))
# 
# d_cs <- left_join(d_cs, d_clean_cs[,c("parcel_id", "lon", "lat")], by = "parcel_id")





### MILLs
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# keep only those that are matched with UML. 
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


leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$CartoDB.DarkMatterNoLabels) %>% # CartoDB.PositronNoLabels  
  setView(lat = -0.493, 
          lng = 107.367, 
          zoom = 4.5) %>% 
  addPolygons(data = island_sf, stroke = FALSE, fill = TRUE, fillColor = "grey", fillOpacity = 0.5) %>% 
  addPolygons(data = st_geometry(d_geo), 
              stroke = TRUE, opacity = 2,  weight = 2, color = "green") %>%   #~cb_ind(d_cs$accu_lucfp)
  addCircleMarkers(data = ibs, radius = 0.001, fillOpacity = 1, fillColor = "red", stroke = FALSE, weight = 0) 


rm(d_clean_cs, d_cs, d_geo, ibs)

#### APE FUNCTION ####
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

res_data <- res_data_list_interact[[1]]
k = 1
K=1
controls_pe = F
SE = "cluster"
CLUSTER = "subdistrict"
rel_price_change = 0.01 # sd/m #
abs_price_change = 1
stddev = FALSE
rounding = 2

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



#### MAIN : INDUS SM LEGAL AND ILLEGAL and EQUALITY TESTS ####

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
  pack_rows("Price elasticity", 1, 2, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = ncol(ape_mat)+1,
              bold = TRUE)

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
  pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "1em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(ape_mat))),
              width = "7em",
              latex_valign = "b") %>% 
  column_spec(column = ncol(ape_mat)+1,
              bold = TRUE)



### COMPARE GROUPS --------------------------------------------------------------------------------
# this is a simpler version of make_APEs function above. 
make_APEs_1regr <- function(res_data, 
                            #SE = "cluster", 
                            stddev = TRUE,
                            CLUSTER = "subdistrict", 
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
    ape_fml_roi <- paste0("((",1+rel_price_change,")^(",linear_ape_fml,") - 1)")#*fv_bar*",pixel_area_ha)
  } else{
    ape_fml_roi <- paste0("(exp(",linear_ape_fml,"*",abs_price_change,") - 1)")#*fv_bar*",pixel_area_ha)
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
