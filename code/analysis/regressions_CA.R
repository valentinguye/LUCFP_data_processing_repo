

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
                   "knitr", "kableExtra",
                   "msm", "car", "fixest", "sandwich", "lmtest", "boot", "multcomp",
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


# # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("leaflet", "leaflet.providers", "png")
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

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

# PIXEL AREA
# to rescale the average partial effects and the predictions from pixel counts to hectares
pixel_area_ha <- (27.8*27.6)/(1e4)

### SET NUMBER OF THREADS USED BY {FIXEST} TO ONE (TO REDUCE R SESSION CRASH)
getFixest_nthreads()

### Set dictionary for names of variables to display in regression tables 
setFixest_dict(c(parcel_id = "grid cell",
                 lucpfap_pixelcount = "All",
                 lucpfip_ha_total = "LUCPFIP (ha)", 
                 lucpfip_pixelcount_total = "Land use change from primary forest to industrial oil palm plantations", 
                 lucpfsmp_ha_total = "LUCPFSMP (ha)", 
                 lucpfsmp_pixelcount_total = "Land use change from primary forest to small or medium-sized oil palm plantations", 
                 lucfip_ha_30th = "LUCFIP (30 pct. canopy density, ha)",
                 lucfip_ha_60th = "LUCFIP (60 pct. canopy density, ha)",
                 lucfip_ha_90th = "LUCFIP (90 pct. canopy density, ha)",
                 lucfip_pixelcount_30th = "LUCFIP (30 pct. canopy density, pixels)",
                 lucfip_pixelcount_60th = "LUCFIP (60 pct. canopy density, pixels)",
                 lucfip_pixelcount_90th = "LUCFIP (90 pct. canopy density, pixels)",
                 # No time dynamics FFB variables 
                 wa_ffb_price_imp1_3ya = "FFB price signal, 3 year average",
                 wa_ffb_price_imp1_3ya_lag1 = "FFB price signal, 3 year average",#(lagged)
                 wa_ffb_price_imp1_4ya = "FFB price signal, 4 year average",
                 wa_ffb_price_imp1_4ya_lag1 = "FFB price signal, 4 year average",#(lagged)
                 wa_ffb_price_imp1_5ya = "FFB price signal, 5 year average",
                 wa_ffb_price_imp1_5ya_lag1 = "FFB price signal, 5 year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_3ya = "FFB price signal y-o-y growth rate, 3 year average",
                 wa_ffb_price_imp1_yoyg_3ya_lag1 = "FFB price signal y-o-y growth rate, 3 year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_4ya = "FFB price signal y-o-y growth rate, 4 year average",
                 wa_ffb_price_imp1_yoyg_4ya_lag1 = "FFB price signal y-o-y growth rate, 4 year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_5ya = "FFB price signal y-o-y growth rate, 5 year average",
                 wa_ffb_price_imp1_yoyg_5ya_lag1 = "FFB price signal y-o-y growth rate, 5 year average",#(lagged)
                 ln_wa_ffb_price_imp1_3ya = "FFB price signal, 3 year average",
                 ln_wa_ffb_price_imp1_3ya_lag1 = "FFB price signal, 3 year average",#(lagged)
                 ln_wa_ffb_price_imp1_4ya = "FFB price signal, 4 year average",
                 ln_wa_ffb_price_imp1_4ya_lag1 = "FFB price signal, 4 year average",#(lagged)
                 ln_wa_ffb_price_imp1_5ya = "FFB price signal, 5 year average",
                 ln_wa_ffb_price_imp1_5ya_lag1 = "FFB price signal, 5 year average",#(lagged)
                 # No time dynamics CPO variables 
                 wa_cpo_price_imp1_3ya = "CPO price signal, 3 year average",
                 wa_cpo_price_imp1_3ya_lag1 = "CPO price signal, 3 year average",#(lagged)
                 wa_cpo_price_imp1_4ya = "CPO price signal, 4 year average",
                 wa_cpo_price_imp1_4ya_lag1 = "CPO price signal, 4 year average",#(lagged)
                 wa_cpo_price_imp1_5ya = "CPO price signal, 5 year average",
                 wa_cpo_price_imp1_5ya_lag1 = "CPO price signal, 5 year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_3ya = "CPO price signal y-o-y growth rate, 3 year average",
                 wa_cpo_price_imp1_yoyg_3ya_lag1 = "CPO price signal y-o-y growth rate, 3 year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_4ya = "CPO price signal y-o-y growth rate, 4 year average",
                 wa_cpo_price_imp1_yoyg_4ya_lag1 = "CPO price signal y-o-y growth rate, 4 year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_5ya = "CPO price signal y-o-y growth rate, 5 year average",
                 wa_cpo_price_imp1_yoyg_5ya_lag1 = "CPO price signal y-o-y growth rate, 5 year average",#(lagged)
                 ln_wa_cpo_price_imp1_3ya = "CPO price signal, 3 year average",
                 ln_wa_cpo_price_imp1_3ya_lag1 = "CPO price signal, 3 year average",#(lagged)
                 ln_wa_cpo_price_imp1_4ya = "CPO price signal, 4 year average",
                 ln_wa_cpo_price_imp1_4ya_lag1 = "CPO price signal, 4 year average",#(lagged)
                 ln_wa_cpo_price_imp1_5ya = "CPO price signal, 5 year average",
                 ln_wa_cpo_price_imp1_5ya_lag1 = "CPO price signal, 5 year average",#(lagged)
                 # SR FFB variables
                 wa_ffb_price_imp1 = "FFB price signal",
                 wa_ffb_price_imp1_lag1 = "FFB price signal",#(lagged)
                 wa_ffb_price_imp1_yoyg = "FFB price signal y-o-y growth rate",
                 wa_ffb_price_imp1_yoyg_lag1 = "FFB price signal y-o-y growth rate",#(lagged)
                 wa_ffb_price_imp1_dev_2pya = "FFB price signal deviation from 2 past year average",
                 wa_ffb_price_imp1_dev_2pya_lag1 = "FFB price signal deviation from 2 past year average",#(lagged)
                 wa_ffb_price_imp1_dev_3pya = "FFB price signal deviation from 3 past year average",
                 wa_ffb_price_imp1_dev_3pya_lag1 = "FFB price signal deviation from 3 past year average",#(lagged)
                 wa_ffb_price_imp1_dev_4pya = "FFB price signal deviation from 4 past year average",
                 wa_ffb_price_imp1_dev_4pya_lag1 = "FFB price signal deviation from 4 past year average",#(lagged)
                 ln_wa_ffb_price_imp1 = "FFB price signal",
                 ln_wa_ffb_price_imp1_lag1 = "FFB price signal",#(lagged)
                 # LR FFB variables
                 wa_ffb_price_imp1_2pya = "FFB price signal, 2 past year average",
                 wa_ffb_price_imp1_2pya_lag1 = "FFB price signal, 2 past year average",#(lagged)
                 wa_ffb_price_imp1_3pya = "FFB price signal, 3 past year average",
                 wa_ffb_price_imp1_3pya_lag1 = "FFB price signal, 3 past year average",#(lagged)
                 wa_ffb_price_imp1_4pya = "FFB price signal, 4 past year average",
                 wa_ffb_price_imp1_4pya_lag1 = "FFB price signal, 4 past year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_2pya = "FFB price signal y-o-y growth rate, 2 past year average",
                 wa_ffb_price_imp1_yoyg_2pya_lag1 = "FFB price signal y-o-y growth rate, 2 past year average",#(lagged)                
                 wa_ffb_price_imp1_yoyg_3pya = "FFB price signal y-o-y growth rate, 3 past year average",
                 wa_ffb_price_imp1_yoyg_3pya_lag1 = "FFB price signal y-o-y growth rate, 3 past year average",#(lagged)
                 wa_ffb_price_imp1_yoyg_4pya = "FFB price signal y-o-y growth rate, 4 past year average",
                 wa_ffb_price_imp1_yoyg_4pya_lag1 = "FFB price signal y-o-y growth rate, 4 past year average",#(lagged)
                 ln_wa_ffb_price_imp1_2pya = "FFB price signal, 2 past year average",
                 ln_wa_ffb_price_imp1_2pya_lag1 = "FFB price signal, 2 past year average",#(lagged)
                 ln_wa_ffb_price_imp1_3pya = "FFB price signal, 3 past year average",
                 ln_wa_ffb_price_imp1_3pya_lag1 = "FFB price signal, 3 past year average",#(lagged)
                 ln_wa_ffb_price_imp1_4pya = "FFB price signal, 4 past year average",
                 ln_wa_ffb_price_imp1_4pya_lag1 = "FFB price signal, 4 past year average",#(lagged)                 
                 # SR CPO variables
                 wa_cpo_price_imp1 = "CPO price signal",
                 wa_cpo_price_imp1_lag1 = "CPO price signal",#(lagged)
                 wa_cpo_price_imp1_yoyg = "CPO price signal y-o-y growth rate",
                 wa_cpo_price_imp1_yoyg_lag1 = "CPO price signal y-o-y growth rate",#(lagged)
                 wa_cpo_price_imp1_dev_2pya = "CPO price signal deviation from 2 past year average",
                 wa_cpo_price_imp1_dev_2pya_lag1 = "CPO price signal deviation from 2 past year average",#(lagged) 
                 wa_cpo_price_imp1_dev_3pya = "CPO price signal deviation from 3 past year average",
                 wa_cpo_price_imp1_dev_3pya_lag1 = "CPO price signal deviation from 3 past year average",#(lagged) 
                 wa_cpo_price_imp1_dev_4pya = "CPO price signal deviation from 4 past year average",
                 wa_cpo_price_imp1_dev_4pya_lag1 = "CPO price signal deviation from 4 past year average",#(lagged) 
                 ln_wa_cpo_price_imp1 = "CPO price signal",
                 ln_wa_cpo_price_imp1_lag1 = "CPO price signal",#(lagged)
                 # LR CPO variables
                 wa_cpo_price_imp1_2pya = "CPO price signal, 2 past year average",
                 wa_cpo_price_imp1_2pya_lag1 = "CPO price signal, 2 past year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_2pya = "CPO price signal y-o-y growth rate, 2 past year average",
                 wa_cpo_price_imp1_yoyg_2pya_lag1 = "CPO price signal y-o-y growth rate, 2 past year average",#(lagged)
                 wa_cpo_price_imp1_3pya = "CPO price signal, 3 past year average",
                 wa_cpo_price_imp1_3pya_lag1 = "CPO price signal, 3 past year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_3pya = "CPO price signal y-o-y growth rate, 3 past year average",
                 wa_cpo_price_imp1_yoyg_3pya_lag1 = "CPO price signal y-o-y growth rate, 3 past year average",#(lagged)
                 wa_cpo_price_imp1_4pya = "CPO price signal, 4 past year average",
                 wa_cpo_price_imp1_4pya_lag1 = "CPO price signal, 4 past year average",#(lagged)
                 wa_cpo_price_imp1_yoyg_4pya = "CPO price signal y-o-y growth rate, 4 past year average",
                 wa_cpo_price_imp1_yoyg_4pya_lag1 = "CPO price signal y-o-y growth rate, 4 past year average",#(lagged)
                 ln_wa_cpo_price_imp1_2pya = "CPO price signal, 2 past year average",
                 ln_wa_cpo_price_imp1_2pya_lag1 = "CPO price signal, 2 past year average",#(lagged)
                 ln_wa_cpo_price_imp1_3pya = "CPO price signal, 3 past year average",
                 ln_wa_cpo_price_imp1_3pya_lag1 = "CPO price signal, 3 past year average",#(lagged)
                 ln_wa_cpo_price_imp1_4pya = "CPO price signal, 4 past year average",
                 ln_wa_cpo_price_imp1_4pya_lag1 = "CPO price signal, 4 past year average",#(lagged)
                 ## interactions 
                 n_reachable_uml_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "reachable mills X FFB price signal",
                 n_reachable_uml_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "reachable mills X CPO price signal",
                 wa_pct_own_loc_gov_imp_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "local gvt mill ownership X FFB price signal",
                 wa_pct_own_loc_gov_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "local gvt mill ownership X CPO price signal",
                 wa_pct_own_nat_priv_imp_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "private mill ownership X FFB price signal",
                 wa_pct_own_nat_priv_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "private mill ownership X CPO price signal",
                 wa_pct_own_for_imp_lag1Xln_wa_ffb_price_imp1_4ya_lag1 = "foreign mill ownership X FFB price signal",
                 wa_pct_own_for_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1 = "foreign mill ownership X CPO price signal",
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
                 wa_pct_own_cent_gov_imp = "Local government mill ownership (pct.)",
                 wa_pct_own_cent_gov_imp_lag1 = "Local government mill ownership (pct., lagged)",
                 wa_pct_own_loc_gov_imp = "Local government mill ownership (pct.)",
                 wa_pct_own_loc_gov_imp_lag1 = "Local government mill ownership (pct., lagged)",
                 wa_pct_own_nat_priv_imp = "Domestic private mill ownership (pct.)",
                 wa_pct_own_nat_priv_imp_lag1 = "Domestic private mill ownership (pct., lagged)",
                 wa_pct_own_for_imp = "Foreign mill ownership (pct.)",
                 wa_pct_own_for_imp_lag1 = "Foreign mill ownership (pct., lagged)", 
                 wa_prex_cpo_imp1 = "Percentage CPO exported",
                 wa_prex_cpo_imp1_lag1 = "Percentage CPO exported",#(lagged)
                 wa_prex_cpo_imp2 = "Percentage CPO exported",
                 wa_prex_cpo_imp2_lag1 = "Percentage CPO exported"#(lagged)
))



# OK SO NOW FIXED EFFECT VARIABLES ARE DEFINED IN add_parcel_variables.R AS geographic_year 
# the only change implied is that now district^year can also be specified as district_year.

## WITHIN CATCHMENT AREA
# They are outputs of merge_lhs_rhs_CA_parcels.R

# 2 hour catchment area
d_CA2 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                  parcel_size/1000,"km_",
                                  "2h_CA.rds")))

# Split them into islands of interest
d_CA2_suma <- d_CA2[d_CA2$island == "Sumatra",]
d_CA2_kali <- d_CA2[d_CA2$island == "Kalimantan",]

rm(d_CA2)

# 4 hour catchment area
d_CA4 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                  parcel_size/1000,"km_",
                                  "4h_CA.rds")))


# Split them into islands of interest
d_CA4_suma <- d_CA4[d_CA4$island == "Sumatra",]
d_CA4_kali <- d_CA4[d_CA4$island == "Kalimantan",]

rm(d_CA4)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### REGRESSION FUNCTION ##### 

# Prefered specifications are passed as default arguments. 
# Currently, the returned object is a fixest object, of *ONE* regression.  
island = "both"
outcome_variable = paste0("lucpf",SIZE,"p_pixelcount_total")
all_producers = TRUE
alt_ca = FALSE
commo = c("cpo")#"ffb", 
x_pya = 3
dynamics = FALSE
log_prices = TRUE
yoyg = FALSE
short_run = "full"
imp = 1
distribution = "quasipoisson"
fe = "parcel_id + district_year"
lag_or_not = "_lag1"
controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml")#"wa_pct_own_loc_gov_imp",
remaining_forest = FALSE
offset = FALSE
interaction_terms = NULL
interacted = "regressors"
spatial_control = FALSE
pya_ov = FALSE
illegal = "all"
weights = FALSE


make_base_reg <- function(island,
                          outcome_variable = "lucpfip_pixelcount_total", # LHS. One of "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total"
                          all_producers = FALSE, # whether industrial and smallholders plantations should be added. 
                          alt_ca = FALSE, # logical, if TRUE, Sumatra's catchment area is based on 4 hour driving times, and Kalimantan on 2 hours.  
                          commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo"), commodities the price signals of which should be included in the RHS
                          x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                          dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                          yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                          log_prices = TRUE, # Logical, should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
                          short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                          imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                          distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                          fe = "parcel_id + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                          lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                          controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml"), # character vectors of names of control variables (don't specify lags in their names)
                          remaining_forest = FALSE, # Logical. If TRUE, the remaining forest is added as a control
                          offset = FALSE, # Logical. Should the log of the remaining forest be added as an offset.  
                          interaction_terms = controls, # may be one or several of the controls specified above. 
                          interacted = "regressors",
                          spatial_control = FALSE, # logical, if TRUE, adds ~30min computation. Should the average of neighbors' outcome variable be added in the RHS. 
                          pya_ov = FALSE, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                          illegal = "all", # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                          weights = FALSE # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
){
  
  ### BASE DATA 
  
  # Catchment area
  if(island == "Sumatra" & alt_ca == FALSE){
    d <- d_CA2_suma
  }
  if(island == "Sumatra" & alt_ca == TRUE){
    d <- d_CA4_suma
  }
  if(island == "Kalimantan" & alt_ca == FALSE){
    d <- d_CA2_kali
  }
  if(island == "Kalimantan" & alt_ca == TRUE){
    d <- d_CA4_kali
  }
  if(island == "both" & alt_ca == FALSE){
    d <- rbind(d_CA2_suma, d_CA2_kali)#, d_CA4_papu
  }
  if(island == "both" & alt_ca == TRUE){
    d <- rbind(d_CA4_suma, d_CA4_kali)#, d_CA2_papu
  }
  
  # create variable of sum of rapid and slow, in case we want to distinguish we the other measure of lucpfip_pixelcount_total, that might include some replacement within plantations
  if(grepl("rapidslow",outcome_variable)){
    d$lucpfip_rapidslow_pixelcount <- d$lucpfip_rapid_pixelcount + d$lucpfip_slow_pixelcount
  }  
  
  # add pixel counts of lucfp from industrial and from small and medium sized plantations. 
  # note that for some observations, this is a sum of 0 and a positive term, and for some others it's a sum of two positive terms. 
  # In the latter case, this grounds on the assumption that the industrial and the S&M plantation maps do not overlap. 
  # Also, we let the industrial outc
  if(grepl("lucpfap", outcome_variable)){
    if(grepl("rapid", outcome_variable)){
      d$lucpfap_pixelcount <- d$lucpfip_rapid_pixelcount + d$lucpfsmp_pixelcount_total
    }else if (grepl("slow", outcome_variable)){
      d$lucpfap_pixelcount <- d$lucpfip_slow_pixelcount + d$lucpfsmp_pixelcount_total
    }else if (grepl("rapidslow", outcome_variable)){
      d$lucpfap_pixelcount <- d$lucpfip_rapidslow_pixelcount + d$lucpfsmp_pixelcount_total
    }else {
      d$lucpfap_pixelcount <- d$lucpfip_pixelcount_total + d$lucpfsmp_pixelcount_total
    }
    outcome_variable <- "lucpfap_pixelcount"
  }
  # For 30% tree cover threshold forest but wont work for now:
  if(grepl("lucfap", outcome_variable)){
    if(grepl("rapid", outcome_variable)){
      d$lucfap_pixelcount <- d$lucfip_rapid_pixelcount + d$lucfsmp_pixelcount_total
    }else if (grepl("slow", outcome_variable)){
      d$lucfap_pixelcount <- d$lucfip_slow_pixelcount + d$lucfsmp_pixelcount_total
    }else if (grepl("rapidslow", outcome_variable)){
      d$lucfap_pixelcount <- d$lucfip_rapidslow_pixelcount + d$lucfsmp_pixelcount_total
    }else {
      d$lucfap_pixelcount <- d$lucfip_pixelcount_30th + d$lucfsmp_pixelcount_total
    }
  }
  
  
  ### SPECIFICATIONS  
  
  ## OFFSET
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
  
  # if, on the other hand, we want to identify the short run effects, controlling for the medium run's. 
  if(dynamics == TRUE){ 
    
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                        paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
      }
      if(yoyg == FALSE){
        if(short_run == "full"){
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
        if(short_run == "full"){
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
  
  # lag controls or not
  if(lag_or_not=="_lag1"){controls <- paste0(controls,lag_or_not)}
  
  # add spatial control or not
  if(spatial_control){controls <- c(controls, "ngb_ov_lag4")}
  
  # add remaining forest or not
  if(remaining_forest){
    if(grepl("lucpf", outcome_variable)){
      controls <- c(controls, "remain_pf_pixelcount")
    }
    if(grepl("lucf", outcome_variable)){
      controls <- c(controls, "remain_f30th_pixelcount")
    }
    
    offset <- FALSE
  }
  
  
  ### SELECT DATA FOR REGRESSION
  
  ## group all the variables necessary in the regression
  # important to do that after outcome_variable, regressors controls etc. have been (re)defined. 
  # (interactions do not need to be in there as they are fully built from the used_vars)
  used_vars <- c(outcome_variable, regressors, controls,
                 "parcel_id", "year", "lat", "lon", "district", "province", "island", "district_year", "province_year")
                 #"n_reachable_ibsuml_lag1", "sample_coverage_lag1", #"pfc2000_total_ha", 
                 #"remain_f30th_pixelcount","remain_pf_pixelcount"
  
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
  # (so remove obs. from all years of grid cells that are not forested in 2000, or years after grid cells get completely deforested)
  if(grepl("lucpf", outcome_variable)){
    d <- d[d$remain_pf_pixelcount > 0,]
  }
  
  if(grepl("lucf", outcome_variable)){
    d <- d[d$remain_f30th_pixelcount > 0,]
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
  
  
  # - are not in an intended RSPO-certified supply base
  d <- d[d$rspo_cert==FALSE,]
  
  # are legal, illegal, or neither
  if(illegal == "no_ill1"){
    d <- d[d$illegal1 == FALSE, ]
  }
  if(illegal == "no_ill2"){
    d <- d[d$illegal2 == FALSE, ]
  }
  if(illegal == "ill1"){
    d <- d[d$illegal1 == TRUE, ]
  }  
  if(illegal == "ill2"){
    d <- d[d$illegal2 == TRUE, ]
  }
  
  # - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
  usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
  names(usable) <- used_vars            
  usable <- bind_cols(usable)
  filter_vec <- base::rowSums(usable)
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d[filter_vec, c(used_vars)]
  if(anyNA(d_nona)){stop()}
  rm(filter_vec, usable)
  
  # - sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
  # probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
  # just remove them if any, so that there is no bug. 
  if(weights == TRUE){
    d_nona <- d_nona[d_nona$sample_coverage_lag1!=0,]
  }
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_clean <- d_nona[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                                d_nona, 
                                family = "poisson"),]
  
  # note that obs2remove has to be the last filtering operation on data, otherwise some units (along fe dimension)
  # may become "only-zero-outcome" units after other data removal.  
  
  
  
  # make the interaction variables in the data
  
  ## INTERACTIONS
  # here we produce the names of the actual interaction variables
  if(length(interaction_terms)>0){
    # unless variables of interest to be interacted are specified, they are all the presently defined regressors
    if(interacted == "regressors"){interacted <- regressors}
    # this function selects the actual controls (lagged if necessary) into the interaction terms
    make_int_term <- function(x){int_term <- controls[grepl(x,controls)] %>% paste0("X",interacted)}
    interaction_vars <- sapply(interaction_terms, FUN = make_int_term)
  }else{interacted <- NULL}
 
  # and here we build the interactions 
  if(length(interaction_terms)>0){
    for(con in interaction_terms){
      actual_ <- controls[grepl(con, controls)]
      for(reg in interacted){
        d_clean[,paste0(actual_,"X",reg)] <- d_clean[,actual_]*d_clean[,reg]
      }
    }
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
  
  if(length(interaction_terms)>0){
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ ",
                                  paste0(regressors, collapse = "+"),
                                  " + ",
                                  paste0(controls, collapse = "+"),
                                  " + ",
                                  paste0(interaction_vars, collapse = "+"),
                                  " | ",
                                  fe))
  }
  
  
  # Estimation
  if(offset == TRUE){
    if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
      if(weights == TRUE){
        var_weights <- d_clean$sample_coverage_lag1/100 
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution,
                                 offset = offset_fml,
                                 glm.iter = 100,
                                 notes = TRUE, 
                                 weights = var_weights)
        
      }else{
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution,
                                 offset = offset_fml,
                                 glm.iter = 100,
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
        var_weights <- d_clean$sample_coverage_lag1/100 
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution, 
                                 glm.iter = 50,
                                 notes = TRUE, 
                                 weights = var_weights)
        
      }else{
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution, 
                                 glm.iter = 50,
                                 notes = TRUE)
      }
    }else{ # no weights allowed in negative binomial
      reg_res <- fixest::fenegbin(fe_model,
                                  data = d_clean, 
                                  notes = TRUE)
    }
  }
  
  rm(d, d_nona)
  toreturn <- list(reg_res, d_clean)
  return(toreturn)
  
}



#### MAIN - RUN REGRESSIONS #### 
# pass arguments from make_base_reg in lapply to get alternative specifications

isl_list <- list("Sumatra", "Kalimantan", "both")


## Elasticities 

# Industrial
res_data_list_ind <- lapply(isl_list, make_base_reg,
                           outcome_variable = "lucpfip_pixelcount_total",
                           interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                           commo = c("cpo"), #
                           offset = FALSE)

names(res_data_list_ind) <- paste0(isl_list, " industrial")
res_list_ind <- lapply(res_data_list_ind, FUN = function(x){x[[1]]})

# Smallholders
res_data_list_sm <- lapply(isl_list, make_base_reg,
                          outcome_variable = "lucpfsmp_pixelcount_total",
                          interaction_terms =  NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp"),
                          commo = c("cpo"),#
                          offset = FALSE) 

names(res_data_list_sm) <- paste0(isl_list, " smallholders")
res_list_sm <- lapply(res_data_list_sm, FUN = function(x){x[[1]]})

# ALL
res_data_list_all <- make_base_reg(island = "all", 
                                   all_producers = TRUE,
                                   outcome_variable = "lucpfip_pixelcount_total", # with all_producers = TRUE, this indicates whether we want primary forest or not, and also which dynamic in industrial lucfp we want. 
                                   interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" "wa_pct_own_nat_priv_imp"),
                                   commo = c("cpo"), #
                                   offset = FALSE) %>% list() # make it a list so that it's at the same level as the others below

names(res_data_list_all) <- "all" 
#res_list_all <- lapply(res_data_list_all, FUN = function(x){x[[1]]})
res_list_all <- res_data_list_all[1][[1]][1]


# For display of elasticities
res_list <- c(res_list_ind, res_list_sm, res_list_all)
rm(res_list_all, res_list_ind, res_list_sm)

# For all the computations of partial effects based on these computations: 
res_data_list <- c(res_data_list_ind, res_data_list_sm, res_data_list_all)


# preview in R
etable(res_list,
       se = "cluster",
       #coefstat = "confint",
       subtitles = c("Sumatra industrial", "Kalimantan industrial", "All industrial", "Sumatra smallholders", "Kalimantan smallholders", "All smallholders", "All"))


### LATEX

#table_title_dyn <- paste0("LUCFP semi-elasticities to medium-run price signals") 
table_title_dyn <- paste0("LUCFP elasticities to medium-run price signals") 
#table_title_dyn <- paste0("LUCFP semi-elasticities to medium-run y-o-y growth rates of price signals") 

etable(res_list, 
       #cluster = oneway_cluster,
       se = "cluster",
       tex = TRUE,
       # file = table_file, 
       # replace = TRUE,
       title = table_title_dyn,
       subtitles = c(" ", "Sumatra", "Kalimantan", "both","Sumatra", "Kalimantan", "both"), # first empty subtitle is for the overall (islands and sizes) column 
       #family = TRUE,
       #drop = c("own", "reachable"),
       #coefstat = "confint",
       sdBelow = TRUE,
       yesNo = "X",
       fitstat = c("sq.cor"),
       convergence = TRUE,
       dict = TRUE, 
       style=list(model = "title:"),
       powerBelow = -7)

#rm(res_list)




#### MAIN - MAKE APEs ####
res_data <- res_data_list[[1]]

# helper function that transforms the list of results into a data frame of average partial effects (APEs) and their standard errors (SEs)
make_APEs <- function(res_data){
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]]
  
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
  
  # average of the regressor of interest
  reg_bar <- mean(d_clean[,coeff_names[1]]) 
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values)
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  linear_ape_fml <- paste0(coeff_names[1])
  i_t <- 1
  while(i_t<=length(interaction_terms)){
    linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[i_t],"*",int_term_avg[[i_t]])
    i_t <- i_t +1
  }
  ape_fml_roi <- paste0("(",linear_ape_fml,")*fv_bar*",pixel_area_ha)
  
  dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                            vcov. = vcov(reg_res, se = "cluster"), 
                            g. = ape_fml_roi, 
                            rhs = 0)
  
  row.names(dM_ape_roi) <- coeff_names[1]
  dM_ape_roi <- as.matrix(dM_ape_roi)
  dM_ape_roi <- dM_ape_roi[,c(1,2,7)]

  ## FORMULA FOR APE OF INTERACTION TERMS 
  dM_ape_list <- list()
  dM_ape_list[[1]] <- dM_ape_roi
  if(length(interaction_terms) >0){
    for(i in 1:length(interaction_terms)){
      
      # there are several parts 
      # 1. the coeff of the interaction effect of interest
      iei <- paste0(interaction_terms[i], "X", coeff_names[1])
      ape_fml1 <- paste0(iei," + ") 
  
      # 2. the linear expression derivative wrt. the regressor of interest 
      ape_fml2 <- paste0("(",linear_ape_fml,")")
      
      # 3. the derivative of the linear expression wrt. the interaction term 
      ape_fml3 <- paste0("*(",
                         coeff_names[coeff_names==interaction_terms[i]]," + ",
                         iei,"*",reg_bar,")")#
      
      # paste everything together
      ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*",fv_bar,"*",pixel_area_ha)
      
      dM <- deltaMethod(object = coef(reg_res), 
                                    vcov. = vcov(reg_res, se = "cluster"), 
                                    g. = ape_fml_it, 
                                    rhs = 0)
      
      row.names(dM) <- interaction_terms[i]
      dM <- as.matrix(dM)
      dM <- dM[,c(1,2,7)]
      dM_ape_list[[i+1]] <- dM
      
    }
  }
  mat <- matrix(ncol = 1, 
                nrow = 3+3*length(interaction_terms), 
                data = unlist(dM_ape_list))  
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, ape_fml_it, dM_ape_list)
  return(mat)
}

rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list, FUN = make_APEs)) %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","SE","p-value"), nrow(ape_mat)/3), "Observations")
# keep ape_mat like this for later comparisons between APEs

# prepare ape_mat for kable
ape_disp <- ape_mat %>% round(digits = 3)
ape_disp[nrow(ape_disp),] <- ape_disp[nrow(ape_disp),] %>% formatC(digits = 0, format = "f")
colnames(ape_disp) <- NULL
ape_disp

options(knitr.table.format = "latex")
kable(ape_disp, booktabs = T, align = "r",
      caption = "Average partial effects on tree replacement (ha)") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "both" = 1#,
                     # "Sumatra" = 1,
                     # "Kalimantan" = 1,
                     # "both" = 1, 
                     # " " = 1
                     ),
                  bold = F,
                  align = "c") %>%
  add_header_above(c(" " = 1,
                     "Industrial plantations" = 3#,
                     #"Smallholder plantations" = 3, 
                     #"All" = 1
                     ),
                  align = "c",
                  strikeout = F) %>%
  pack_rows("CPO medium-run price", 1, 3, 
            italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n domestic private ownership", 4, 6, # domestic private ownership  "Interaction with \n domestic private ownership"
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n foreign ownership", 7, 9,
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n # of reachable mills", 10, 12,
  #           italic = TRUE, bold = TRUE)  %>%
  pack_rows(start_row =  nrow(df), end_row = nrow(df),  latex_gap_space = "0.5em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "14em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(df))),
              width = "5em",
              latex_valign = "b")

rm(ape_disp)


#### MAIN - COMPARE APEs #### 

# First, rerun everything to be sure that all grousp are on the same specification and also bc 
# we need additional regressions that were not run above, it's the by island - all producers
isl_list <- list("Sumatra", "Kalimantan", "both")

## Elasticities 

# Industrial
res_data_list_ind <- lapply(isl_list, make_base_reg,
                            outcome_variable = "lucpfip_pixelcount_total",
                            interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                            commo = c("cpo"), #
                            offset = FALSE)

names(res_data_list_ind) <- paste0(isl_list, " industrial")
res_list_ind <- lapply(res_data_list_ind, FUN = function(x){x[[1]]})

# Smallholders
res_data_list_sm <- lapply(isl_list, make_base_reg,
                           outcome_variable = "lucpfsmp_pixelcount_total",
                           interaction_terms =  NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp"),
                           commo = c("cpo"),#
                           offset = FALSE) 

names(res_data_list_sm) <- paste0(isl_list, " smallholders")
res_list_sm <- lapply(res_data_list_sm, FUN = function(x){x[[1]]})

# ALL
res_data_list_all <- lapply(isl_list, make_base_reg,
                            all_producers = TRUE,
                            outcome_variable = "lucpfip_pixelcount_total", # with all_producers = TRUE, this indicates whether we want primary forest or not, and also which dynamic in industrial lucfp we want. 
                            interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" "wa_pct_own_nat_priv_imp"),
                            commo = c("cpo"), #
                            offset = FALSE) 

names(res_data_list_all) <- paste0(isl_list, " all") 
res_list_all <- lapply(res_data_list_all, FUN = function(x){x[[1]]})

# For all the computations of partial effects based on these computations: 
res_data_list <- c(res_data_list_ind, res_data_list_sm, res_data_list_all)


rm(ape_mat)
ape_mat <- bind_cols(lapply(res_data_list, FUN = make_APEs)) %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","SE","p-value"), nrow(ape_mat)/3), "Observations")
# keep ape_mat like this for comparisons between APEs


compare_APEs_across_groups <- function(group1, group2, m0 = 0, alpha = 0.05, alternative = "two.sided") { 
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

comp_ape_mat <- matrix(ncol = 5, nrow = 2, data = NA)
colnames(comp_ape_mat) <- c("all", "Sumatra", "Kalimantan", "industrial", "smallholders") 
row.names(comp_ape_mat) <- c("Sumatra = Kalimantan", "industrial = smallholders")

for(SIZE in c("all", "industrial", "smallholders")){
  comp_ape_mat["Sumatra = Kalimantan",SIZE] <- compare_APEs_across_groups(group1 = paste0("Sumatra ",SIZE), 
                                                                          group2 = paste0("Kalimantan ",SIZE)) 
}
for(ISL in c("all", "Sumatra", "Kalimantan")){
  comp_ape_mat["industrial = smallholders",ISL] <- compare_APEs_across_groups(group1 = paste0(ISL," industrial"), 
                                                                              group2 = paste0(ISL," smallholders")) 
}

comp_ape_disp <- comp_ape_mat

comp_ape_disp <- comp_ape_mat %>% formatC(digits = 4, format = "f")
comp_ape_disp[comp_ape_disp=="   NA"] <- "" 
colnames(comp_ape_disp) <- NULL
comp_ape_disp

options(knitr.table.format = "latex")
kable(comp_ape_disp, booktabs = T, align = "r",
      caption = "Comparisons of CPO medium-run price signal APEs on LUCFP across groups") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Overall" = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "Industrial plantations" = 1,
                     "Smallholder plantations" = 1
                      ),
                    bold = F,
                    align = "c") %>%
  pack_rows("Ho", 1, 2, 
            bold = TRUE)  %>%
  # pack_rows("Interaction with \n domestic private ownership", 4, 6, # domestic private ownership  "Interaction with \n domestic private ownership"
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n foreign ownership", 7, 9,
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n # of reachable mills", 10, 12,
  #           italic = TRUE, bold = TRUE)  %>%
  pack_rows(start_row =  nrow(df), end_row = nrow(df),  latex_gap_space = "0.5em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "14em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(df))),
              width = "5em",
              latex_valign = "b")

rm(ape_disp)

#### LEGAL - RUN REGRESSIONS #### 
# Reiterate steps from tables 1-2 but do not print LateX for table 1. 

# infrastructure to store results
res_data_list_ill <- list()
elm <- 1

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def))

for(ILL in ill_status){
   res_data_list_ill[[elm]] <- make_base_reg(island = "all", 
                                            all_producers = TRUE,
                                            outcome_variable = "lucpfip_rapidslow_pixelcount",
                                            illegal = ILL,
                                            commo = c("cpo"), #"ffb",
                                            interaction_terms = c("n_reachable_uml"),# "wa_pct_own_for_imp", "n_reachable_uml"),
                                            offset = FALSE)
  elm <- elm + 1
}

ov_list <- list("lucpfip_rapidslow_pixelcount", "lucpfsmp_pixelcount_total")
isl_list <- list("Sumatra", "Kalimantan", "both")
for(OV in ov_list){
  for(ISL in isl_list){
    for(ILL in ill_status){
      if(!(OV == "lucpfsmp_pixelcount_total" & ISL == "Kalimantan" & ILL == "ill2")){# because in this case there is not enough variation
             res_data_list_ill[[elm]] <- make_base_reg(island = ISL, 
                                               outcome_variable = OV,
                                               illegal = ILL,
                                               commo = c("cpo"), #"ffb",
                                               interaction_terms = c("n_reachable_uml"),# "wa_pct_own_for_imp", "n_reachable_uml"),
                                               offset = FALSE) 
      }

      elm <- elm + 1
    }
  }
}
# elm should be 15 here

# res_list_ill <- lapply(res_data_list_ill, FUN = function(x){x[[1]]})
# 
# 
# # preview in R
# etable(res_list_ill,
#        se = "cluster",
#        #coefstat = "confint",
#        subtitles = c("All legal", "All illegal", 
#                      "Sumatra industrial legal", "Sumatra industrial illegal", 
#                      "Kalimantan industrial legal", "Kalimantan industrial illegal", 
#                      "All industrial legal", "All industrial illegal", 
#                      "Sumatra smallholders legal", "Sumatra smallholders illegal", 
#                      "Kalimantan smallholders legal", "Kalimantan smallholders illegal", 
#                      "All smallholders legal", "All smallholders illegal"))



#### LEGAL - MAKE APEs ####
res_data <- res_data_list_ill[[1]]
make_APEs <- function(res_data){
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]]
  
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
  
  # average of the regressor of interest
  reg_bar <- mean(d_clean[,coeff_names[1]]) 
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values)
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  linear_ape_fml <- paste0(coeff_names[1])
  i_t <- 1
  while(i_t<=length(interaction_terms)){
    linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[i_t],"*",int_term_avg[[i_t]])
    i_t <- i_t +1
  }
  ape_fml_roi <- paste0("(",linear_ape_fml,")*fv_bar*",pixel_area_ha)
  
  dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                            vcov. = vcov(reg_res, se = "cluster"), 
                            g. = ape_fml_roi, 
                            rhs = 0)
  
  row.names(dM_ape_roi) <- coeff_names[1]
  dM_ape_roi <- as.matrix(dM_ape_roi)
  dM_ape_roi <- dM_ape_roi[,c(1,2,7)]
  
  ## FORMULA FOR APE OF INTERACTION TERMS 
  dM_ape_list <- list()
  dM_ape_list[[1]] <- dM_ape_roi
  if(length(interaction_terms) >0){
    for(i in 1:length(interaction_terms)){
      
      # there are several parts 
      # 1. the coeff of the interaction effect of interest
      iei <- paste0(interaction_terms[i], "X", coeff_names[1])
      ape_fml1 <- paste0(iei," + ") 
      
      # 2. the linear expression derivative wrt. the regressor of interest 
      ape_fml2 <- paste0("(",linear_ape_fml,")")
      
      # 3. the derivative of the linear expression wrt. the interaction term 
      ape_fml3 <- paste0("*(",
                         coeff_names[coeff_names==interaction_terms[i]]," + ",
                         iei,"*",reg_bar,")")#
      
      # paste everything together
      ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*",fv_bar,"*",pixel_area_ha)
      
      dM <- deltaMethod(object = coef(reg_res), 
                        vcov. = vcov(reg_res, se = "cluster"), 
                        g. = ape_fml_it, 
                        rhs = 0)
      
      row.names(dM) <- interaction_terms[i]
      dM <- as.matrix(dM)
      dM <- dM[,c(1,2,7)]
      dM_ape_list[[i+1]] <- dM
      
    }
  }
  mat <- matrix(ncol = 1, 
                nrow = 3+3*length(interaction_terms), 
                data = unlist(dM_ape_list))  
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, ape_fml_it, dM_ape_list)
  return(mat)
}

rm(df)
df <- bind_cols(lapply(res_data_list_ill[lengths(res_data_list_ill)>0], FUN = make_APEs)) %>% as.matrix()
# prepare df for kable
# row.names(df) <- paste0("APE - ",names(res_data_list))


row.names(df) <- c(rep(c("Estimate","SE","p-value"), nrow(df)/3), "Observations")
df <- df %>% round(digits = 3)
df[nrow(df),] <- df[nrow(df),] %>% formatC(digits = 0, format = "f")
colnames(df) <- NULL
df

# remove here the models that did not converge.
# with current specification it's 12th model, Kalimantan smallholders illegal
# df <- df[, c(1:11,13,14)]

options(knitr.table.format = "latex")
kable(df, booktabs = T, align = "r",
      caption = "Average partial effects legal and illegal LUCFP (ha)") %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "legal" = 1, 
                     "illegal" = 1, 
                     "legal" = 1, 
                     "illegal" = 1, 
                     "legal" = 1, 
                     "illegal" = 1, 
                     "legal" = 1, 
                     "illegal" = 1, 
                     "legal" = 1, 
                     "illegal" = 1, 
                     "legal" = 1, 
                     #"illegal" = 1, 
                     "legal" = 1, 
                     "illegal" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     " " = 2,
                     "Sumatra" = 2,
                     "Kalimantan" = 2,
                     "both" = 2,
                     "Sumatra" = 2,
                     "Kalimantan" = 1,
                     "both" = 2),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "All" = 2,
                     "Industrial plantations" = 6,
                     "Smallholder plantations" = 5),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("CPO medium-run price", 1, 3, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction with \n # of reachable mills", 4, 6, # domestic private ownership  "Interaction with \n domestic private ownership"
            italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n foreign ownership", 7, 9,
  #           italic = TRUE, bold = TRUE)  %>%
  # pack_rows("Interaction with \n # of reachable mills", 10, 12,
  #           italic = TRUE, bold = TRUE)  %>%
  pack_rows(start_row =  nrow(df), end_row = nrow(df),  latex_gap_space = "0.5em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "12em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(df))),
              width = "4em",
              latex_valign = "b")

rm(df)

#### LUCFIP DYNAMICS - RUN REGRESSIONS #### 
# Reiterate steps from tables 1-2 but do not print LateX for table 1. 


# infrastructure to store results
res_data_list_dyn <- list()
elm <- 1

ov_list <- list("lucpfip_rapid_pixelcount", "lucpfip_slow_pixelcount")
isl_list <- list("Sumatra", "Kalimantan", "both")
for(OV in ov_list){
  for(ISL in isl_list){
    res_data_list_dyn[[elm]] <- make_base_reg(island = ISL, 
                                             outcome_variable = OV,
                                             commo = c("cpo"), #"ffb",
                                             interaction_terms = c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp", "n_reachable_uml"),
                                             offset = FALSE)
    elm <- elm + 1
  }
}
# elm should be 7 here
elm
names(res_data_list_dyn) <- c("Sumatra rapid", "Kalimantan rapid", "All rapid", "Sumatra slow", "Kalimantan slow", "All slow")

res_list_dyn <- lapply(res_data_list_dyn, FUN = function(x){x[[1]]})

# preview in R
etable(res_list_dyn,
       se = "cluster",
       #coefstat = "confint",
       subtitles = c("Sumatra rapid", "Kalimantan rapid", "All rapid", "Sumatra slow", "Kalimantan slow", "All slow"))


#### LUCFIP DYNAMICS - MAKE APEs #### 
res_data <- res_data_list_dyn[[1]]
make_APEs <- function(res_data){
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]]
  
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
  
  # average of the regressor of interest
  reg_bar <- mean(d_clean[,coeff_names[1]]) 
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values)
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  linear_ape_fml <- paste0(coeff_names[1])
  i_t <- 1
  while(i_t<=length(interaction_terms)){
    linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[i_t],"*",int_term_avg[[i_t]])
    i_t <- i_t +1
  }
  ape_fml_roi <- paste0("(",linear_ape_fml,")*fv_bar*",pixel_area_ha)
  
  dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                            vcov. = vcov(reg_res, se = "cluster"), 
                            g. = ape_fml_roi, 
                            rhs = 0)
  
  row.names(dM_ape_roi) <- coeff_names[1]
  dM_ape_roi <- as.matrix(dM_ape_roi)
  dM_ape_roi <- dM_ape_roi[,c(1,2,7)]
  
  ## FORMULA FOR APE OF INTERACTION TERMS 
  dM_ape_list <- list()
  dM_ape_list[[1]] <- dM_ape_roi
  if(length(interaction_terms) >0){
    for(i in 1:length(interaction_terms)){
      
      # there are several parts 
      # 1. the coeff of the interaction effect of interest
      iei <- paste0(interaction_terms[i], "X", coeff_names[1])
      ape_fml1 <- paste0(iei," + ") 
      
      # 2. the linear expression derivative wrt. the regressor of interest 
      ape_fml2 <- paste0("(",linear_ape_fml,")")
      
      # 3. the derivative of the linear expression wrt. the interaction term 
      ape_fml3 <- paste0("*(",
                         coeff_names[coeff_names==interaction_terms[i]]," + ",
                         iei,"*",reg_bar,")")#
      
      # paste everything together
      ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*",fv_bar,"*",pixel_area_ha)
      
      dM <- deltaMethod(object = coef(reg_res), 
                        vcov. = vcov(reg_res, se = "cluster"), 
                        g. = ape_fml_it, 
                        rhs = 0)
      
      row.names(dM) <- interaction_terms[i]
      dM <- as.matrix(dM)
      dM <- dM[,c(1,2,7)]
      dM_ape_list[[i+1]] <- dM
      
    }
  }
  mat <- matrix(ncol = 1, 
                nrow = 3+3*length(interaction_terms), 
                data = unlist(dM_ape_list))  
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, ape_fml_it, dM_ape_list)
  return(mat)
}

df <- bind_rows(lapply(res_data_list_dyn, FUN = make_APEs)) %>% as.matrix()

# prepare df for kable
row.names(df) <- c(rep(c("Estimate","SE","p-value"), nrow(df)/3), "Observations")
df <- df %>% round(digits = 3)
df[nrow(df),] <- df[nrow(df),] %>% formatC(digits = 0, format = "f")
colnames(df) <- NULL
df


# remove here the models that did not converge.
options(knitr.table.format = "latex")
kable(df, booktabs = T, align = "r",
      caption = "Average partial effects on rapid and slow industrial LUCFP (ha)") %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "both" = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "both" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "Rapid transition to industrial plantations" = 3,
                     "Slow transition to industrial plantations" = 3),
                   align = "c",
                   strikeout = F) %>%
  pack_rows("CPO medium-run price", 1, 3, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction with \n domestic private ownership", 4, 6, # domestic private ownership  "Interaction with \n domestic private ownership"
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction with \n foreign ownership", 7, 9,
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Interaction with \n # of reachable mills", 10, 12,
            italic = TRUE, bold = TRUE)  %>%
  pack_rows(start_row =  nrow(df), end_row = nrow(df),  latex_gap_space = "0.5em", hline_before = TRUE) %>% 
  column_spec(column = 1,
              width = "12em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(df))),
              width = "4em",
              latex_valign = "b")



rm(df)


#### LEGAL, LUCFIP DYNAMICS, COMMODITY, PRICE DYNAMICS - COMPARE APEs ####

## infrastructure to store results
ape_mat_list <- list()
elm <- 1

comp_ape_mat <- matrix(ncol = 9, nrow = 6, data = NA)
colnames(comp_ape_mat) <- c("Sumatra_a", "Kalimantan_a", "both_a", "Sumatra_i", "Kalimantan_i", "both_i", "Sumatra_sm", "Kalimantan_sm", "both_sm") 
row.names(comp_ape_mat) <- c("Sumatra = Kalimantan", 
                             "industrial = smallholders", 
                             "legal = illegal", 
                             "rapid = slow",
                             "FFB = CPO",
                             "short-run = medium-run")


### COMPARE APEs ACROSS GROUPS
make_APEs_1regr <- function(res_data){
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]]
  
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
  
  # average of the regressor of interest
  reg_bar <- mean(d_clean[,coeff_names[1]]) 
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values)
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  linear_ape_fml <- paste0(coeff_names[1])
  i_t <- 1
  while(i_t<=length(interaction_terms)){
    linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[i_t],"*",int_term_avg[[i_t]])
    i_t <- i_t +1
  }
  ape_fml_roi <- paste0("(",linear_ape_fml,")*fv_bar*",pixel_area_ha)
  
  dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                            vcov. = vcov(reg_res, se = "cluster"), 
                            g. = ape_fml_roi, 
                            rhs = 0)
  
  row.names(dM_ape_roi) <- coeff_names[1]
  dM_ape_roi <- as.matrix(dM_ape_roi)
  dM_ape_roi <- dM_ape_roi[,c(1,2,7)]
  
  dM_ape_list <- list()
  dM_ape_list[[1]] <- dM_ape_roi
  
  ## FORMULA FOR APE OF INTERACTION TERMS 

  # if(length(interaction_terms) >0){
  #   for(i in 1:length(interaction_terms)){
  #     
  #     # there are several parts 
  #     # 1. the coeff of the interaction effect of interest
  #     iei <- paste0(interaction_terms[i], "X", coeff_names[1])
  #     ape_fml1 <- paste0(iei," + ") 
  #     
  #     # 2. the linear expression derivative wrt. the regressor of interest 
  #     ape_fml2 <- paste0("(",linear_ape_fml,")")
  #     
  #     # 3. the derivative of the linear expression wrt. the interaction term 
  #     ape_fml3 <- paste0("*(",
  #                        coeff_names[coeff_names==interaction_terms[i]]," + ",
  #                        iei,"*",reg_bar,")")#
  #     
  #     # paste everything together
  #     ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*",fv_bar,"*",pixel_area_ha)
  #     
  #     dM <- deltaMethod(object = coef(reg_res), 
  #                       vcov. = vcov(reg_res, se = "cluster"), 
  #                       g. = ape_fml_it, 
  #                       rhs = 0)
  #     
  #     row.names(dM) <- interaction_terms[i]
  #     dM <- as.matrix(dM)
  #     dM <- dM[,c(1,2,7)]
  #     dM_ape_list[[i+1]] <- dM
  #     
  #   }
  # }
  mat <- matrix(ncol = 1, 
                nrow = 3+3*length(interaction_terms), 
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

isl_list <- list("Sumatra", "Kalimantan", "both")
size_list <- list("a","i","sm")#"s","m"
dyn_list <- list("rapid", "slow")
# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def))

## Main estimates
for(ISL in isl_list){
  for(SIZE in size_list){
    
    # make the regression
    res_data <- make_base_reg(island = ISL,
                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount_total"), # or can be  lucpf",SIZE,"p_rapidslow_pixelcount"
                              interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                              commo = c("cpo"), #
                              offset = FALSE)
    # make the APE
    ape_mat_list[[elm]] <- make_APEs_1regr(res_data = res_data)
    names(ape_mat_list)[elm] <- paste0(ISL,"_",SIZE)
    rm(res_data)
    elm <- elm + 1
  }
}
# these loops yield 9 estimates, elm should be 10
elm 

## Adding estimations with distinction between legal and illegal plantations
for(ISL in isl_list){#3x
  for(SIZE in size_list){#3x
    for(ILL in ill_status){#2
      if(!(SIZE == "sm" & ISL == "Kalimantan" & ILL == "ill2")){# because in this case there is not enough variation
        # make the regression
        res_data <- make_base_reg(island = ISL,
                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount_total"),
                                  illegal = ILL,
                                  interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                                  commo = c("cpo"), #
                                  offset = FALSE)
        # make the APE
        ape_mat_list[[elm]] <- make_APEs_1regr(res_data = res_data)
        names(ape_mat_list)[elm] <- paste0(ISL,"_",SIZE,"_",ILL)
        rm(res_data)
        elm <- elm + 1
      }
    }
  }
}
# these loops yield 17 estimates, elm should be 27
elm 


## Adding estimations with distinction between rapid and low lucfp 
for(ISL in isl_list){#3x
  SIZE <- "i" # here size is fixed to industrial, as lucfp dynamics are available for these plantations only
  for(DYN in dyn_list){#2x
    
    # make the regression
      res_data <- make_base_reg(island = ISL,
                                outcome_variable = paste0("lucpf",SIZE,"p_",DYN,"_pixelcount"),
                                interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                                commo = c("cpo"), #
                                offset = FALSE)
      # make the APE
      ape_mat_list[[elm]] <- make_APEs_1regr(res_data = res_data)
      names(ape_mat_list)[elm] <- paste0(ISL,"_",SIZE,"_",DYN)
      rm(res_data)
      elm <- elm + 1  
  }
}
# these loops yield 6 estimates, elm should be 33
elm 

length(ape_mat_list)
names(ape_mat_list)
ape_mat_list_saved <- ape_mat_list[[1]]
rm(ape_mat)
ape_mat <- bind_cols(ape_mat_list) %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","SE","p-value"), nrow(ape_mat)/3), "Observations")
# keep ape_mat like this for later comparisons between APEs

# fill the matrix of p values of equality tests BY ROW
for(SIZE in c("a", "i", "sm")){
  comp_ape_mat["Sumatra = Kalimantan",paste0("both_",SIZE)] <- compare_APEs_across_groups(group1 = paste0("Sumatra_",SIZE), 
                                                                                          group2 = paste0("Kalimantan_",SIZE)) 
}
for(ISL in c("both", "Sumatra", "Kalimantan")){
  comp_ape_mat["industrial = smallholders",paste0(ISL,"_a")] <- compare_APEs_across_groups(group1 = paste0(ISL,"_i"), 
                                                                                           group2 = paste0(ISL,"_sm")) 
}
for(ISL in isl_list){
  for(SIZE in size_list){
    if(!(ISL=="Kalimantan" & SIZE == "sm")){
    comp_ape_mat["legal = illegal",paste0(ISL,"_",SIZE)] <- compare_APEs_across_groups(group1 = paste0(ISL,"_",SIZE,"_no_ill2"), 
                                                                                       group2 = paste0(ISL,"_",SIZE,"_ill2")) 
    }
  }
}
for(ISL in isl_list){
  for(DYN in dyn_list){
    comp_ape_mat["rapid = slow",paste0(ISL,"_i")] <- compare_APEs_across_groups(group1 = paste0(ISL,"_i","_rapid"), 
                                                                                group2 = paste0(ISL,"_i","_slow")) 
  }
}



### COMPARE APEs WITHIN GROUPS

# This function, as of now, is written to compare the first and the second coefficients of a fixest estimation result
compare_coeff_within_groups <- function(res_data, m0 = 0, alternative = "two.sided") { 
  reg_res <- res_data[[1]]
  coeff1 <- reg_res$coefficients[1] 
  coeff2 <- reg_res$coefficients[2]
  # # n1 <- ape_mat["Observations","Sumatra industrial"]
  # # n2 <- ape_mat["Observations","Sumatra smallholders"]
  sigma1 <- vcov(reg_res, se = "cluster")[names(coeff1),names(coeff1)]
  sigma2 <- vcov(reg_res, se = "cluster")[names(coeff2),names(coeff2)]
  cov12 <- vcov(reg_res, se = "cluster")[names(coeff1),names(coeff2)]
  
  statistic <- (coeff1 - coeff2 - m0) / sqrt(sigma1 + sigma2 - 2*cov12)
  
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

## Adding estimations with FFB *and* CPO prices  
for(ISL in isl_list){
  for(SIZE in size_list){
    # make the regression
    res_data <- make_base_reg(island = ISL,
                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount_total"), # or can be  lucpf",SIZE,"p_rapidslow_pixelcount"
                              interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                              commo = c("ffb","cpo"), #
                              offset = FALSE)
    # make the APE
    comp_ape_mat["FFB = CPO",paste0(ISL,"_",SIZE)] <- compare_coeff_within_groups(res_data = res_data)
    rm(res_data)
  }
}

## Adding estimations with short-run *and* medium-run CPO prices  
for(ISL in isl_list){
  for(SIZE in size_list){
    # make the regression
    res_data <- make_base_reg(island = ISL,
                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount_total"), # or can be  lucpf",SIZE,"p_rapidslow_pixelcount"
                              interaction_terms = NULL, #c("wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "n_reachable_uml"),#, "wa_pct_own_for_imp" ),
                              dynamics = TRUE,
                              commo = c("cpo"), #
                              offset = FALSE)
    # make the APE
    comp_ape_mat["short-run = medium-run",paste0(ISL,"_",SIZE)] <- compare_coeff_within_groups(res_data = res_data)
    rm(res_data)
  }
}




comp_ape_disp <- comp_ape_mat

comp_ape_disp <- comp_ape_mat %>% formatC(digits = 4, format = "f")
comp_ape_disp[comp_ape_disp=="   NA"] <- "" 
colnames(comp_ape_disp) <- NULL
comp_ape_disp

options(knitr.table.format = "latex")
kable(comp_ape_disp, booktabs = T, align = "r",
      caption = "p-values from equality tests of LUCFP responsiveness to CPO prices") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Ho" = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "Both" = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "Both" = 1,
                     "Sumatra" = 1,
                     "Kalimantan" = 1,
                     "Both" = 1),
                    bold = F,
                    align = "c") %>%  
  add_header_above(c(" " = 1,
                     "All plantations" = 3,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 3
                    ),
                  bold = F,
                  align = "c") %>%
  pack_rows("Equality of APEs \n across groups", 1, 4, 
            italic = TRUE, bold = TRUE)  %>%
  pack_rows("Equality of coefficients \n within groups", 5, 6, 
            italic = TRUE, bold = TRUE)  %>%
  column_spec(column = 1,
              width = "12em",
              latex_valign = "b") %>% 
  column_spec(column = c(2:(ncol(comp_ape_disp))),
              width = "5em",
              latex_valign = "b")

#### TABLE 5 - COMPARISONS OF APEs ACROSS MAIN GROUPS #### 

#### TABLE 5 - COMPARISONS OF APEs ACROSS OTHER GROUPS #### 




####-------------------------------------------------------------BELOW IS PREVIOUS VERSION OF DISPLAYING RESULTS--------------------------------------------------------------------------####
# we do not delete it now to compare the CA results with the CR in the leaflet currently

##### TABLES OF ELASTICITIES ##### 
# pass arguments from make_base_reg in lapply to get alternative specifications

ov_list <- list("lucpfip_pixelcount_total", "lucpfsmp_pixelcount_total")
isl_list <- list("Sumatra", "Kalimantan")
### WITH DISTINCTIONS SHORT/medium TERM & FFB/CPO PRICE SIGNALS 


## Elasticities 

# Industrial
res_list_dst_ind <- lapply(isl_list, make_base_reg,
                           outcome_variable = "lucpfip_pixelcount_total",
                           dynamics = TRUE,
                           commo = c("ffb","cpo"), #
                           offset = FALSE) 

# Smallholders
res_list_dst_sm <- lapply(isl_list, make_base_reg,
                          outcome_variable = "lucpfsmp_pixelcount_total",
                          dynamics = TRUE,
                          commo = c("ffb", "cpo"),#
                          offset = FALSE) 

res_list_dst <- append(res_list_dst_ind, res_list_dst_sm)
rm(res_list_dst_ind, res_list_dst_sm)
# preview in R 
etable(res_list_dst, 
       se = "cluster", 
       #coefstat = "confint",
       subtitles = c("Sumatra", "Kalimantan","Sumatra", "Kalimantan"))


# In LateX
# title

#table_title_dyn <- paste0("LUCFP semi-elasticities to short and medium-run price signals") 

table_title_dyn <- paste0("LUCFP elasticities to short and medium-run price signals") 


#table_title_dyn <- paste0("LUCFP semi-elasticities to short and medium-run y-o-y growth rates of price signals") 


# LateX table
etable(res_list_dst, 
       #cluster = oneway_cluster,
       se = "cluster",
       tex = TRUE,
       # file = table_file, 
       # replace = TRUE,
       title = table_title_dyn,
       subtitles = c("Sumatra", "Kalimantan", "Sumatra",  "Kalimantan"),
       family = TRUE,
       drop = c("own", "reachable"),
       #coefstat = "confint",
       sdBelow = F,
       yesNo = "X",
       fitstat = c("sq.cor"),
       dict = TRUE, 
       powerBelow = -7)

rm(res_list_dst)
# ILLEGAL vs LEGAL LUCFP
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def))
ov_list <- list("lucpfip_pixelcount_total", "lucpfsmp_pixelcount_total")
isl_list <- list("Sumatra", "Kalimantan")
res_list_dst_ill <- list()
elm <- 1
for(OV in ov_list){
  for(ISL in isl_list){
    for(ILL in ill_status){
      res_list_dst_ill[[elm]] <- make_base_reg(island = ISL, 
                                               outcome_variable = OV,
                                               illegal = ILL,
                                               dynamics = TRUE,
                                               commo = c("ffb","cpo"), #
                                               offset = FALSE)
      elm <- elm + 1
    }
  }
}


# preview in R 
etable(res_list_dst_ill, 
       se = "cluster", 
       subtitles = c("Sumatra, legal", "Sumatra, illegal", 
                     "Kalimantan, legal", "Kalimantan, illegal", 
                     "Sumatra, legal", "Sumatra, illegal", 
                     "Kalimantan, legal"))#, "Kalimantan, illegal" this one has dependent variable constant (for both ill def). 

# In LateX
# title

#table_title_dyn <- paste0("Legal and illegal LUCFP semi-elasticities to short and medium-run price signals") 

table_title_dyn <- paste0("Legal and illegal LUCFP elasticities to short and medium-run price signals") 

#table_title_dyn <- paste0("Legal and illegal LUCFP semi-elasticities to short and medium-run y-o-y growth rates of price signals") 


# LateX table
etable(res_list_dst_ill, 
       #cluster = oneway_cluster,
       se = "cluster",
       tex = TRUE,
       # file = table_file, 
       # replace = TRUE,
       title = table_title_dyn,
       subtitles = c("Sumatra, legal", "Sumatra, illegal", 
                     "Kalimantan, legal", "Kalimantan, illegal", 
                     "Sumatra, legal", "Sumatra, illegal", 
                     "Kalimantan, legal"),#, "Kalimantan, illegal" this one has dependent variable constant.        
       drop = c("own", "reachable"),
       #coefstat = "confint",
       sdBelow = TRUE,
       yesNo = "X",
       fitstat = c("sq.cor"),
       dict = TRUE, 
       powerBelow = -7)

### INTERACTIONS 

## Elasticities 

ov_list <- list("lucpfip_pixelcount_total", "lucpfsmp_pixelcount_total")
isl_list <- list("Sumatra", "Kalimantan")
# Industrial
res_int_list_dst_ind <- lapply(isl_list, make_base_reg,
                               outcome_variable = "lucpfip_pixelcount_total",
                               dynamics = FALSE,
                               interaction_terms = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"),#n_reachable_uml
                               commo = c("cpo"), #"ffb",
                               offset = FALSE) 

# Smallholders
res_int_list_dst_sm <- lapply(isl_list, make_base_reg,
                              outcome_variable = "lucpfsmp_pixelcount_total",
                              dynamics = FALSE,
                              interaction_terms = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"),#n_reachable_uml
                              commo = c("cpo"),#"ffb",
                              offset = FALSE) 

res_int_list_dst <- append(res_int_list_dst_ind, res_int_list_dst_sm)

# preview in R 
etable(res_int_list_dst, 
       se = "cluster", 
       subtitles = c("Sumatra", "Kalimantan","Sumatra", "Kalimantan"))

# In LateX
# title

#table_title_dyn <- paste0("LUCPFP semi-elasticities to short and medium term price signals") 

table_title_dyn <- paste0("LUCPFP elasticities to short and medium term price signals") 

#table_title_dyn <- paste0("LUCPFP semi-elasticities to short and medium term y-o-y growth rates of price signals") 


# LateX table
etable(res_list_dst, 
       #cluster = oneway_cluster,
       se = "cluster",
       tex = TRUE,
       # file = table_file, 
       # replace = TRUE,
       title = table_title_dyn,
       subtitles = c("Sumatra", "Kalimantan", "Sumatra",  "Kalimantan"),
       family = TRUE,
       drop = c("reachable"), # "own", "reachable"
       sdBelow = FALSE,
       yesNoFixef = "X",
       fitstat = c("sq.cor"),
       dict = TRUE, 
       powerBelow = -7)


##### TABLES OF MARGINAL EFFECT  ##### 

ov_list <- list("lucpfip_pixelcount_total", "lucpfsmp_pixelcount_total")
isl_list <- list("Sumatra", "Kalimantan")
commo_list <- list(c("ffb"), "cpo")#, "cpo"
res_list <- list()
elm <- 1
for(OV in ov_list){
  for(ISL in isl_list){
    for(COMMO in commo_list){
      res_list[[elm]] <- make_base_reg(island = ISL, 
                                       outcome_variable = OV,
                                       dynamics = FALSE,
                                       yoyg = FALSE,
                                       illegal = "no_ill2",
                                       commo = COMMO,
                                       remaining_forest = FALSE,
                                       offset = FALSE)
      elm <- elm + 1
    }
  }
}

# preview in R 
etable(res_list, 
       se = "cluster", 
       coefstat = "confint")

# Latex table 
# title
#table_title <- paste0("LUCFP semi-elasticities to medium-run y-o-y growth rates of price signals") 

table_title <- paste0("LUCFP elasticities to medium-run price signals")
# LateX table
etable(res_list,
       #cluster = oneway_cluster,
       se = "cluster",
       tex = TRUE,
       # file = table_file,
       # replace = TRUE,
       title = table_title,
       subtitles = c("Sumatra", "Sumatra", "Kalimantan", "Kalimantan",
                     "Sumatra", "Sumatra", "Kalimantan", "Kalimantan"),
       family = TRUE,
       drop = c("own", "reachable"),
       #coefstat = "confint",
       sdBelow = TRUE,
       yesNo = "X",
       fitstat = c("sq.cor"),
       dict = TRUE,
       powerBelow = -7)

res <- res_list[[2]]
# keep only groups with significant and positive coefficients for policy insights
cond <- lapply(res_list, FUN = function(res){abs(summary(res, se = "cluster")$coeftable[1,"t value"])>1.645 & res$coefficients[1] > 0})
res_list_signif <- res_list[unlist(cond)]


# helper function that transforms the list of results into a data frame of 
make_price_change <- function(res){
  
  coeff <- res$coefficients[[1]]
  price_change <- -(1/coeff) 
  col <- as.matrix(c(0.1,0.5,1)*price_change)
  colnames(col) <- "colname"
  return(col)
}
df <- bind_cols(lapply(res_list_signif, FUN = make_price_change)) %>% as.matrix()

# prepare df for kable
row.names(df) <- paste0(c(10, 50, 100),"% LUCFP reduction")
df <- df*100
df <- df %>% round(digits = 0)
df <- df %>% as.data.frame()
df <- apply(df, c(1,2), paste0,"%")
colnames(df) <- NULL


options(knitr.table.format = "latex") 
kable(df, booktabs = T, align = "r", 
      caption = "Price distortions against deforestation-based products to achieve different reductions in LUCFP ") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "FFB" = 1, "CPO" = 1,
                     #"FFB" = 1, #"CPO" = 1,
                     #"FFB" = 1, "CPO" = 1,
                     #"FFB" = 1, 
                     "CPO" = 1
  ), 
  bold = F, 
  align = "c") %>% 
  add_header_above(c(" " = 1, 
                     "Sumatra" = 3
                     #"Kalimantan" = 1,
                     #"Sumatra" = 2,
                     #"Kalimantan" = 1
  ), 
  align = "c", 
  strikeout = F) %>% 
  
  add_header_above(c(" " = 1, 
                     "Industrial plantations" = 2, 
                     "Smallholder plantations" = 1), 
                   bold = T, 
                   align = "c") %>% 
  column_spec(column = c(2:(ncol(df)+1)),
              width = "5em", 
              latex_valign = "b") 
# 
# footnote(general = c("Sample mills are the 470 IBS palm oil mills we matched with the UML database."),
#          threeparttable = TRUE, 
#          escape = TRUE) 



pixel_area_ha <- (27.8*27.6)/(1e4)
# helper function that transforms the list of results into a data frame of partial effects. 
make_PE <- function(res, variation){
  # res is an element of res_list_main, so it's a fixest object. 
  
  # Compute the appropriate fitted values first. 
  
  coeff <- res$coefficients[[1]]
  avg_fv <- mean(res$fitted.values)
  ape <- variation*coeff*avg_fv
  
  if(grepl("pixelcount",as.character(res$fml)[2])){
    ape <- ape*pixel_area_ha
  }
  
  
} 

apes <- sapply(res_list_main, make_PE, variation = 1)


# save fitted values and linear predictors (which include FEs)
d_clean$fitted.values <- fe_reg$fitted.values
d_clean$linear.predictors <- fe_reg$linear.predictors

## AVERAGE FITTED VALUE 
distr_avg_fv <- ddply(d_clean, "district", summarise, 
                      distr_avg_fv = mean(fitted.values))

distr_lvl <- merge(distr_lvl, distr_avg_fv, by = "district", all = TRUE)
# replace NAs with island average
distr_lvl$distr_avg_fv <- replace(x = distr_lvl$distr_avg_fv, 
                                  list = is.na(distr_lvl$distr_avg_fv), 
                                  values = mean(distr_avg_fv$distr_avg_fv))
distr_lvl$distr_avg_fv_ha <- distr_lvl$distr_avg_fv*(27.8*27.6)/(1e4)




### DEMAND FUNCTIONS 

ffb_demand_avg <- function(tax, fitted_value_avg, aggr_factor, coeff_ffb){
  fitted_value_avg*aggr_factor/exp(coeff_ffb*tax)
}
ffb_demand_avg(tax = 0, 
               fitted_value_avg = elmts[["Sumatra"]]["fitted_value_avg"],
               aggr_factor = 1,
               coeff_ffb = elmts[["Sumatra"]]["coeff_ffb"])

cpo_demand_avg <- function(tax, fitted_value_avg, aggr_factor, coeff_cpo){
  fitted_value_avg*aggr_factor/exp(coeff_cpo*tax)
}
cpo_demand_avg(tax = 0, 
               fitted_value_avg = elmts[["Sumatra"]]["fitted_value_avg"],
               aggr_factor = 1,
               coeff_cpo = elmts[["Sumatra"]]["coeff_cpo"])








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
                               alt_ca = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                               outcome_variable = "lucpfip_pixelcount_total", # LHS. One of "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total"p
                               commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo")), commodities the price signals of which should be included in the RHS
                               log_prices = TRUE, # Logical, should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
                               x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                               dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                               yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                               short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                               imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                               distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                               fe = "parcel_id + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                               remaining_forest = FALSE, # Logical. If TRUE, the remaining forest is added as a control
                               offset = FALSE, # Logical. Should the log of the remaining forest be added as an offset.  
                               lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                               controls = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml"), # character vectors of names of control variables (don't specify lags in their names)
                               spatial_control = FALSE, # logical, if TRUE, adds ~30min computation. Should the average of neighbors' outcome variable be added in the RHS. 
                               pya_ov = FALSE, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                               weights = FALSE, # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
                               # additional argument to this function: 
                               variable = "cpo",
                               cluster = "cluster"
){
  
  ### Get coefficient and SE for a given specification passed to make_base_reg.
  
  # call make_base_reg
  reg_res <- make_base_reg(island = island,
                           outcome_variable = outcome_variable,
                           commo = commo, 
                           log_prices = log_prices,
                           x_pya = x_pya, 
                           dynamics = dynamics, 
                           yoyg = yoyg, 
                           short_run = short_run, 
                           imp = imp, 
                           distribution = distribution, 
                           fe = fe, 
                           remaining_forest = remaining_forest, 
                           offset = offset,
                           lag_or_not = lag_or_not, 
                           controls = controls,
                           spatial_control = spatial_control, 
                           pya_ov = pya_ov, 
                           weights = weights)
  
  # compute SE the way specified by cluster
  reg_summary <- summary(reg_res, se = cluster)
  # get coeff and SE
  coeff_and_se <- reg_summary$coeftable[match(variable, commo), 1:2]
  
  ### make indicator variables that will be used to label specifications. 
  ind_var <- data.frame(# outcome variable
    "offset" = FALSE,
    "larger_forest_def" = FALSE,
    # sample
    "CA_2h" = FALSE,
    "CA_4h" = FALSE,
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
    "remaining_forest" = FALSE,
    "ngb_ov_lag4" = FALSE,
    # weights
    "weights" = FALSE,
    # standard errors
    "oneway_cluster" = FALSE, 
    "twoway_cluster" = FALSE
  )
  
  # change the OV indicator variable to TRUE 
  if(offset){ind_var[,"offset"] <- TRUE}
  if(grepl("0th", outcome_variable)){ind_var[,"larger_forest_def"] <- TRUE}
  
  # sample  
  # change the CA indicator variable to TRUE for the corresponding CA
  # Catchment radius
  if(island == "Sumatra"){
    travel_time <- 2
  }else{
    travel_time <- 4}
  
  if(alt_ca){  
    if(island == "Sumatra"){
      travel_time <- 4
    }else{
      travel_time <- 2
    }
  }
  ind_var[,grepl(paste0("CA_", travel_time,"h"), colnames(ind_var))] <- TRUE
  
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
  if(fe == "parcel_id"){ind_var[,"unit_fe"] <- TRUE}
  if(fe == "parcel_id + year"){ind_var[,"tw_fe"] <- TRUE}
  if(fe == "parcel_id + province_year"){ind_var[,"unit_provyear_fe"] <- TRUE}
  if(fe == "parcel_id + district_year"){ind_var[,"unit_distryear_fe"] <- TRUE}
  
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
  if(pya_ov){ind_var[,"remaining_forest"] <- TRUE}
  # pya outcome spatial
  if(spatial_control){ind_var[,"ngb_ov_lag4"] <- TRUE}
  
  # weights indicator variable
  if(weights){ind_var[,"weights"] <- TRUE}
  if(distribution == "negbin"){ind_var[,"weights"] <- FALSE} # turn it back to FALSE if negbin distribution
  
  # clustering
  if(cluster == "cluster"){ind_var[,"oneway_cluster"] <- TRUE}
  if(cluster == "twoway"){ind_var[,"twoway_cluster"] <- TRUE}
  
  ### Bind together coeff, SE, and specification labels
  spec_df <- cbind(coeff_and_se, ind_var)
  return(spec_df)
  
}

### COMPUTE THE DATASETS FOR EACH ISLAND, AND COMMODITY OF INTEREST
for(ISL in c("Sumatra", "Kalimantan")){
  for(OV in c("lucpfip_pixelcount_total", "lucpfsmp_pixelcount_total")){
    VAR <- "cpo"
    
    #if(!file.exists(file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",VAR,"_",OV)))){
    
    reg_stats_indvar_list <- list()
    i <- 1
    ### Add to the list the specifications that are to be compared in the chart. 
    # in make_spec_chart_df, specify only arguments over which it's looping. Other arguments are set to default. 
    
    # First loops over critical parameters
    for(IMP in c(1, 2)){
      #for(DISTR in c("quasipoisson", "negbin")){#
      for(XPYA in c(2, 3, 4)){
        for(LAG in c("_lag1", "")){
          for(WGH in c(TRUE, FALSE)){
            for(OFF in c(TRUE, FALSE)){
              reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                               outcome_variable = OV,
                                                               offset = OFF,
                                                               imp = IMP,
                                                               x_pya = XPYA,
                                                               lag_or_not = LAG,
                                                               weights = WGH,
                                                               variable = VAR)
              i <- i+1
            }
          }
        }
      }
    }
    #}
    
    # Then add particular departures from the preferred specification (which arguments are already set by default)
    
    ## For an alternative catchment radius (with and without weights)
    for(WGH in c(TRUE, FALSE)){
      for(OFF in c(TRUE, FALSE)){
        reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                         outcome_variable = OV,
                                                         offset = OFF,
                                                         alt_ca = TRUE,
                                                         weights = WGH,
                                                         variable = VAR)
        i <- i+1
      }
    } 
    ## For an alternative outcome variable
    if(OV == "lucpfip_pixelcount_total"){ 
      for(OFF in c(TRUE, FALSE)){
        reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                         outcome_variable = "lucfip_pixelcount_30th",
                                                         offset = OFF,
                                                         variable = VAR)
        i <- i+1
      }
    }
    # if(OV == "lucpfsmp_pixelcount_total"){
    # for(OFF in c(TRUE, FALSE)){
    #   reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
    #                                                    outcome_variable = "lucfsmp_pixelcount_30th",
    #                                                    offset = OFF,
    #                                                    variable = VAR)
    #   i <- i+1
    # }
    # }
    
    ## For alternative distributionnal assumptions
    for(DISTR in c("poisson", "negbin")){
      reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                       outcome_variable = OV,
                                                       distribution = DISTR,
                                                       variable = VAR)
      i <- i+1
    }
    
    ## For alternative fixed effects
    for(FE in c("parcel_id", "parcel_id + year", "parcel_id + province_year")){
      reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                       outcome_variable = OV,
                                                       fe = FE,
                                                       variable = VAR)
      i <- i+1
    }
    
    ## For an alternative standard error computation (two-way clustering)
    reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = OV,
                                                     cluster = "twoway",
                                                     variable = VAR)
    i <- i+1
    
    ## For alternative control sets
    # Without ownership control
    reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = OV,
                                                     control = c("n_reachable_uml"),
                                                     variable = VAR)
    i <- i+1
    
    # Without remoteness/competition control
    reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = OV,
                                                     control = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"),
                                                     variable = VAR)
    i <- i+1
    
    # With FFB price signal control
    reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = OV,
                                                     commo = c("ffb", "cpo"),
                                                     variable = VAR)
    i <- i+1
    
    
    # With percentage exported control
    reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = OV,
                                                     control = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"),
                                                     variable = VAR)
    i <- i+1
    
    
    # With past remaining forest control 
    reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                     outcome_variable = OV,
                                                     remaining_forest = TRUE,
                                                     variable = VAR)
    i <- i+1
    
    
    # With average of neighbors' 4-year-lagged outcome   
    long_to_compute <- file.path(paste0("temp_data/reg_results/spec_chart_df_spatial_",ISL,"_",VAR))
    if(file.exists(long_to_compute)){
      reg_stats_indvar_list[[i]] <- readRDS(long_to_compute)
    }else{
      reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                       outcome_variable = OV,
                                                       spatial_control = TRUE,
                                                       variable = VAR)
    }
    i <- i+1
    
    
    
    
    # convert to dataframe to be able to chart
    reg_stats_indvar <- bind_rows(reg_stats_indvar_list)
    
    if(sum(duplicated(reg_stats_indvar))==0 & nrow(reg_stats_indvar)+1 == i){
      # save it 
      if(length(ISL) ==1){
        saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",VAR,"_",OV)))
      }else{
        saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_all_",VAR,"_",OV)))
      }
    }else{print(paste0("SOMETHING WENT WRONG in spec_chart_df_",ISL,"_",VAR,"_",OV))}
    
    #}else{reg_stats_indvar <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",VAR,"_",OV)))}
  }
}

### PLOTTING 
### GIVE HERE THE ISLAND, THE OUTCOME AND THE COMMODITY FOR WHICH YOU WANT THE SPEC CHART TO BE PLOTTED
ISL <- "Kalimantan"
OV <-   "lucpfsmp_pixelcount_total"
VAR <- "ffb"
# make labels for the chart
schart_labels <- list("Dependent variable:" = c("Relative to remaining forest",
                                                "Larger forest definition"),
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
                                      "Remaining forest", 
                                      "Neighbors' outcomes, 4-year lagged"),
                      "Weights" = "", 
                      "Standard errors:" = c("Grid cell cluster", 
                                             "Two-way cluster")
) 

# remove rows where SEs could not be computed
problems <- reg_stats_indvar[is.na(reg_stats_indvar$`Std. Error`),]
reg_stats_indvar <- reg_stats_indvar[!is.na(reg_stats_indvar$`Std. Error`),]
rownames(reg_stats_indvar) <- seq(1, nrow(reg_stats_indvar))

# find position of model to highlight in original data frame
a <- reg_stats_indvar
a <- a[a$imp1 &
         #a$CA_2h &  
         a$quasipoisson &
         a$pya_4 &
         a$lag_or_not &
         a$unit_distryear_fe &
         a$two_commo == FALSE & 
         a$control_own & 
         a$n_reachable_uml_control  &
         a$prex_cpo_control == FALSE &
         a$remaining_forest == FALSE & 
         a$weights == FALSE & 
         a$oneway_cluster, ] 

if(ISL == "Sumatra"){
  model_idx <- a[a$CA_2h, ] %>% rownames() 
}else{  
  model_idx <- a[a$CA_4h, ] %>% rownames() 
}

# this is our baseline model(s)
a[model_idx,]

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
       ylab = "Semi-elasticity estimates",
       lwd.est = 5.8,
       lwd.symbol = 1
       #pch.est=20
)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#### MAIN - APEs OF INTERACTION TERMS 
# make_APE_interactions <- function(res_data, int_term){
#   reg_res <- res_data[[1]]
#   d_clean <- res_data[[2]]
#   
#   # data POINTS that deltaMethod will grab 
#   
#   # averages of the interaction terms -the controls that we interacted with the regressor of interest)  
#   int_terms <- !(grepl(pattern = "X", names(coef(reg_res))))
#   int_terms[1] <- FALSE
#   
#   int_terms <- names(coef(reg_res))[int_terms]
#   int_term_avg <- list()
#   for(i in 1:length(int_terms)){
#     int_term_avg[[i]] <- mean(d_clean[,int_terms[i]])
#   }
#   
#   # average of the regressor of interest
#   reg_bar <- mean(d_clean[,names(coef(reg_res))[1]])
#   
#   # average fitted values
#   fv_bar <- mean(reg_res$fitted.values)
#   
#   # names of coefficients of interactions (all those which names include the name of the first regressor)
#   varn <- coef(reg_res) %>% names()
#   varn <- varn[grepl(varn[1], varn)]
#   
#   # Formula for APE of interaction term
#   # there are several parts 
#   # 1. the coeff of the interaction effect of interest
#   #  int_term <- "n_reachable_uml_lag1"
#   iei <- varn[grepl(int_term, varn)]
#   ape_fml1 <- paste0(iei," + ") 
#   
#   # 2. the linear expression derivative wrt. regressor of interest 
#   ape_fml2 <- paste0("(",varn[1])
#   i_t <- 1
#   while(i_t<=length(int_terms)){
#     ape_fml2 <- paste0(ape_fml2," + ", varn[i_t+1],"*",int_term_avg[[i_t]])
#     i_t <- i_t +1
#   }
#   ape_fml2 <- paste0(ape_fml2,")")
#   
#   # 3. the derivative of the linear expression wrt. the interaction term 
#   ape_fml3 <- paste0("*(",
#                      names(coef(reg_res))[names(coef(reg_res))==int_term]," + ",
#                      iei,"*",reg_bar,")")#
#   
#   # paste everything together
#   ape_fml <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*",fv_bar,"*",pixel_area_ha)
#   
#   dM <- deltaMethod(object = coef(reg_res), 
#                     vcov. = vcov(reg_res, se = "cluster"), 
#                     g. = ape_fml, 
#                     rhs = 0)
#   
#   row.names(dM) <- paste0("APE")
#   
#   rm(reg_res, d_clean, int_terms, int_term_avg, varn, ape_fml)
#   return(dM)
# }
# 
# ### INTERACTION WITH DOMESTIC PRIVATE OWNERSHIP
# 
# rm(df)
# df <- bind_rows(lapply(res_data_list, FUN = make_APE_interactions, int_term = "wa_pct_own_nat_priv_imp_lag1")) %>% as.matrix()
# # prepare df for kable
# # row.names(df) <- paste0("APE - ",names(res_data_list))
# df %>% class()
# df <- t(df)
# df <- df[c(1,2,7),]
# row.names(df)[3] <- "p-value"
# df <- df %>% round(digits = 2)
# df <- df %>% as.data.frame()
# colnames(df) <- NULL
# 
# 
# options(knitr.table.format = "latex")
# kable(df, booktabs = T, align = "r",
#       caption = "Average partial effects of interactions between medium-run CPO price signal and ownership shift from public to domestic private, on LUCFP (ha)") %>%
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c(" " = 1,
#                      " " = 1, 
#                      "Sumatra" = 1,
#                      "Kalimantan" = 1,
#                      "both" = 1,
#                      "Sumatra" = 1,
#                      "Kalimantan" = 1,
#                      "both"),
#                    bold = F,
#                    align = "c") %>%
#   add_header_above(c(" " = 1,
#                      "All" = 1,
#                      "Industrial plantations" = 3,
#                      "Smallholder plantations" = 3),
#                    align = "c",
#                    strikeout = F) %>%
#   column_spec(column = c(2:(ncol(df))),
#               width = "5em",
#               latex_valign = "b")
# 
# ###  INTERACTION WITH FOREIGN OWNERSHIP
# rm(df)
# df <- bind_rows(lapply(res_data_list, FUN = make_APE_interactions, int_term = "wa_pct_own_for_imp_lag1")) %>% as.matrix()
# # prepare df for kable
# # row.names(df) <- paste0("APE - ",names(res_data_list))
# df %>% class()
# df <- t(df)
# df <- df[c(1,2,7),]
# row.names(df)[3] <- "p-value"
# df <- df %>% round(digits = 2)
# df <- df %>% as.data.frame()
# colnames(df) <- NULL
# 
# 
# options(knitr.table.format = "latex")
# kable(df, booktabs = T, align = "r",
#       caption = "Average partial effects of interactions between medium-run CPO price signal and ownership shift from public to foreign, on LUCFP (ha)") %>%
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c(" " = 1,
#                      " " = 1, 
#                      "Sumatra" = 1,
#                      "Kalimantan" = 1,
#                      "both" = 1,
#                      "Sumatra" = 1,
#                      "Kalimantan" = 1,
#                      "both"),
#                    bold = F,
#                    align = "c") %>%
#   add_header_above(c(" " = 1,
#                      "All" = 1,
#                      "Industrial plantations" = 3,
#                      "Smallholder plantations" = 3),
#                    align = "c",
#                    strikeout = F) %>%
#   column_spec(column = c(2:(ncol(df))),
#               width = "5em",
#               latex_valign = "b")
# 
# 
# ###  INTERACTION WITH NUMBER OF REACHABLE UML
# rm(df)
# df <- bind_rows(lapply(res_data_list, FUN = make_APE_interactions, int_term = "n_reachable_uml_lag1")) %>% as.matrix()
# # prepare df for kable
# # row.names(df) <- paste0("APE - ",names(res_data_list))
# df %>% class()
# df <- t(df)
# df <- df[c(1,2,7),]
# row.names(df)[3] <- "p-value"
# df <- df %>% round(digits = 2)
# df <- df %>% as.data.frame()
# colnames(df) <- NULL
# 
# 
# options(knitr.table.format = "latex")
# kable(df, booktabs = T, align = "r",
#       caption = "Average partial effects of interactions between medium-run CPO price signal and number of reachable mills, on LUCFP (ha)") %>%
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c(" " = 1,
#                      " " = 1, 
#                      "Sumatra" = 1,
#                      "Kalimantan" = 1,
#                      "both" = 1,
#                      "Sumatra" = 1,
#                      "Kalimantan" = 1,
#                      "both"),
#                    bold = F,
#                    align = "c") %>%
#   add_header_above(c(" " = 1,
#                      "All" = 1,
#                      "Industrial plantations" = 3,
#                      "Smallholder plantations" = 3),
#                    align = "c",
#                    strikeout = F) %>%
#   column_spec(column = c(2:(ncol(df))),
#               width = "5em",
#               latex_valign = "b")
# 
# 
# 
# 
# 
# rm(res_data, res_data_list, res_data_list_all, res_data_list_ind, res_data_list_sm, 
#    res_list, res_list_all, res_list_ind, res_list_sm)
# 


