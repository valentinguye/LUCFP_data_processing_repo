
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


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### REGRESSION FUNCTION ##### 

# Prefered specifications are passed as default arguments. 
# Currently, the returned object is a fixest object, of *ONE* regression.  


make_base_reg <- function(island,
                            alt_cr = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                            outcome_variable = "lucpfip_pixelcount_total", # LHS. One of "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total"p
                            commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo")), commodities the price signals of which should be included in the RHS
                            x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                            dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                            yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                            short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                            imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                            distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                            fe = "parcel_id + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                            lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                            controls = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml"), # character vectors of names of control variables (don't specify lags in their names)
                            spatial_control = FALSE, # logical, if TRUE, adds ~30min computation. Should the average of neighbors' outcome variable be added in the RHS. 
                            pya_ov = FALSE, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                            weights = FALSE # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
){
  
  ### SPECIFICATIONS  
  
  ## variable of interests (called regressors here)
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
  
  # add pya outcome variable 
  if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}
  
  # lag controls or not
  if(lag_or_not=="_lag1"){controls <- paste0(controls,lag_or_not)}
  
  # add spatial control or not
  if(spatial_control){controls <- c(controls, "ngb_ov_lag4")}
  
  
  ### DATA FOR REGRESSIONS
  
  # Catchment radius
  if(island == "Sumatra"){
    catchment_radius <- 3e4
  }else{
    catchment_radius <- 5e4}
  
  if(alt_cr){  
    if(island == "Sumatra"){
    catchment_radius <- 5e4
    }else{
      catchment_radius <- 3e4
    }
  }
  
  # Read in full final data set, output of merge_lhs_rhs_parcels.R
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
  d_norspo <- d[d$rspo_cert==FALSE,]
  
  # - have no NA on any of the variables used (otherwise they get removed by {fixest})
  used_vars <- c("parcel_id", "year", "lat", "lon", "district", "province", "island", "district_year", "province_year",
                 "n_reachable_ibsuml_lag1", "sample_coverage_lag1", 
                 outcome_variable, regressors, controls)
  
  filter_vec <- base::rowSums(!is.na(d_norspo[,used_vars]))
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d_norspo[filter_vec, used_vars]
  if(anyNA(d_nona)){stop()}
  
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
  
  # note that obs2remove has to the last filtering operation on data, otherwise some units (along fe dimension)
  # may become "only-zero-outcome" units after other data removal.  
  

  
  ### REGRESSIONS
  
  # Formula
  fe_model <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                paste0(regressors, collapse = "+"),
                                " + ",
                                paste0(controls, collapse = "+"),
                                " | ",
                                fe))

  # run regression for each fixed-effect model 
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    if(weights == TRUE){
      var_weights <- d_clean$sample_coverage_lag1/100 
      reg_res <- fixest::feglm(fe_model,
                              data = d_clean, 
                              family = distribution, 
                              notes = TRUE, 
                              weights = var_weights)
      
    }else{
      reg_res <- fixest::feglm(fe_model,
                              data = d_clean, 
                              family = distribution, 
                              notes = TRUE)
    }
  }else{ # no weights allowed in negative binomial
    reg_res <- fixest::fenegbin(fe_model,
                               data = d_clean, 
                               family = distribution, 
                               notes = TRUE)
  }
  
  rm(d)
  return(reg_res)
}




##### REGRESSION TABLES ##### 
island_list <- list("Sumatra", "Kalimantan")#c("Sumatra", "Kalimantan", "Papua"), 
# pass arguments from make_base_reg in lapply to get alternative specifications
# WITHOUT dynamics and CPO only
result_list_nodyn <- lapply(island_list, make_base_reg) %>% unlist(recursive = FALSE)
# WITH dynamics and both FFB and CPO 
result_list_dyn <- lapply(island_list, make_base_reg,                           
                          dynamics = TRUE,
                          commo = c("ffb", "cpo")) %>% unlist(recursive = FALSE)
# SHAPE TABLES WRT. SPECIFICATIONS PASSED ABOVE
# NOT DONE YET ### ### ### ### ### ### ### ### 
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
                                      "neighbor outcome, 4-year lagged"),
                      "Weights" = "", 
                      "Standard errors:" = c("Grid cell cluster", 
                                             "Two-way cluster")
) 



## We now produce the "data" argument for the function schart, i.e. we run one regression, and bind in a dataframe 
# the coefficients, the SE, and indicator variables for the correponding specifications. 
make_spec_chart_df <- function(island,
                              alt_cr = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                              outcome_variable = "lucpfip_pixelcount_total", # LHS. One of "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total"p
                              commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo")), commodities the price signals of which should be included in the RHS
                              x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                              dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                              yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                              short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                              imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                              distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                              fe = "parcel_id + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
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
                            x_pya = x_pya, 
                            dynamics = dynamics, 
                            yoyg = yoyg, 
                            short_run = short_run, 
                            imp = imp, 
                            distribution = distribution, 
                            fe = fe, 
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
    "ngb_ov_lag4" = FALSE,
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
  if(pya_ov){ind_var[,"pya_outcome_control"] <- TRUE}
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


### GIVE HERE THE ISLAND AND COMMODITY FOR WHICH YOU WANT THE SPEC CHART TO BE EDITED
ISL <- "Kalimantan"
VAR <- "cpo"

reg_stats_indvar_list <- list()

### Add to the list the specifications that are to be compared in the chart. 
# in make_spec_chart_df, specify only arguments over which it's looping. Other arguments are set to default. 

# First loops over critical parameters
for(IMP in c(1, 2)){
  #for(DISTR in c("quasipoisson", "negbin")){#
      for(XPYA in c(2, 3, 4)){
        for(LAG in c("_lag1", "")){
            for(WGH in c(TRUE, FALSE)){
            reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                             outcome_variable = OV,
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
#}

# Then add particular departures from the preferred specification (which arguments are already set by default)

# For an alternative catchment radius (with and without weights)
for(WGH in c(TRUE, FALSE)){
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   alt_cr = TRUE,
                                                   weights = WGH,
                                                   variable = VAR)
  i <- i+1
}

# For an alternative outcome variable
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 outcome_variable = "lucfip_pixelcount_30th",
                                                 variable = VAR)
i <- i+1

# For alternative distributionnal assumptions
for(DISTR in c("poisson", "negbin")){
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   distribution = DISTR,
                                                   variable = VAR)
  i <- i+1
}

# For alternative fixed effects
for(FE in c("parcel_id", "parcel_id + year", "parcel_id + province_year")){
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   fe = FE,
                                                   variable = VAR)
  i <- i+1
}

## For alternaive control sets
# Without ownership control
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 control = c("n_reachable_uml"),
                                                 variable = VAR)
i <- i+1

# Without remoteness/competition control
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 control = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"),
                                                 variable = VAR)
i <- i+1

# With FFB price signal control
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 commo = c("ffb", "cpo"),
                                                 variable = VAR)
i <- i+1


# With percentage exported control
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 control = c("wa_pct_own_loc_gov_imp","wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"),
                                                 variable = VAR)
i <- i+1


# With past year average outcome 
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 pya_ov = TRUE,
                                                 variable = VAR)
i <- i+1

# With average of neighbors' 4-year-lagged outcome   
long_to_compute <- file.path(paste0("temp_data/reg_results/spec_chart_df_spatial_",ISL,"_",VAR))
if(file.exists(long_to_compute)){
  reg_stats_indvar_list[[i]] <- readRDS(long_to_compute)
}else{
  reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                   spatial_control = TRUE,
                                                   variable = VAR)
}
i <- i+1


## For an alternative standard error computation (two-way clustering)
reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
                                                 cluster = "twoway",
                                                 variable = VAR)
i <- i+1

# convert to dataframe to be able to chart
reg_stats_indvar <- bind_rows(reg_stats_indvar_list)

if(sum(duplicated(reg_stats_indvar))==0 & nrow(reg_stats_indvar)+1 == i){
  # save it 
  if(length(ISL) ==1){
    saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",VAR)))
  }else{
    saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_all_",VAR)))
  }
}else{stop("SOMETHING WENT WRONG")}



# find position of model to highlight in original data frame
a <- reg_stats_indvar
model_idx <- a[a$lucpfip_pixelcount_total & 
                 a$CR_30km & 
                 a$imp1 &
                 a$quasipoisson &
                 a$pya_4 &
                 a$lag_or_not &
                 a$unit_distryear_fe &
                 # a$two_commo & 
                 a$control_own & 
                 a$n_reachable_uml_control &
                 a$prex_cpo_control == FALSE &
                 a$pya_outcome_control == FALSE &
                 a$ngb_ov_lag4 == FALSE &
                 a$weights == FALSE &
                 a$oneway_cluster, ] %>% rownames()
model_idx
rm(a)

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
