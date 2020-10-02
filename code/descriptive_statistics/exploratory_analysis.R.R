
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
# PB: on n'arrive pas à installed leaflet, kableExtra (qui n'est peut être pas nécessaire) et velox (qui n'est pas nécessaire)
neededPackages = c("dplyr", "data.table",  "stringr", "sjmisc", 
                   "foreign", "readxl", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "parallel", "foreach","doParallel")
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
troublePackages <- c("plyr",  
                     "tidyr", "dplyr", "ggplot2", "corrplot", "lfe", "multiwayvcov", "lmtest", "stargazer", "scales", "shiny", "DT", "openssl", "tictoc", "shinycssloaders", "kableExtra", "rio", "zip", "rlang", 
                     "ExPanDaR", "packrat", "rsconnect")
# Attempt to load packages from user's default libraries.
lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 

### NEW FOLDERS USED IN THIS SCRIPT 

### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### IBS ####  
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
any(duplicated(ibs[,c("firm_id", "year")]))


# LOG PRICES
ibs <- dplyr::mutate(ibs, 
                     ln_ffb_price_imp1 = log(ffb_price_imp1),
                     ln_ffb_price_imp2 = log(ffb_price_imp2),
                     ln_cpo_price_imp1 = log(cpo_price_imp1),
                     ln_cpo_price_imp2 = log(cpo_price_imp2),
                     ln_pko_price_imp1 = log(pko_price_imp1),
                     ln_pko_price_imp2 = log(pko_price_imp2))

# pre-select variables before exploring the dataset with ExPanD 
# (this is the same set as the one selected in wa_at_parcels.R + geo_sample and is_mill + geographic variables)
ibs <- ibs[, c("firm_id", "year", "is_mill", "geo_sample", "analysis_sample",
               "trase_code", "uml_id", "mill_name", "parent_co", "lat", "lon",
               "island_name", "district_name", "kec_name", "village_name", 
               "min_year","est_year", "est_year_imp", "max_year", 
               "ffb_price_imp1", "ffb_price_imp2", "in_ton_ffb_imp1", "in_ton_ffb_imp2", "in_val_ffb_imp1", "in_val_ffb_imp2",
               "cpo_price_imp1","cpo_price_imp2", "out_ton_cpo_imp1", "out_ton_cpo_imp2", "out_val_cpo_imp1", "out_val_cpo_imp2",
               "prex_cpo_imp1", "prex_cpo_imp2",
               "pko_price_imp1","pko_price_imp2", "out_ton_pko_imp1", "out_ton_pko_imp2", "out_val_pko_imp1", "out_val_pko_imp2",
               "prex_pko_imp1", "prex_pko_imp2",
               "export_pct_imp", "revenue_total", "workers_total_imp3",
               "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp", 
               "iv2_imp1", "iv2_imp2", "iv3_imp1", "iv3_imp2", "iv4_imp1", "iv4_imp2", 
               "concentration_10", "concentration_30", "concentration_50",
               "ln_ffb_price_imp1", "ln_ffb_price_imp2","ln_cpo_price_imp1", "ln_cpo_price_imp2","ln_pko_price_imp1", "ln_pko_price_imp2")]
ExPanD(df = ibs, cs_id = "firm_id", ts_id = "year")



##### PARCEL PANEL - LUCPFIP #### 
parcel_size <- 3000
catchment_radius <- 3e4

parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
names(parcels)

# number of distinct parcels in the different samples 
# (values reported here in code for 30km catchment radius)
length(unique(parcels$parcel_id))
length(unique(parcels$parcel_id[parcels$any_fc2000_30th]))
length(unique(parcels$parcel_id[parcels$any_pfc2000_total]))

# We moght want to restrict our sample to those parcels that have at least one year some outcome variable positive, 
# because units that have a null sum of outcomes over time do not contribute to the estimation (Wooldridge 2002)

# (~15 minutes)
any_outcome <- ddply(parcels, "parcel_id", summarise, 
                      any_lucfip_30th = sum(lucfip_ha_30th, na.rm = TRUE)>0, 
                      any_lucpfip_total = sum(lucpfip_ha_total, na.rm = TRUE)>0)

any_outcome$any_lucpfip_total %>% sum()
any_outcome$any_lucfip_30th %>% sum()

parcels <- merge(parcels, any_outcome, by = "parcel_id")


# select here the df you want to explore in ExPanD

parcels1 <- parcels[,c("parcel_id", "year", "island", "district",
                      "lucpfip_ha_intact", "lucpfip_ha_total",
                      "lucfip_ha_90th", "lucfip_ha_60th", "lucfip_ha_30th", 
                      "wa_est_year_imp", 
                      "wa_ffb_price_imp1",
                      "wa_cpo_price_imp1", 
                      "wa_pko_price_imp1",
                      "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                      paste0("wa_concentration_",catchment_radius/1000), 
                      paste0("n_reachable_uml_",catchment_radius/1000,"km"),
                      paste0("sample_coverage_",catchment_radius/1000,"km"))]



parcels1 <- parcels[parcels$any_lucfip_30th,] # this removes also parcels that have NA in 1998-2000 and then only zeros.)
length(unique(parcels1$parcel_id))

ExPanD(df = parcels1, cs_id = "parcel_id", ts_id = "year")

