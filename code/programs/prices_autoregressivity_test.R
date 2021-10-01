

### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table",
                   "foreign", "readstata13", "readxl",
                   "raster", "rgdal",  "sp", "spdep", "sf",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "msm", "car", "fixest", "sandwich", "lmtest", "boot", "multcomp",
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


d_50 <- readRDS(file.path(paste0("temp_data/store_final_data_server/panel_parcels_ip_final_",
                                 parcel_size/1000,"km_",
                                 "50CR.rds")))
# Split them into islands of interest
d_50_suma <- d_50[d_50$island == "Sumatra",]
d_50_kali <- d_50[d_50$island == "Kalimantan",]
#d_50_papu <- d_50[d_50$island == "Papua",]


d_30 <- readRDS(file.path(paste0("temp_data/store_final_data_server/panel_parcels_ip_final_",
                                 parcel_size/1000,"km_",
                                 "30CR.rds")))
d_30 <- d_30[,names(d_50)]
# Split them into islands of interest
d_30_suma <- d_30[d_30$island == "Sumatra",]
d_30_kali <- d_30[d_30$island == "Kalimantan",]
#d_30_papu <- d_30[d_30$island == "Papua",]

rm(d_30, d_50)

parcels <- rbind(d_30_suma, d_50_kali)
rm(d_30_suma, d_50_suma, d_30_kali, d_50_kali)

## short to long lags
for(lag in c(1:10)){
  parcels <- dplyr::arrange(parcels, lonlat, year)
  parcels <- DataCombine::slide(parcels,
                                Var = "wa_cpo_price_imp1", 
                                TimeVar = "year",
                                GroupVar = "lonlat",
                                NewVar = paste0("wa_cpo_price_imp1_lag",lag),
                                slideBy = -lag, 
                                keepInvalid = TRUE)
  # parcels <- dplyr::arrange(parcels, lonlat, year)
  
}

model <- as.formula(paste0("wa_cpo_price_imp1 ~ wa_cpo_price_imp1_lag1 + wa_cpo_price_imp1_lag2 + wa_cpo_price_imp1_lag3 + wa_cpo_price_imp1_lag4",
                    " + wa_cpo_price_imp1_lag5  + wa_cpo_price_imp1_lag6 + wa_cpo_price_imp1_lag7 + ",
                    "wa_cpo_price_imp1_lag8 + wa_cpo_price_imp1_lag9 + wa_cpo_price_imp1_lag10 | lonlat + year"))

res4 <- feols(model, data = parcels, 
             notes = TRUE)

model_4pya <- as.formula("wa_cpo_price_imp1 ~ cpo_price_imp1_yoyg_4ya_lag1 | lonlat")

res <- feols(model_4pya, data = parcels, 
                notes = TRUE)