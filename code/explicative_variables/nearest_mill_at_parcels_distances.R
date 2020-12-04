
### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("data.table", "dplyr", "readstata13", "readxl",
                   "rgdal", "sf")
#install.packages("sf", source = TRUE)
# library(sf)
# 
# neededPackages = c("tidyverse","data.table", "readxl","foreign", "data.table", "readstata13", "here",
#                    "rgdal", "raster", "velox","sp", "lwgeom", "rnaturalearth", 
#                    "rlist", "parallel", "foreach", "iterators", "doParallel" )
# 

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
# troublePackages <- c("")
# # Attempt to load packages from user's default libraries. 
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE 
parcel_size <- 3000

# read the sample panel of IBS geolocalized mills
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# make it a cross section
ibs_cs <- ibs[!duplicated(ibs$firm_id),]

# IBS for analysis
ibsana <- ibs_cs[ibs_cs$analysis_sample==1,]


#IBS-UML
# keep only those that are matched with UML. 
ibsuml <- ibs_cs[ibs_cs$uml_matched_sample==1 & ibs_cs$analysis_sample==1,]

ibs_desacentro <- ibs_cs[ibs_cs$uml_matched_sample==0 & ibs_cs$analysis_sample==1,]


# UML
# read the most complete version of UML we have. 
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
# this is already a cross section


# join a firm_id to UML
uml <- left_join(uml, ibsuml[,c("trase_code", "firm_id")], by = "trase_code")

# add ibs not matched with uml (coordinates are desa centroids) 
uml <- dplyr::select(uml, firm_id, lon, lat)
all <- rbind(uml, ibs_desacentro[,c("firm_id","lon","lat")])

# make spatial
all <- st_as_sf(all, coords = c("lon", "lat"), remove = FALSE, crs = 4326)
all <- st_transform(all, crs = indonesian_crs)

### MATCH EACH PARCEL WITH ITS NEAREST MILL
# parcel_size <- 3000
# catchment_radius <- 3e4

for(catchment_radius in c(1e4, 3e4, 5e4)){
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  # this is a data frame that has vocation to get merged with a main data set, so remove extra variables
  parcels <- dplyr::select(parcels, lonlat, year, idncrs_lon, idncrs_lat)
  
  # add three cross sections to the panel, for years 1998-2000. 
  for(y in c(2000, 1999, 1998)){
    cs <- parcels[!duplicated(parcels$lonlat),]
    cs[,"year"] <- y
    parcels <- rbind(cs,parcels)
  }
  rm(cs)
  
  # keep only one cross-section, no matter which. 
  parcels_centro <- parcels[!duplicated(parcels$lonlat),]
  
  # turn it into a sf object (idncrs_lon idncrs_lat are already expressed in indonesian crs)
  parcels_centro <- st_as_sf(parcels_centro, coords = c("idncrs_lon", "idncrs_lat"), remove = TRUE, crs = indonesian_crs)
  
  nearest_mill_idx <- st_nearest_feature(parcels_centro, all)
  
  parcels_centro$nearest_firm_id <- all$firm_id[nearest_mill_idx]
  
  parcels_centro <- st_drop_geometry(parcels_centro)
  
  parcels <- left_join(parcels[,c("lonlat", "year")], parcels_centro[,c("lonlat", "nearest_firm_id")], by = "lonlat")
  
  rm(parcels_centro)
  # select variables we want to "transfer" from mills to their nearest parcels
  # These variables are transfered only if the parcel is nearest to an ibs mill. Not if nearer is a UML mill not matched with ibs. 
  ibsvars <- dplyr::select(ibs,
                           firm_id, year,
                           trase_code, uml_id, mill_name, parent_co,
                           kec_name, village_name, desa_id,desa_code_2000, est_year_imp, 
                           ffb_price_imp1, ffb_price_imp2,
                           cpo_price_imp1, cpo_price_imp2, prex_cpo_imp1, prex_cpo_imp2,
                           pct_own_cent_gov_imp, pct_own_loc_gov_imp, pct_own_nat_priv_imp, pct_own_for_imp,
                           concentration_30, concentration_50
                           )
  
  # let the names like in IBS, they all got the wa_ prefix in wa_at_parcels_distances so there wont be duplicated names. 
  # just rename firm_id
  names(ibsvars)[names(ibsvars)=="firm_id"] <- "nearest_firm_id"
  
  parcels <- left_join(parcels, ibsvars, by = c("nearest_firm_id", "year"))
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/nm_panel_parcels_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
 
  rm(parcels, nearest_mill_idx) 
}

rm(ibs, ibsuml, ibsvars, ibsana, uml, all, ibs_cs, ibs_desacentro)


