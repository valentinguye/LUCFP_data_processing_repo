### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#   Some mills have a missing island name, although they are geo-referenced. 
#   This script simply aims at identifying these mills and the island corresponding to their coordinates
#   in order to make these changes in add_variables.do
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "readstata13", 
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



### NEW FOLDERS USED IN THIS SCRIPT 

### INDONESIAN CRS 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### IBS YEARS 
years <- seq(1998, 2015, 1)

ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# ibs <- ibs[!is.na(ibs$lat),]

ibs <- st_as_sf(ibs, coords = c("lon", "lat"), crs = 4326)
ibs <- st_transform(ibs, crs = indonesian_crs)

ibs_mi <- ibs[ibs$island_name=="",]

island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

# make the operation faster by using island bbox (other wise the island polygons make 
# the computation very long)
# (and this also includes parcel centroids in the sea)
island_sf_prj_bbox <- sapply(island_sf_prj$geometry, function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)

sgbp <- st_within(ibs_mi$geometry, island_sf_prj_bbox)

sgbp[lengths(sgbp)>1] 

# island_sf_prj features are in this order : 1 Sumatra; 2 Papua; 3 Kalimantan
ibs_mi$island <- unlist(sgbp)
ibs_mi$island <- replace(ibs_mi$island, ibs_mi$island == 1, "Sumatra")
unique(ibs_mi$island)
ibs_mi$island <- replace(ibs_mi$island, ibs_mi$island == 3, "Kalimantan")

ibs_mi$firm_id[ibs_mi$island == "Sumatra"]
display <- ibs_mi[ibs_mi$island == "Sumatra","firm_id"]
View(display)
