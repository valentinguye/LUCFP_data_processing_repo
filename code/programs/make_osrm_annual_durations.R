
# This script repeats roughly twice the same actions, once for parcels within supply sheds of IBS mills only, 
# and once for all UML mills. 

### SELECT RASTER CELLS TO A DATAFRAME
# NOW WE CONSIDER CATCHMENT AREAS BASED ON TRAVEL TIME BETWEEN PLANTATIONS AND MILLS


##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("data.table", "dplyr", "readstata13", "readxl",
                   "raster", "rgdal", "sp", "sf", "osrm", "osrmr",
                   "doParallel", "foreach", "parallel")

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


### SET OSRM SERVER 
# See the notes below and in Evernote on OSRM
options(osrm.server = paste0(osrmr:::server_address(TRUE), "/"), osrm.profile = "driving")
map_name = "indonesia-latest.osrm --port 5000 --max-table-size=1000000000"
osrm_path = "C:/Users/GUYE/osrm"

### NEW FOLDERS USED IN THIS SCRIPT 
# final outputs of this script are stored in input_data because their reproducibility relies on user having a local OSRM server. 
dir.create("input_data/processed_parcels")
dir.create("input_data/local_osrm_outputs")


### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### IBS YEARS
years <- seq(from = 1998, to = 2015, by = 1)


for(island in c("Sumatra", "Kalimantan", "Papua")){
  for(travel_time in c(2,4,6)){
    ### PREPARE IBS DATA 
  
    ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))  
    # keep only geolocalized mills and useful variables
    ibs <- ibs[ibs$analysis_sample == 1,c("firm_id", "year", "island_name", "lon", "lat")]
    ibs <- ibs[ibs$island_name == island,]
    ibs <- st_as_sf(ibs, coords =  c("lon", "lat"), remove = TRUE, crs = 4326)
    ibs <- as(ibs, "Spatial")
  
    ### PREPARE PARCELS
    # Import a parcel panel outputed from make_osrm_CA.R for a given island and a given maximum travel time. 
    parcels_centro <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_CA_total.rds")))
    # keep only one cross-section, no matter which. 
    parcels_centro <- parcels_centro[!duplicated(parcels_centro),]
    # turn it into a sf object (lon lat are already expressed in indonesian crs)
    parcels_centro <- st_as_sf(parcels_centro, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
    parcels_centro <- st_transform(parcels_centro, crs = 4326)
    parcels_centro <- dplyr::select(parcels_centro, parcel_id, geometry)
    parcels_centro <- as(parcels_centro, "Spatial")
    
    osrmr::run_server(osrm_path = osrm_path, map_name = map_name)
    
    for(t in years){
      ibs_cs <- ibs[ibs$year == t,]
      
      dur_list <-  osrmTable(src = parcels_centro, dst = ibs_cs) 
  
      saveRDS(dur_list, "input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_",t)
    }  
    
    osrmr::quit_server()
  
    rm(dur_list,ibs, ibs_cs, parcels_centro, years)
  }
}




