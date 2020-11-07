
# This script has in two main parts. 
# The first part selects cells within catchment area (CA) for outcome variables, not annually. 
# The second part does the same but annually (the annual variation is whether mills are there or not), for explanatory variables (computed in wa_at_parcels_durations.R)

# The first part is repeated for the IBS and the UML sets of mills.  



# Framed below, are notes on the procedure to set this script up and running. 

### -------------------------------------------------------------------------------------------------------------------------------------------------###
# We want a local instance of OSRM because it enables to compute large duration matrices efficiently. 
# Slicing the inputs and the task and using the public server is easily coded, but very inefficient as many requests need to be sent (it's in terms of days of running times in our case). 
# 
# Download OSRM Release build (February 01 2018 13:27:05 CEST) from here http://build.project-osrm.org/
# Put the folder content in a folder at "C:/Users/GUYE/osrm" (or where ever)
# download indonesia-latest.osrm.pbf at https://download.geofabrik.de/asia/indonesia.html on 29/10/2020
# put the downloaded data in the osrm folder. 
# open command prompt and run: 
#     cd C:\Users\GUYE\osrm
#     osrm-extract.exe indonesia-latest.osm.pbf -p car.lua
#     osrm-contract.exe indonesia-latest.osrm
#     osrm-routed.exe indonesia-latest.osrm --port 5000 --max-table-size=1000000000 
# 
# The first line (after cd...) takes a bit of time bc it compiles the data for indonesia, the routes etc. 
# The last one is useless but if it returns  "running and waiting for requests" this means the installation of the local OSRM instance is achieved. 
# 
# Now, in R, first thing is to install {osrm} and {osrmr}. 
# The former executes osrmTable, and for it to work, it needs to be in an old version, like 3.0.0, hence 
# remotes::install_version("osrm", version= "3.0.0"). 
# otherwise there is an error from no compatibility with the local OSRM version downloaded. See there https://github.com/rCarto/osrm/issues/52
# 
# {osrmr} enables to open the shell and get the local osrm server "running and waiting for requests". Typically, we pass arguments to osrmr::run_server() so that it prompts
# cd C:\Users\GUYE\osrm
# osrm-routed.exe indonesia-latest.osrm --port 5000 --max-table-size=1000000000 
# from R. 
# 
# (The  --max-table-size=1000000000 is allowing for a large duration table to be returned from the local server.)
# 
# Hence R code 
# options(osrm.server = paste0(osrmr:::server_address(TRUE), "/"), osrm.profile = "driving")
#     map_name = "indonesia-latest.osrm --port 5000 --max-table-size=1000000000"
#     osrm_path = "C:/Users/GUYE/osrm"
#     Sys.setenv("OSRM_PATH"=osrm_path) 
# 
#     osrmr::run_server(map_name = map_name)
#     
#     dur_list <-  osrmTable(src = m.df_wide_lonlat_sp, dst = mills_sp)
# 
#     osrmr::quit_server()
# 
# Where the inputs of osrmTable need to be in Spatial class in the 3.0.0 version of osrm needed here. 
# 
# 
# # this procedure is inspired from 
# https://phabi.ch/2020/05/06/run-osrm-in-docker-on-windows/
# http://build.project-osrm.org/
# https://reckoningrisk.com/OSRM-server/

### -------------------------------------------------------------------------------------------------------------------------------------------------###




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
dir.create("input_data/local_osrm_outputs")


### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### Prepare polygons of three Indonesian islands of interest ###
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)


parcel_size <- 3000




#### OSRM GENERAL DURATIONS #### 

for(island in c("Sumatra", "Kalimantan")){ #"Papua" does not work with IBS mills (probably not enough routes)
  ### IBS
  ## 0. Preparing

  ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
  # keep only a cross section of those that are geolocalized mills, on island of interest
  ibs <- ibs[ibs$analysis_sample==1,]
  ibs <- ibs[!duplicated(ibs$firm_id),]
  ibs <- ibs[ibs$island_name == island,]
  #turn into an sf object.
  ibs <- st_as_sf(ibs,	coords	=	c("lon",	"lat"), remove = FALSE, crs=4326)
  # keep only the geometry, we do not need mills attributes here.
  ibs_geom <- st_geometry(ibs)
  # set CRS and project
  ibs_prj <- st_transform(ibs_geom, crs = indonesian_crs)
  #define big catchment areas to have a large AOI.
  ibs_ca <- st_buffer(ibs_prj, dist = 80000) # 80km. The point is not to be too restrictive here.
  # work with squares rather than with circles
  for(i in 1:length(ibs_ca)){
    ibs_ca[i] <- st_as_sfc(st_bbox(ibs_ca[i]))
  }
  total_ibs_ca <- st_union(st_geometry(ibs_ca))
  # coerce to a SpatialPolygon
  total_ibs_ca_sp <- as(total_ibs_ca, "Spatial")
  # keep ibs_prj we need it below
  rm(total_ibs_ca, ibs_ca, ibs_geom)


  ## 1. Masking one typical raster
  # (they are all the same dimensions across the measured outcomes - lucpfip, lucfip, lucfsp, etc.)
  # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
  parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_total")
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))

  mask(x = parcels_brick, mask = total_ibs_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)

  rm(parcels_brick)

  # Turn the masked raster to a sf data frame
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))

  # note the na.rm = TRUE. Due to the mask, this keeps in the df object only the parcel within 80km of a mill at least one year.
  ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)

  ibs_msk_df <- ibs_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
  ibs_msk_df <- st_as_sf(ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  ibs_msk_df <- st_transform(ibs_msk_df, crs = 4326)
  ibs_msk_df$lon <- st_coordinates(ibs_msk_df)[,"X"]
  ibs_msk_df$lat <- st_coordinates(ibs_msk_df)[,"Y"]

  ## 2. REQUESTING THE OSRM DURATIONS
  # give row names to source and destinations data to retrieve them more surely
  row.names(ibs_msk_df) <- mutate(ibs_msk_df, lonlat = paste0(lon, lat))$lonlat
  row.names(ibs) <- ibs$firm_id

  # transform to lon-lat and to Spatial is necessary for osrm v.3.0.0
  #ibs_msk_df <- st_transform(ibs_msk_df, crs = 4326)
  ibs_msk_df_sp <- as(ibs_msk_df, "Spatial")
  ibs_sp <- as(ibs, "Spatial")

  # See the notes below and in Evernote on OSRM
  osrmr::run_server(osrm_path = osrm_path, map_name = map_name)

  dur_list <-  osrmTable(src = ibs_msk_df_sp, dst = ibs_sp)

  osrmr::quit_server()

  saveRDS(dur_list, file = file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_IBS")))

  rm(dur_list, ibs_msk_df, ibs_msk_df_sp, ibs_sp)



  ### UML
  
  ## 0. Preparing
  
  uml <- read_xlsx(file.path("input_data/uml/mills_20200129.xlsx"))
  uml <- uml %>% as.data.frame()
  uml$latitude <- as.numeric(uml$latitude)
  uml$longitude <- as.numeric(uml$longitude)
  uml$lat <- uml$latitude
  uml$lon <- uml$longitude
  uml <- st_as_sf(uml,	coords	=	c("longitude",	"latitude"), crs = 4326)
  uml_geom <- st_geometry(uml)
  uml_prj <- st_transform(uml_geom, crs = indonesian_crs)
  
  # there is no island column in this dataset, hence we select mills on the specific island geographically
  select_within <- st_within(x = uml_prj, y = island_sf_prj[island_sf_prj$shape_des == island,])
  uml_prj <- uml_prj[lengths(select_within)>0,]
  uml <- uml[lengths(select_within)>0,]
  
  #define big catchment areas to have a large AOI.
  uml_ca <- st_buffer(uml_prj, dist = 80000)
  # work with squares rather than with circles
  for(i in 1:length(uml_ca)){
    uml_ca[i] <- st_as_sfc(st_bbox(uml_ca[i]))
  }
  total_uml_ca <- st_union(st_geometry(uml_ca))
  # coerce to a SpatialPolygon
  total_uml_ca_sp <- as(total_uml_ca, "Spatial")
  # keep uml_prj we need it below
  rm(total_uml_ca, uml_ca, uml_geom)
  
  
  ## 1. Masking one typical raster 
  # (they are all the same dimensions across the measured outcomes - lucpfip, lucfip, lucfsp, etc.)
  # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
  parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_total")
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
  
  mask(x = parcels_brick, mask = total_uml_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(parcels_brick)
  
  # Turn the masked raster to a sf data frame
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")))
  
  # note the na.rm = TRUE. Due to the mask, this keeps in the df object only the parcel within 80km of a mill at least one year. 
  uml_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  uml_msk_df <- uml_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
  uml_msk_df <- st_as_sf(uml_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  uml_msk_df <- st_transform(uml_msk_df, crs = 4326)
  uml_msk_df$lon <- st_coordinates(uml_msk_df)[,"X"]
  uml_msk_df$lat <- st_coordinates(uml_msk_df)[,"Y"]
  
  ## 2. REQUESTING THE OSRM DURATIONS 
  # give row names to source and destinations data to retrieve them more surely
  row.names(uml_msk_df) <- mutate(uml_msk_df, lonlat = paste0(lon, lat))$lonlat
  row.names(uml) <- uml$trase_code
  
  # transform to lon-lat and to Spatial is necessary for osrm v.3.0.0
  #uml_msk_df <- st_transform(uml_msk_df, crs = 4326) 
  uml_msk_df_sp <- as(uml_msk_df, "Spatial")
  uml_sp <- as(uml, "Spatial")

  # See the notes below and in Evernote on OSRM
  osrmr::run_server(osrm_path = osrm_path, map_name = map_name)
  
  dur_list <-  osrmTable(src = uml_msk_df_sp, dst = uml_sp)
  
  osrmr::quit_server()
  
  saveRDS(dur_list, file = file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_UML")))
  
  rm(dur_list, uml_msk_df, uml_msk_df_sp, uml_sp)
  
}

island <- "Sumatra"
travel_time <- 6
t <- 1998


#### OSRM ANNUAL DURATIONS ####  
# this part cannot be executed right after the part above, because make_osrm_CA.R needs to be executed first, and it 
# requires the outputs from the part above. 
### IBS YEARS
years <- seq(from = 1998, to = 2015, by = 1)

for(island in c("Sumatra", "Kalimantan")){ # , "Papua"
  
    ### PREPARE IBS DATA 
    ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))  
    # keep only geolocalized mills and useful variables
    ibs <- ibs[ibs$analysis_sample == 1,c("firm_id", "year", "min_year", "max_year", "island_name", "lon", "lat")]
    # keep all unique records of mills
    #ibs <- ibs[!duplicated(ibs$firm_id),] # this line is for the computation of n_reachable_ibs_imp, not for the main workflow
    ibs <- ibs[ibs$island_name == island,]
    ibs <- st_as_sf(ibs, coords =  c("lon", "lat"), remove = TRUE, crs = 4326)
    
    for(travel_time in c(2,4,6)){

    ### PREPARE PARCELS
    # Import a parcel panel outputed from make_osrm_CA.R for a given island and a given maximum travel time. 
    parcels_centro <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_CA_total.rds")))
    # keep only one cross-section, no matter which. 
    parcels_centro <- parcels_centro[!duplicated(parcels_centro$parcel_id),]
    # turn it into a sf object
    parcels_centro <- st_as_sf(parcels_centro, coords = c("lon", "lat"), remove = FALSE, crs = 4326)

    # give row names to source and destinations data to retrieve them more surely
    row.names(parcels_centro) <- mutate(parcels_centro, lonlat = paste0(lon, lat))$lonlat

    # transform to Spatial is necessary for osrm v.3.0.0
    parcels_centro_sp <- as(parcels_centro, "Spatial")

    osrmr::run_server(osrm_path = osrm_path, map_name = map_name)
    
    for(t in years){#c(2013,2014,2015)
      # this line makes sure that mills that just miss in IBS for a period of time but then appear again, are counted as reachable.
      # ibs_cs <- ibs[ibs$min_year <= t & ibs$max_year >=t,]  # this line is for the computation of n_reachable_ibs_imp, not for the main workflow
      ibs_cs <- ibs[ibs$year == t,]
      
      row.names(ibs_cs) <- ibs_cs$firm_id
      ibs_cs_sp <- as(ibs_cs, "Spatial")
      
      dur_list <-  osrmTable(src = parcels_centro_sp, dst = ibs_cs_sp) 
  
      saveRDS(dur_list, file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_",t)))
    
      rm(dur_list, ibs_cs, ibs_cs_sp)
    }  
    
    osrmr::quit_server()
  
    rm(parcels_centro, parcels_centro_sp)
  }
  rm(ibs)  
}
rm(years)







### ALTERNATIVE OSRM WITH PUBLIC SERVER (BUT TAKES A WHILE !)
# dur_mat <- matrix(nrow = nrow(m.df_wide_lonlat), ncol = nrow(mills))
# 
# if(nrow(mills) > 200){
#   # chope the calculation for it to work with the free OSRM SERVER 
#   n <- nrow(mills)
#   third <- trunc(n/3)
#   first_3rd <- 1:third
#   snd_3rd <- (third+1):(third*2)
#   last_3rd <- (third*2+1):n
#   
#   for(srci in 1:3){#nrow(m.df_wide_lonlat)
#     #for(desti in 1:nrow(mills)){
#     dur_mat[srci, first_3rd] <- osrmTable(src = m.df_wide_lonlat[srci,], dst = mills[first_3rd,])$durations# %>% as.data.frame()
#     dur_mat[srci, snd_3rd] <- osrmTable(src = m.df_wide_lonlat[srci,], dst = mills[snd_3rd,])$durations# %>% as.data.frame()
#     dur_mat[srci, last_3rd] <- osrmTable(src = m.df_wide_lonlat[srci,], dst = mills[last_3rd,])$durations# %>% as.data.frame()
#     #names(durations[[i]]) <- colnames(durations[[i]])
#     #}
#   }
# }else{
#   for(srci in 1:nrow(m.df_wide_lonlat)){#
#     dur_mat[srci,] <- osrmTable(src = m.df_wide_lonlat[srci,], dst = mills)$durations
#   }
# }
# osrmTable(src = m.df_wide_lonlat_sp[1,], dst = mills_sp)$durations
# # ou alors : 
# slice_size <- 10000/500
# i <- 1
# while(slice < nrow(m.df_wide_lonlat)){
#   dur_mat[i:(i+slice_size),]  <- osrmTable(src = m.df_wide_lonlat[i:(i+slice_size),], dst = mills)$durations
#   
# }
