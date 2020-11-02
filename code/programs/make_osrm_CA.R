
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

# preliminary: read the island shapefile
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)


to_panel_within_IBS_CA <- function(island, parcel_size){
  ### Function description

  # to_panel_within_IBS_CA converts the raster bricks of annual layers of parcels, for different outcomes in cell values, to a panel data frame.
  # This is executed for each forest, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 0. Preparing the mill data (here, geolocalized IBS mills)
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. Selecting parcels within 80km to a mill at least one year
  # 3. Requesting the local OSRM server for the driving durations between all mill-parcel pairs 
  # 4. Switching from raster cells to data frame cells, selecting parcels within travel times and reshaping to long format. 
  
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
  
  
  ## 1. Masking.
  # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
  parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_total")
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))

  mask(x = parcels_brick, mask = total_ibs_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(parcels_brick, total_ibs_ca_sp)
  
  
  ## 2. Selecting parcels within 80km to a mill at least one year
  # (i.e. the parcel is present in the data frame in all years even if it is within 80km of a mill only since 2014)
  
  # Turn the masked raster to a sf data frame
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))
  
  # note the na.rm = TRUE 
  ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  ibs_msk_df <- ibs_msk_df %>% dplyr::rename(lon = x, lat = y)
  ibs_msk_df <- st_as_sf(ibs_msk_df, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
  
  island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
  ibs_msk_df$parcel_id <- paste0(island_id, c(1:nrow(ibs_msk_df))) %>% as.numeric()
  
  
  ## 3. REQUESTING THE OSRM DURATIONS 
  # transform to lon-lat and to Spatial is necessary for osrm v.3.0.0
  ibs_msk_df_lonlat <- st_transform(ibs_msk_df, crs = 4326) 
  ibs_msk_df_lonlat_sp <- as(ibs_msk_df_lonlat, "Spatial")
  ibs_sp <- as(ibs, "Spatial")
  
  # See the notes below and in Evernote on OSRM
  osrmr::run_server(osrm_path = osrm_path, map_name = map_name)
  
  dur_list <-  osrmTable(src = ibs_msk_df_lonlat_sp, dst = ibs_sp)
  
  osrmr::quit_server()
  
  saveRDS(dur_list, "input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_IBS")

  rm(dur_list, ibs_msk_df, ibs_msk_df_lonlat, ibs_msk_df_lonlat_sp)
  
 
  
  ## 4. Selecting parcels within travel times, and reshaping, for all different outcomes (lucpfip, lucfsmp, lucfip, lucfsmp)
  
  select_within_TT <- function(forest, size){
    
    if(forest %in% c("intact", "degraded", "total")){
      lucf <- "lucpf"
    }
    if(forest %in% c("30th", "60th", "90th")){
      lucf <- "lucf"
    }
    
    # Turn the masked raster to a sf data frame
    parcels_brick_name <- paste0("parcel_",lucf,size,"p_",island,"_",parcel_size/1000,"km_",forest)
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))

    # note the na.rm = TRUE 
    ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
    
    ibs_msk_df <- ibs_msk_df %>% dplyr::rename(lon = x, lat = y)
    
    # read in the duration matrix, common to all parcel layers (outcomes) 
    dur_mat <- readRDS("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_IBS")
    dur_mat <- dur_mat$durations
    
    # Durations are in minutes, and we impose a restriction of no more than 2, 4 or 6 hours of travel between plantation and mill.
    for(travel_time in c(2,4,6)){

      # keep only the parcels that satisfy the travel time condition for at least one mill in the whole period.
      dur_mat_log <- dur_mat/(60) < travel_time
      
      ibs_msk_TT_df <- ibs_msk_df[base::rowSums(dur_mat_log, na.rm = TRUE)>0,]
      # na.rm = TRUE is in case of mills for which no duration could be computed. 
      # which is the case for 4 mills in Sumatra, 3 of which are not UML matched but desa centroid located.
      
      # dur_mat[is.na(dur_mat[1:3,])]
      # ibs[ibs$firm_id %in% c(4096, 4097,51432,54268),"uml_matched_sample"]
      # ibs[ibs$firm_id==51432,]
      # ibs[ibs$uml_matched_sample ==0,] %>% nrow()
      # is.na(dur_mat[1:3,]) %>% sum()
      
      # vector of the names in the wide format of our time varying variables
      # the column names are the layer names in the parcels_brick + the layer index, separated by "." see examples in raster::as.data.frame
      # the layer index (1-18) indeed corresponds to the year (2001-2018) because, in aggregate_lucpfip, at the end, 
      # rasterlist is ordered along 2001-2018 and this order is preserved when the layers are bricked. 
      varying_vars <- paste0(parcels_brick_name, "_IBS_masked.", seq(from = 1, to = 18))
      
      # reshape to long
      ibs_msk_TT_long_df <- stats::reshape(ibs_msk_TT_df,
                                           varying = varying_vars,
                                           v.names = paste0(lucf,size,"p_pixelcount_",forest),
                                           sep = ".",
                                           timevar = "year",
                                           idvar = c("parcel_id"), # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                                           ids = "parcel_id",
                                           direction = "long")
                    
      rm(varying_vars, ibs_msk_TT_df)
      # replace the indices from the raster::as.data.frame with actual years.
      years <- seq(from = 2001, to = 2018, by = 1)
      ibs_msk_TT_long_df <- mutate(ibs_msk_TT_long_df, year = years[year])
      
      ibs_msk_TT_long_df <- dplyr::arrange(ibs_msk_TT_long_df, parcel_id, year)
      
      # Add columns of converted pixel counts to hectares.
      ibs_msk_TT_long_df$inha <- ibs_msk_TT_long_df[,grepl("pixelcount",colnames(ibs_msk_TT_long_df))]*(27.8*27.6)/(1e4)
      
      names(ibs_msk_TT_long_df[names(ibs_msk_TT_long_df)=="inha"]) <- paste0(lucf,size,"p_ha_",forest)
      
      
      saveRDS(ibs_msk_TT_long_df,
              file = file.path(paste0("temp_data/processed_parcels/",lucf,size,"p_panel_",
                                      island,"_",
                                      parcel_size/1000,"km_",
                                      travel_time,"h_IBS_CA_", 
                                      forest,".rds")))
      rm(dur_mat_log)
    }

    rm(dur_mat, ibs_msk_df_lonlat)
  }
  
  # Loop it over outcomes
  forestS <- c("intact", "degraded", "total")# "30th", "60th", "90th",
  for(forest in forestS){
    for(size in c("i","s", "m")){
      select_within_TT(forest = forest,
                   size = size)
      
      removeTmpFiles(h=0)
      
    }
  }  
}




to_panel_within_UML_CA <- function(island, parcel_size){
  ### Function description
  
  # to_panel_within_UML_CA converts the raster bricks of annual layers of parcels, for different outcomes in cell values, to a panel data frame.
  # This is executed for each forest, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 0. Preparing the mill data (here, UML mills)
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. Selecting parcels within 80km to a mill at least one year
  # 3. Requesting the local OSRM server for the driving durations between all mill-parcel pairs 
  # 4. Switching from raster cells to data frame cells, selecting parcels within travel times and reshaping to long format. 
  
  ## 0. Preparing
  
  uml <- read_xlsx(file.path("input_data/uml/mills_20200129.xlsx"))
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
  
  
  ## 1. Masking.
  # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
  parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_total")
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
  
  mask(x = parcels_brick, mask = total_uml_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(parcels_brick, total_uml_ca_sp)
  
  
  ## 2. Selecting parcels within 80km to a mill at least one year
  # (i.e. the parcel is present in the data frame in all years even if it is within 80km of a mill only since 2014)
  
  # Turn the masked raster to a sf data frame
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")))
  
  # note the na.rm = TRUE 
  uml_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  uml_msk_df <- uml_msk_df %>% dplyr::rename(lon = x, lat = y)
  uml_msk_df <- st_as_sf(uml_msk_df, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
  
  island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
  uml_msk_df$parcel_id <- paste0(island_id, c(1:nrow(uml_msk_df))) %>% as.numeric()
  
  
  ## 3. REQUESTING THE OSRM DURATIONS 
  # transform to lon-lat and to Spatial is necessary for osrm v.3.0.0
  uml_msk_df_lonlat <- st_transform(uml_msk_df, crs = 4326) 
  uml_msk_df_lonlat_sp <- as(uml_msk_df_lonlat, "Spatial")
  uml_sp <- as(uml, "Spatial")
  
  # See the notes below and in Evernote on OSRM
  osrmr::run_server(osrm_path = osrm_path, map_name = map_name)
  
  dur_list <-  osrmTable(src = uml_msk_df_lonlat_sp, dst = uml_sp)
  
  osrmr::quit_server()
  
  saveRDS(dur_list, "input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_UML")
  
  rm(dur_list, uml_msk_df, uml_msk_df_lonlat, uml_msk_df_lonlat_sp)



  ## 4. Selecting parcels within travel times, and reshaping, for all different outcomes (lucpfip, lucfsmp, lucfip, lucfsmp)
  
  select_within_TT <- function(forest, size){
    
    if(forest %in% c("intact", "degraded", "total")){
      lucf <- "lucpf"
    }
    if(forest %in% c("30th", "60th", "90th")){
      lucf <- "lucf"
    }
    
    # Turn the masked raster to a sf data frame
    parcels_brick_name <- paste0("parcel_",lucf,size,"p_",island,"_",parcel_size/1000,"km_",forest)
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")))
    
    # note the na.rm = TRUE 
    uml_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
    
    uml_msk_df <- uml_msk_df %>% dplyr::rename(lon = x, lat = y)
    
    # read in the duration matrix, common to all parcel layers (outcomes) 
    dur_mat <- readRDS("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_UML")
    dur_mat <- dur_mat$durations
    
    # Durations are in minutes, and we impose a restriction of no more than 2, 4 or 6 hours of travel between plantation and mill.
    for(travel_time in c(2,4,6)){
      
      # keep only the parcels that satisfy the travel time condition for at least one mill in the whole period.
      dur_mat_log <- dur_mat/(60) < travel_time
      
      uml_msk_TT_df <- uml_msk_df[base::rowSums(dur_mat_log, na.rm = TRUE)>0,]
      # na.rm = TRUE is in case of mills for which no duration could be computed. 
      # which is the case for 4 mills in Sumatra, 3 of which are not UML matched but desa centroid located. 
      
      # vector of the names in the wide format of our time varying variables
      # the column names are the layer names in the parcels_brick + the layer index, separated by "." see examples in raster::as.data.frame
      # the layer index (1-18) indeed corresponds to the year (2001-2018) because, in aggregate_lucpfip, at the end, 
      # rasterlist is ordered along 2001-2018 and this order is preserved when the layers are bricked. 
      varying_vars <- paste0(parcels_brick_name, "_UML_masked.", seq(from = 1, to = 18))
      
      # reshape to long
      uml_msk_TT_long_df <- stats::reshape(uml_msk_TT_df,
                                           varying = varying_vars,
                                           v.names = paste0(lucf,size,"p_pixelcount_",forest),
                                           sep = ".",
                                           timevar = "year",
                                           idvar = c("parcel_id"), # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                                           ids = "parcel_id",
                                           direction = "long")
      
      rm(varying_vars, uml_msk_TT_df)
      # replace the indices from the raster::as.data.frame with actual years.
      years <- seq(from = 2001, to = 2018, by = 1)
      uml_msk_TT_long_df <- mutate(uml_msk_TT_long_df, year = years[year])
      
      uml_msk_TT_long_df <- dplyr::arrange(uml_msk_TT_long_df, parcel_id, year)
      
      # Add columns of converted pixel counts to hectares.
      uml_msk_TT_long_df$inha <- uml_msk_TT_long_df[,grepl("pixelcount",colnames(uml_msk_TT_long_df))]*(27.8*27.6)/(1e4)
      
      names(uml_msk_TT_long_df[names(uml_msk_TT_long_df)=="inha"]) <- paste0(lucf,size,"p_ha_",forest)
      
      
      saveRDS(uml_msk_TT_long_df,
              file = file.path(paste0("temp_data/processed_parcels/",lucf,size,"p_panel_",
                                      island,"_",
                                      parcel_size/1000,"km_",
                                      travel_time,"h_UML_CA_", 
                                      forest,".rds")))
      rm(dur_mat_log)
    }
    
    rm(dur_mat, uml_msk_df_lonlat)
  }

  # Loop it over outcomes
  forestS <- c("intact", "degraded", "total")# "30th", "60th", "90th",
  for(forest in forestS){
    for(size in c("i","s", "m")){
      select_within_TT(forest = forest,
                       size = size)
      
      removeTmpFiles(h=0)
      
    }
  }  
}





forest <- "total"
size <- "s"
island <- "Kalimantan"
parcel_size <- 3000
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### For each island and for each aggregation factor, extract panels of parcels within catchment areas of IBS or UML mills.
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  to_panel_within_IBS_CA(island = Island,
                     parcel_size = PS)
    
}

PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  to_panel_within_UML_CA(island = Island,
                     parcel_size = PS)
  
}

#### Gather the lucfip variables for each parcel_size and catchment_radius combinations. ####
PS <- 3000  
sampleS <- c("IBS", "UML")
travel_timeS <- c(2,4,6)
for(size in c("i","s", "m")){
  for(sample in sampleS){
    for(TT in travel_timeS){
      # Repeat for primary and normal forest. 
      ## PRIMARY
      # For each Island, join columns of lucfip variable for different forest definitions 
      pf_df_list <- list()
      IslandS <- c("Sumatra", "Kalimantan", "Papua")
      for(Island in IslandS){
        
        df_intact   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpf",size,"p_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_intact.rds")))
        df_degraded <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpf",size,"p_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_degraded.rds")))
        df_total    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpf",size,"p_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_total.rds")))
        
        df_degraded <- dplyr::select(df_degraded, -lon, -lat)
        df <- inner_join(df_intact, df_degraded, by = c("parcel_id", "year"))
        
        df_total <- dplyr::select(df_total, -lon, -lat)
        pf_df_list[[match(Island, IslandS)]] <- inner_join(df, df_total, by = c("parcel_id", "year"))
      }
      
      # stack the three Islands together
      indo_df <- bind_rows(pf_df_list)
      
      indo_df <- dplyr::select(indo_df, parcel_id, year,
                               everything())
      
      saveRDS(indo_df, file.path(paste0("input_data/processed_parcels/lucpf",size,"p_panel_",PS/1000,"km_",TT,"h_",sample,"_CA.rds")))
      
      rm(indo_df, pf_df_list)
      
      
      
      ## "NORMAL" forest
      df_list <- list()
      IslandS <- c("Sumatra", "Kalimantan", "Papua")
      for(Island in IslandS){
        
        df_90th <- readRDS(file.path(paste0("temp_data/processed_parcels/lucf",size,"p_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_90th.rds")))
        df_60th <- readRDS(file.path(paste0("temp_data/processed_parcels/lucf",size,"p_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_60th.rds")))
        df_30th <- readRDS(file.path(paste0("temp_data/processed_parcels/lucf",size,"p_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_30th.rds")))
        
        df_60th <- dplyr::select(df_60th, -lon, -lat)
        df <- inner_join(df_90th, df_60th, by = c("parcel_id", "year"))
        
        df_30th <- dplyr::select(df_30th, -lon, -lat)
        df_list[[match(Island, IslandS)]] <- inner_join(df, df_30th, by = c("parcel_id", "year"))
      }
      
      # stack the three Islands together
      indo_df <- bind_rows(df_list)
      
      indo_df <- dplyr::select(indo_df, parcel_id, year,
                               everything())
      
      saveRDS(indo_df, file.path(paste0("input_data/processed_parcels/lucf",size,"p_panel_",PS/1000,"km_",TT,"h_",sample,"_CA.rds")))
      
      rm(indo_df, df_list)
    }
  }
}