### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
#   PREPARE BASELINE (2000) FOREST EXTENT MAPS AND PARCEL VARIABLES  
# 
#   This is necessary to compute descriptive tables in lucfp_des_stats.R
#   and to compute parcel-level variables of baseline forest extents. 
# 
#   Input: - island polygons from prepare_island_sf
#           ---> temp_data/processed_indonesia_spatial/island_sf
# 
#         RASTER BASELINE FOREST EXTENT PART
#          - GFC data from prepare_gfc.R
#           ---> temp_data/processed_lu/gfc_data_Indonesia_30th.tif and also for 60th and 90th
#
#          - 2000 industrial plantations from prepare_lucfip.R
#           ---> temp_data/processed_lu/austin_ioppm_2000_Sumatra_aligned.tif and also for Kalimantan and Papua
# 
#          - 2000 primary forest from prepare_lucpfip.R
#           ---> temp_data/processed_lu/margono_primary_forest_Sumatra_aligned.tif and also for Kalimantan and Papua
#  
#         DATAFRAME BASELINE FOREST EXTENT PART 
#          - rasters of baseline forest extents prepare in this script, previous part
#           ---> temp_data/processed_lu/gfc_fc2000_Sumatra_30th_prj.tif and also for Kalimantan and Papua
#           ---> margono_total_primary_forest_Sumatra_aligned.tif and also for Kalimantan and Papua
#
#          - IBS and UML dataframes 
#           ---> temp_data/processed_mill_geolocalization/IBS_UML_panel.dta (from merging_geolocalization_works.do)
#           ---> input_data/uml/mills_20200129.xlsx
# 
#   Main Output: 
#          - A cross-section dataframe of parcels, with variables for pixelcounts and areas (ha) of forest extents in 2000 
#             for two forest extents definitions: 
#               1. 30% tree canopy density outside industrial oil palm plantations
#               2. Total primary forest                                  
#       temp_data/processed_parcels/baseline_fc_cs_3km_10km_IBS_CR.rds and also for 30km and 50km, and UML
# 
# 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
# PB: on n'arrive pas à installed leaflet, kableExtra (qui n'est peut être pas nécessaire) et velox (qui n'est pas nécessaire)
neededPackages = c("dplyr", "data.table", 
                   "foreign", "readxl", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
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

# # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("")
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

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

### Prepare polygons of three Indonesian islands of interest ###
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




##### 1. FORET COVER IN 2000 ##### 

### Prepare 30, 60, 90% forest cover outside industrial plantations in 2000
prepare_fc2000 <- function(island){
  
  aoi <- island_sf[island_sf$shape_des == island,]
  aoi <- st_as_sfc(st_bbox(aoi))
  aoi_sp <- as(aoi, "Spatial")
  rm(aoi)
  # it is not projected # 
  
  ### Prepare forest cover at 30, 60 and 90% canopy closure thresholds
  
  thresholdS <- c(30, 60, 90)
  for(th in thresholdS){
    # Import gfc data prepare (Downloaded, extracted, thresholded with gfcanalysis package in prepare_gfc.R)
    thed_gfc_data <- brick(file.path(paste0("temp_data/processed_lu/gfc_data_Indonesia_",th,"th.tif"))) 
    
    # Select the 2000 forest cover layer: this is the first one, as described in the gfcanalysis package documentation. 
    # It is already coded in 1 (forest, i.e. canopy closure > threshold) and 0 (non forest i.e. canopy closure < threshold) 
    
    fc2000 <- thed_gfc_data[[1]]
    
    rm(thed_gfc_data)
    
    ## Crop to extent of island 
    crop(fc2000, aoi_sp, 
         filename = file.path(paste0("temp_data/processed_lu/gfc_fc2000_",island,"_",th,"th.tif")),
         datatype = "INT1U",
         overwrite = TRUE)
    
    
    ## Project fc2000 layer
    
    # This is necessary because the projected data are the one used for analysis
    
    fc2000 <- raster(file.path(paste0("temp_data/processed_lu/gfc_fc2000_",island,"_",th,"th.tif")))
    
    beginCluster() # this uses by default detectCores() - 1 
    
    projectRaster(fc2000,
                  method = "ngb",
                  crs = indonesian_crs, 
                  filename = file.path(paste0("temp_data/processed_lu/gfc_fc2000_",island,"_",th,"th_prj.tif")), 
                  datatype = "INT1U",
                  overwrite = TRUE)
    
    endCluster()
    
    rm(fc2000) 
    removeTmpFiles(h=0)
    
    print(paste0("complete ", "temp_data/processed_lu/gfc_fc2000_",island,"_",th,"th_prj.tif"))
    
  }

  ### Overlay 2000 forest cover and 2000 industrial plantations.
  
  ioppm2000 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif")))
  
  # overlay function
  overlay_maps <- function(rs){rs[[1]]*(1-rs[[2]])}
  # multiplies a cell of 2000 forest cover (rs[[1]]) by 0 (i.e. "removes" it) if it it is a plantation in 2000 (rs[[2]]) or if is not a plantation in 2015 (rs[[3]])
  
  ## For each threshold, overlay forest fc2000 map with 2000 industrial plantation maps in a clusterR setting 
  th <- 30
  while(th < 100){
    # call the fc2000 layer for threshold th
    fc2000 <- raster(file.path(paste0("temp_data/processed_lu/gfc_fc2000_",island,"_",th,"th_prj.tif")))
    
    # stack fc2000 with plantation maps (necessary for clusterR)
    rs <- stack(fc2000, ioppm2000)
    
    # run the computation in parallel with clusterR, as cells are processed one by one independently.
    beginCluster() # uses by default detectedCores() - 1
    clusterR(rs,
             fun = calc, # note we use calc but this is equivalent to using overlay 
             # (but more appropriate to the input being a stack)
             args = list(overlay_maps),
             filename = file.path(paste0("temp_data/processed_lu/gfc_fc2000_outside_ip_",island,"_",th,"th.tif")),
             datatype = "INT1U",
             overwrite = TRUE )
    endCluster()
    rm(fc2000)
    
    print(paste0("complete ", "temp_data/processed_lu/gfc_fc2000_outside_ip_",island,"_",th,"th.tif"))
    
    th <- th + 30
  }
}
### Execute the function
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/gfc_fc2000_outside_ip_",island,"_90th.tif")))){
    
    prepare_fc2000(island)
  }
}


##### PRIMARY FOREST 2000 #####

### Prepare intact, degraded and total primary forest cover in 2000
prepare_pfc2000 <- function(island){
  pfc2000 <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif")))
  
  beginCluster()
  
  # intact primary forest
  m_intact <- c(0,0,0,
                0,1,1, 
                1,2,0)
  rclmat_intact <- matrix(m_intact, ncol = 3, byrow = TRUE)
  
  raster::clusterR(pfc2000,
                   fun = raster::reclassify,
                   args = list(rcl = rclmat_intact, right = TRUE, include.lowest = TRUE), 
                   filename = file.path(paste0("temp_data/processed_lu/margono_intact_primary_forest_",island,"_aligned.tif")), 
                   datatype = "INT1U",
                   progress = "text",
                   overwrite = TRUE)
  
  # degraded primary forest 
  m_degraded <- c(0,0,0,
                  0,1,0, 
                  1,2,1)
  rclmat_degraded <- matrix(m_degraded, ncol = 3, byrow = TRUE)
  
  raster::clusterR(pfc2000,
                   fun = raster::reclassify,
                   args = list(rcl = rclmat_degraded, right = TRUE, include.lowest = TRUE), 
                   filename = file.path(paste0("temp_data/processed_lu/margono_degraded_primary_forest_",island,"_aligned.tif")), 
                   datatype = "INT1U",
                   progress = "text",
                   overwrite = TRUE)
  
  # total primary forest
  m_total <- c(0,0,0,
               0,1,1, 
               1,2,1)
  rclmat_total <- matrix(m_total, ncol = 3, byrow = TRUE)
  
  raster::clusterR(pfc2000,
                   fun = raster::reclassify,
                   args = list(rcl = rclmat_total, right = TRUE, include.lowest = TRUE), 
                   filename = file.path(paste0("temp_data/processed_lu/margono_total_primary_forest_",island,"_aligned.tif")), 
                   datatype = "INT1U",
                   progress = "text",
                   overwrite = TRUE)
  
  endCluster()
}
### Execute it 
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(island in IslandS){
  prepare_pfc2000(island) 
}  



##### MAKE BASELINE FOREST EXTENT VARIABLES AT PARCEL LEVEL ##### 

# We add a 2000 forest extent variable to be able to then restrict the sample to parcels that have a minimum of 
# forest extent in 2000 (e.g. exclude cities and hence remove parcels that have 0 tree cover in 2000). 
# We do this for two kinds of 2000 forest extents: 30% tree canopy density outside industrial plantations
# and 30% tree canopy density in primary forest. 
# The spatial procedure is exactly the same as in prepare_lucfp codes, because: 
# 1. code is already written, hence not effort intensive eventhough shorter code could be possible from there
# 2. This enables to compute grid cells with exact same area as those derived from raster and 
# used in final analysis, i.e. 3002.4m x 3008.4m rather than sharp 3km x 3km grid cells computed like:
# parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",PS/1000,"km_",CR/1000,"km_",sample,"_CR.rds")))
# parcels <- parcels[!duplicated(parcels$parcel_id),]
# parcels <- st_as_sf(parcels, coords = c("lon", "lat"), crs = indonesian_crs)   
# parcel_points <- st_geometry(parcels)
# parcel_buffers <- st_buffer(parcel_points, dist = PS/2)
# parcel_squares <- sapply(parcel_buffers, FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
# st_area(parcel_squares[1])
# res(parcel_raster)[1]*res(parcel_raster)[2]

# So, as in prepare_lucpfip, aggregate high resolution raster to parcel size, 
# and then convert to data frame while selecting parcels that are within catchment radius. 

#### AGGREGATE THE PIXELS TO A GIVEN PARCEL SIZE. ####

aggregate_baseline_fc <- function(island, parcel_size){
   
  # what we aggregate
  fc2000 <- raster(file.path(paste0("temp_data/processed_lu/gfc_fc2000_outside_ip_",island,"_30th.tif")))
  pfc2000 <- raster(file.path(paste0("temp_data/processed_lu/margono_total_primary_forest_",island,"_aligned.tif")))
  
  # define output file name
  output_filename_fc2000 <- file.path(paste0("temp_data/processed_lu/parcel_fc2000_",island,"_30th_",parcel_size/1000,"km.tif"))
  output_filename_pfc2000 <- file.path(paste0("temp_data/processed_lu/parcel_pfc2000_",island,"_total_",parcel_size/1000,"km.tif"))
  
  # aggregate fc2000 from the ~30m cells to parcel_size cells with sum function.
  raster::aggregate(fc2000, fact = c(parcel_size/res(fc2000)[1], parcel_size/res(fc2000)[2]),
                    expand = FALSE,
                    fun = sum,
                    na.rm = FALSE, # NA cells are in margins, see the NOTES part of prepare_lucpfip.R. 
                    # If FALSE, aggregations at margins that use NA are discarded because the sum would 
                    # be spurious as it would count all NA as 0s while it is not necessary the case.
                    filename = output_filename_fc2000,
                    datatype = "INT4U", # because the sum may go up to ~ 10 000 with parcel_size = 3000,
                    # but to more than 65k with parcel_size = 10000 so INT4U will be necessary;
                    overwrite = TRUE, 
                    progess = "text")
  
  raster::aggregate(pfc2000, fact = c(parcel_size/res(pfc2000)[1], parcel_size/res(pfc2000)[2]),
                    expand = FALSE,
                    fun = sum,
                    na.rm = FALSE, # NA cells are in margins, see the NOTES part of prepare_lucpfip.R.
                    # If FALSE, aggregations at margins that use NA are discarded because the sum would
                    # be spurious as it would count all NA as 0s while it is not necessary the case.
                    filename = output_filename_pfc2000,
                    datatype = "INT4U", # because the sum may go up to ~ 10 000 with parcel_size = 3000,
                    # but to more than 65k with parcel_size = 10000 so INT4U will be necessary;
                    overwrite = TRUE,
                    progess = "text")
}

PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
    aggregate_baseline_fc(island = Island,
                          parcel_size = PS)
  
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### CONVERT TO A DATAFRAME THE PARCELS WITHIN A GIVEN CATCHMENT RADIUS  ####

# Note on the use of mask and then st_is_within_distance:
# Parcels that are near the coast and hence consist of less ground than others are not discarded by the use 
# of island polygons when making the mask. This is just a particular fixed feature.  
# Mills are selected, as belonging to the island, and then the mask is made by extending large (60km) CR around them. 
# Then, parcels that remain (i.e. are not masked) are not discarded by st_is_within_distance,
# even if parts of their areas are in the sea, as long as their centroids are less than a given distance. 


to_panel_within_CR <- function(island, parcel_size, catchment_radius){
  
  ### Function description
  # This repeats roughly twice the same actions, once for parcels within catchment radiuses of IBS mills only, 
  # and once for all UML mills. 
  
  # raster_to_df converts the raster bricks of annual layers of parcels to a panel dataframe.
  # This is executed for each pf_type, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment radius.
  # 3. reshaping the values in these parcels to a long format panel dataframe
    
    ### IBS
    input_filename_fc2000 <- file.path(paste0("temp_data/processed_lu/parcel_fc2000_",island,"_30th_",parcel_size/1000,"km.tif"))
    input_filename_pfc2000 <- file.path(paste0("temp_data/processed_lu/parcel_pfc2000_",island,"_total_",parcel_size/1000,"km.tif"))
    
    fc2000 <- raster(input_filename_fc2000)
    pfc2000 <- raster(input_filename_pfc2000)
    ## 1. Masking.
    # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
    
    # Make the mask
    mills <- read.dta13(file.path("temp_data/processed_mill_geolocalization/IBS_UML_panel.dta"))
    # keep only a cross section of those that are geolocalized mills, on island of interest
    mills <- mills[mills$analysis_sample==1,]
    mills <- mills[!duplicated(mills$firm_id),]
    mills <- mills[mills$island_name == island,]
    #turn into an sf object.
    mills <- st_as_sf(mills,	coords	=	c("lon",	"lat"), crs=4326)
    # keep only the geometry, we do not need mills attributes here.
    mills <- st_geometry(mills)
    # set CRS and project
    mills_prj <- st_transform(mills, crs = indonesian_crs)
    #define big catchment areas to have a large AOI.
    mills_ca <- st_buffer(mills_prj, dist = 60000)
    # work with squares rather than with circles
    for(i in 1:length(mills_ca)){
      mills_ca[i] <- st_as_sfc(st_bbox(mills_ca[i]))
    }
    total_ca <- st_union(st_geometry(mills_ca))
    # coerce to a SpatialPolygon
    total_ca_sp <- as(total_ca, "Spatial")
    # keep mills_prj we need it below
    rm(total_ca, mills_ca, mills)
    
    # MASK
    mask(x = fc2000, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_IBS_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(fc2000)
    
    mask(x = pfc2000, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_IBS_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(pfc2000)
    
    
    ## 2. Selecting parcels within a given distance to a mill at least one year
    # (i.e. the parcel is present in the dataframe in all years even if it is within say 50km of a mill only since 2014)
    
    # Turn the masked rasters to a sf dataframe
    convert_to_df <- function(in_fc_type_file_name){
        
      m.parcels <- raster(file.path(paste0("temp_data/processed_lu/",in_fc_type_file_name,".tif")))
  
      m.df <- raster::as.data.frame(m.parcels, na.rm = TRUE, xy = TRUE, centroids = TRUE)
      
      m.df <- m.df %>% dplyr::rename(lon = x, lat = y)
      m.df <- st_as_sf(m.df, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
      
      # Remove here parcels that are not within the catchment area of a given size (defined by catchment radius)
      # coordinates of all mills (crs is indonesian crs, unit is meter)
      within <- st_is_within_distance(m.df, mills_prj, dist = catchment_radius)
      m.df <- m.df %>% dplyr::filter(lengths(within) >0)
      m.df <- m.df %>% st_drop_geometry()
      
      rm(within, m.parcels)

      ## 3. Reshaping to long format
      # make parcel id
      island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
      m.df$parcel_id <- paste0(island_id, c(1:nrow(m.df))) %>% as.numeric()

      return(m.df)
    }
    
  # execute it 
    # for fc2000
    in_fc_type_file_name <- paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_IBS_masked")

    m.df <- convert_to_df(in_fc_type_file_name)

    # rename the pixel count variable
    names(m.df)[names(m.df)==in_fc_type_file_name] <- "fc2000_30th_pixelcount"

    out_fc_type_file_name <- paste0("fc2000_cs_",
                                    island,"_30th_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_IBS_CR.tif")
    
    saveRDS(m.df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    rm(m.df)
    
    # for pfc2000
    in_fc_type_file_name <- paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_IBS_masked")
    
    m.df <- convert_to_df(in_fc_type_file_name)
    
    # rename the pixel count variable
    names(m.df)[names(m.df)==in_fc_type_file_name] <- "pfc2000_total_pixelcount"
    
    out_fc_type_file_name <- paste0("pfc2000_cs_",
                                    island,"_total_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_IBS_CR.tif")
    
    saveRDS(m.df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    rm(m.df)
    
    
    ### UML 
    input_filename_fc2000 <- file.path(paste0("temp_data/processed_lu/parcel_fc2000_",island,"_30th_",parcel_size/1000,"km.tif"))
    input_filename_pfc2000 <- file.path(paste0("temp_data/processed_lu/parcel_pfc2000_",island,"_total_",parcel_size/1000,"km.tif"))
    
    fc2000 <- raster(input_filename_fc2000)
    pfc2000 <- raster(input_filename_pfc2000)
    
    ## 1. Masking.
    # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
    
    # Make the mask
    mills <- read_xlsx(file.path("input_data/uml/mills_20200129.xlsx"))
    mills$latitude <- as.numeric(mills$latitude)
    mills$longitude <- as.numeric(mills$longitude)
    mills$lat <- mills$latitude
    mills$lon <- mills$longitude
    mills <- st_as_sf(mills,	coords	=	c("longitude",	"latitude"), crs = 4326)
    mills <- st_geometry(mills)
    mills_prj <- st_transform(mills, crs = indonesian_crs)
    
    # there is no island column in this dataset, hence we select mills on the specific island geographically
    select_within <- st_within(x = mills_prj, y = island_sf_prj[island_sf_prj$shape_des == island,])
    mills_prj <- mills_prj[lengths(select_within)>0,]
    
    #define big catchment areas to have a large AOI.
    mills_ca <- st_buffer(mills_prj, dist = 60000)
    # work with squares rather than with circles
    for(i in 1:length(mills_ca)){
      mills_ca[i] <- st_as_sfc(st_bbox(mills_ca[i]))
    }
    total_ca <- st_union(st_geometry(mills_ca))
    # coerce to a SpatialPolygon
    total_ca_sp <- as(total_ca, "Spatial")
    # keep mills_prj we need it below
    rm(total_ca, mills_ca, mills)
    
    # MASK
    mask(x = fc2000, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_UML_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(fc2000)
    
    mask(x = pfc2000, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_UML_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(pfc2000)
    
    ## 2. Selecting parcels within a given distance to a mill at least one year
    # (i.e. the parcel is present in the dataframe in all years even if it is within say 50km of a mill only since 2014)
    
    # Turn the masked rasters to a sf dataframe
    convert_to_df <- function(in_fc_type_file_name){
      
      m.parcels <- raster(file.path(paste0("temp_data/processed_lu/",in_fc_type_file_name,".tif")))
      
      m.df <- raster::as.data.frame(m.parcels, na.rm = TRUE, xy = TRUE, centroids = TRUE)
      
      m.df <- m.df %>% dplyr::rename(lon = x, lat = y)
      m.df <- st_as_sf(m.df, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
      
      # Remove here parcels that are not within the catchment area of a given size (defined by catchment radius)
      # coordinates of all mills (crs is indonesian crs, unit is meter)
      within <- st_is_within_distance(m.df, mills_prj, dist = catchment_radius)
      m.df <- m.df %>% dplyr::filter(lengths(within) >0)
      m.df <- m.df %>% st_drop_geometry()
      
      rm(within, m.parcels)
      
      ## 3. Reshaping to long format
      # make parcel id
      island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
      m.df$parcel_id <- paste0(island_id, c(1:nrow(m.df))) %>% as.numeric()
      
      return(m.df)
    }
    
    # execute it 
    # for fc2000
    in_fc_type_file_name <- paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_UML_masked")
    
    m.df <- convert_to_df(in_fc_type_file_name)
    
    # rename the pixel count variable
    names(m.df)[names(m.df)==in_fc_type_file_name] <- "fc2000_30th_pixelcount"
    
    out_fc_type_file_name <- paste0("fc2000_cs_",
                                    island,"_30th_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_UML_CR.tif")
    
    saveRDS(m.df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    
    # for pfc2000
    in_fc_type_file_name <- paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_UML_masked")
    
    m.df <- convert_to_df(in_fc_type_file_name)
    
    # rename the pixel count variable
    names(m.df)[names(m.df)==in_fc_type_file_name] <- "pfc2000_total_pixelcount"
    
    out_fc_type_file_name <- paste0("pfc2000_cs_",
                                    island,"_total_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_UML_CR.tif")
    
    saveRDS(m.df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
}

PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
    to_panel_within_CR(island = Island,
                       parcel_size = PS,
                       catchment_radius = CR)
    
    CR <- CR + 20000
  }
}




#### CONVERT TO A DATAFRAME THE PARCELS WITHIN A GIVEN CATCHMENT AREA  ####

# Note on the use of mask and then st_is_within_distance:
# Parcels that are near the coast and hence consist of less ground than others are not discarded by the use 
# of island polygons when making the mask. This is just a particular fixed feature.  
# Mills are selected, as belonging to the island, and then the mask is made by extending large (60km) CR around them. 
# Then, parcels that remain (i.e. are not masked) are not discarded by st_is_within_distance,
# even if parts of their areas are in the sea, as long as their centroids are less than a given distance. 


to_panel_within_CA <- function(island, parcel_size){
  
  ### Function description
  # This repeats roughly twice the same actions, once for parcels within catchment radiuses of IBS mills only, 
  # and once for all UML mills. 
  
  # raster_to_df converts the raster bricks of annual layers of parcels to a panel dataframe.
  # This is executed for each pf_type, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment area
  # 3. reshaping the values in these parcels to a long format panel dataframe
  
  ### IBS
  input_filename_fc2000 <- file.path(paste0("temp_data/processed_lu/parcel_fc2000_",island,"_30th_",parcel_size/1000,"km.tif"))
  input_filename_pfc2000 <- file.path(paste0("temp_data/processed_lu/parcel_pfc2000_",island,"_total_",parcel_size/1000,"km.tif"))
  
  fc2000 <- raster(input_filename_fc2000)
  pfc2000 <- raster(input_filename_pfc2000)
  ## 1. Masking.

  # Make the mask
  ibs <- read.dta13(file.path("temp_data/processed_mill_geolocalization/IBS_UML_panel.dta"))
  # keep only a cross section of those that are geolocalized mills, on island of interest
  ibs <- ibs[ibs$analysis_sample==1,]
  ibs <- ibs[!duplicated(ibs$firm_id),]
  ibs <- ibs[ibs$island_name == island,]
  #turn into an sf object.
  ibs <- st_as_sf(ibs,	coords	=	c("lon",	"lat"), crs=4326)
  # keep only the geometry, we do not need mills attributes here.
  ibs <- st_geometry(ibs)
  # set CRS and project
  ibs_prj <- st_transform(ibs, crs = indonesian_crs)
  #define big catchment areas to have a large AOI.
  ibs_ca <- st_buffer(ibs_prj, dist = 80000)
  # work with squares rather than with circles
  for(i in 1:length(ibs_ca)){
    ibs_ca[i] <- st_as_sfc(st_bbox(ibs_ca[i]))
  }
  total_ibs_ca <- st_union(st_geometry(ibs_ca))
  # coerce to a SpatialPolygon
  total_ibs_ca_sp <- as(total_ibs_ca, "Spatial")
  # keep ibs_prj we need it below
  rm(total_ibs_ca, ibs_ca, ibs)
  
  # MASK
  mask(x = fc2000, mask = total_ibs_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_IBS_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(fc2000)
  
  mask(x = pfc2000, mask = total_ibs_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_IBS_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(pfc2000)
  
  
  ## 2. Selecting parcels within a given distance to a mill at least one year

  # read in the duration matrix, common to all parcel layers (outcomes) 
  dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_IBS")))
  dur_mat <- dur_mat$durations
  
  # Turn the masked rasters to a sf dataframe
  convert_to_df <- function(in_fc_type_file_name, travel_time){
    
    m.parcels <- raster(file.path(paste0("temp_data/processed_lu/",in_fc_type_file_name,".tif")))
    
    ibs_msk_df <- raster::as.data.frame(m.parcels, na.rm = TRUE, xy = TRUE, centroids = TRUE)
    
    # the whole point of the below is to merge more safely the duration matrix and the parcel data frame. 
    ibs_msk_df <- ibs_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
    ibs_msk_df <- st_as_sf(ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
    ibs_msk_df <- st_transform(ibs_msk_df, crs = 4326)
    ibs_msk_df$lon <- st_coordinates(ibs_msk_df)[,"X"]
    ibs_msk_df$lat <- st_coordinates(ibs_msk_df)[,"Y"]
    ibs_msk_df <- mutate(ibs_msk_df, lonlat = paste0(lon, lat))
    ibs_msk_df <- st_drop_geometry(ibs_msk_df)
    
    # Besides, make IDs 
    island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
    ibs_msk_df$parcel_id <- paste0(island_id, c(1:nrow(ibs_msk_df))) %>% as.numeric()
    

    ## keep only the parcels that satisfy the travel time condition for at least one mill in the whole period.
    dur_mat_log <- dur_mat/(60) < travel_time
    dur_mat_log <- as.data.frame(dur_mat_log)
    
    # The row.names have been computed the same way the ibs_msk_df$lonlat column was is computed above. 
    colnames(dur_mat_log) <- paste0("firm_id_",colnames(dur_mat_log))
    dur_mat_log$lonlat <- row.names(dur_mat_log)
    
    ibs_msk_TT_df <- merge(ibs_msk_df, dur_mat_log, by = "lonlat", all = FALSE)
    ibs_msk_TT_df <- ibs_msk_TT_df[base::rowSums(ibs_msk_TT_df[,grepl("firm_id",colnames(ibs_msk_TT_df))], na.rm = TRUE)>0,]
    # na.rm = TRUE is in case of mills for which no duration could be computed. 
    # which is the case for 4 mills in Sumatra, 3 of which are not UML matched but desa centroid located.
    
    ibs_msk_TT_df <- ibs_msk_TT_df[,!grepl("firm_id",colnames(ibs_msk_TT_df))]
    ibs_msk_TT_df <- dplyr::select(ibs_msk_TT_df, -lonlat)

    rm(m.parcels, ibs_msk_df, dur_mat_log)

    return(ibs_msk_TT_df)
  }
  
  # execute it 
  # for fc2000
  in_fc_type_file_name <- paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_IBS_masked")
  
  for(travel_time in c(2,4,6)){
    ibs_msk_df <- convert_to_df(in_fc_type_file_name = in_fc_type_file_name, 
                          travel_time = travel_time)
    
    # rename the pixel count variable
    names(ibs_msk_df)[names(ibs_msk_df)==in_fc_type_file_name] <- "fc2000_30th_pixelcount"
    
    out_fc_type_file_name <- paste0("fc2000_cs_",
                                    island,"_30th_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_IBS_CA.tif")
    
    saveRDS(ibs_msk_df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    rm(ibs_msk_df)
    
    # for pfc2000
    in_fc_type_file_name <- paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_IBS_masked")
    
    ibs_msk_df <- convert_to_df(in_fc_type_file_name = in_fc_type_file_name, 
                          travel_time = travel_time)
    
    # rename the pixel count variable
    names(ibs_msk_df)[names(ibs_msk_df)==in_fc_type_file_name] <- "pfc2000_total_pixelcount"
    
    out_fc_type_file_name <- paste0("pfc2000_cs_",
                                    island,"_total_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_IBS_CA.tif")
    
    saveRDS(ibs_msk_df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    rm(ibs_msk_df)
  }
  
  
  ### UML 
  input_filename_fc2000 <- file.path(paste0("temp_data/processed_lu/parcel_fc2000_",island,"_30th_",parcel_size/1000,"km.tif"))
  input_filename_pfc2000 <- file.path(paste0("temp_data/processed_lu/parcel_pfc2000_",island,"_total_",parcel_size/1000,"km.tif"))
  
  fc2000 <- raster(input_filename_fc2000)
  pfc2000 <- raster(input_filename_pfc2000)
  
  ## 1. Masking.

  # Make the mask
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
  
  
  # MASK
  mask(x = fc2000, mask = total_uml_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_UML_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(fc2000)
  
  mask(x = pfc2000, mask = total_uml_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_UML_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(pfc2000)
  
  ## 2. Selecting parcels within a given distance to a mill at least one year
  
  # read in the duration matrix, common to all parcel layers (outcomes) 
  dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_UML")))
  dur_mat <- dur_mat$durations
  
  # Turn the masked rasters to a sf dataframe
  convert_to_df <- function(in_fc_type_file_name, travel_time){
    
    m.parcels <- raster(file.path(paste0("temp_data/processed_lu/",in_fc_type_file_name,".tif")))
    
    uml_msk_df <- raster::as.data.frame(m.parcels, na.rm = TRUE, xy = TRUE, centroids = TRUE)
    
    # the whole point of the below is to merge more safely the duration matrix and the parcel data frame. 
    uml_msk_df <- uml_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
    uml_msk_df <- st_as_sf(uml_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
    uml_msk_df <- st_transform(uml_msk_df, crs = 4326)
    uml_msk_df$lon <- st_coordinates(uml_msk_df)[,"X"]
    uml_msk_df$lat <- st_coordinates(uml_msk_df)[,"Y"]
    uml_msk_df <- mutate(uml_msk_df, lonlat = paste0(lon, lat))
    uml_msk_df <- st_drop_geometry(uml_msk_df)
    
    # Besides, make IDs 
    island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
    uml_msk_df$parcel_id <- paste0(island_id, c(1:nrow(uml_msk_df))) %>% as.numeric()
    
    
    ## keep only the parcels that satisfy the travel time condition for at least one mill in the whole period.
    dur_mat_log <- dur_mat/(60) < travel_time
    dur_mat_log <- as.data.frame(dur_mat_log)
    
    # The row.names have been computed the same way the uml_msk_df$lonlat column was is computed above. 
    colnames(dur_mat_log) <- paste0("firm_id_",colnames(dur_mat_log))
    dur_mat_log$lonlat <- row.names(dur_mat_log)
    
    uml_msk_TT_df <- merge(uml_msk_df, dur_mat_log, by = "lonlat", all = FALSE)
    uml_msk_TT_df <- uml_msk_TT_df[base::rowSums(uml_msk_TT_df[,grepl("firm_id",colnames(uml_msk_TT_df))], na.rm = TRUE)>0,]
    # na.rm = TRUE is in case of mills for which no duration could be computed. 
    # which is the case for 4 mills in Sumatra, 3 of which are not UML matched but desa centroid located.
    
    uml_msk_TT_df <- uml_msk_TT_df[,!grepl("firm_id",colnames(uml_msk_TT_df))]
    uml_msk_TT_df <- dplyr::select(uml_msk_TT_df, -lonlat)
    
    rm(m.parcels, uml_msk_df, dur_mat_log)
    
    return(uml_msk_TT_df)
  }
  
  # execute it 
  # for fc2000
  in_fc_type_file_name <- paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_UML_masked")
  
  for(travel_time in c(2,4,6)){
    uml_msk_df <- convert_to_df(in_fc_type_file_name = in_fc_type_file_name, 
                          travel_time = travel_time)
    
    # rename the pixel count variable
    names(uml_msk_df)[names(uml_msk_df)==in_fc_type_file_name] <- "fc2000_30th_pixelcount"
    
    out_fc_type_file_name <- paste0("fc2000_cs_",
                                    island,"_30th_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_UML_CA.tif")
    
    saveRDS(uml_msk_df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    rm(uml_msk_df)
    
    # for pfc2000
    in_fc_type_file_name <- paste0("parcel_pfc2000_",island,"_total_",parcel_size/1000,"km_UML_masked")
    
    uml_msk_df <- convert_to_df(in_fc_type_file_name = in_fc_type_file_name, 
                          travel_time = travel_time)
    
    # rename the pixel count variable
    names(uml_msk_df)[names(uml_msk_df)==in_fc_type_file_name] <- "pfc2000_total_pixelcount"
    
    out_fc_type_file_name <- paste0("pfc2000_cs_",
                                    island,"_total_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_UML_CA.tif")
    
    saveRDS(uml_msk_df, file = file.path(paste0("temp_data/processed_lu/",out_fc_type_file_name)))
    rm(uml_msk_df)
  }
}

PS <- 3000
IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
for(Island in IslandS){
  to_panel_within_CA(island = Island,
                     parcel_size = PS)
    
}
rm(to_panel_within_CA)




#### Gather the fc2000 variables for each parcel_size and catchment_radius combinations. ####

PS <- 3000  
sampleS <- c("IBS", "UML")
for(sample in sampleS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
    
    # For each Island, join columns of baseline forest extent pixelcounts
    df_list <- list()
    IslandS <- c("Sumatra", "Kalimantan", "Papua")
    for(Island in IslandS){
      
      df_fc2000   <- readRDS(file.path(paste0("temp_data/processed_lu/fc2000_cs_",
                                              Island,"_30th_",
                                              PS/1000,"km_",
                                              CR/1000,"km_",
                                              sample,"_CR.tif")))
      
      df_pfc2000 <- readRDS(file.path(paste0("temp_data/processed_lu/pfc2000_cs_",
                                             Island,"_total_",
                                             PS/1000,"km_",
                                             CR/1000,"km_",
                                             sample,"_CR.tif")))
      
      df_pfc2000 <- dplyr::select(df_pfc2000, -lon, -lat)
      df_list[[match(Island, IslandS)]] <- inner_join(df_fc2000, df_pfc2000, by = "parcel_id")
    }
    
    # stack the three Islands together
    indo_df <- bind_rows(df_list)
    
    ### Add columns of converted pixel counts to hectares.
    pixel_area <- (27.8*27.6)/(1e4)
    # fc2000
    indo_df <- mutate(indo_df, fc2000_30th_ha = fc2000_30th_pixelcount*pixel_area) 
    # pfc2000
    indo_df <- mutate(indo_df, pfc2000_total_ha = pfc2000_total_pixelcount*pixel_area) 
    
    indo_df <- dplyr::select(indo_df, 
                             parcel_id,
                             fc2000_30th_ha,
                             pfc2000_total_ha, 
                             fc2000_30th_pixelcount,
                             pfc2000_total_pixelcount,
                             everything())
    
    
    ### Add categorical variables to distinguish parcels that had a positive forest extent in 2000. 
    indo_df$any_fc2000_30th <- indo_df$fc2000_30th_ha > 0
    indo_df$any_pfc2000_total <- indo_df$pfc2000_total_ha > 0
    
    
    saveRDS(indo_df, file.path(paste0("temp_data/processed_parcels/baseline_fc_cs_",PS/1000,"km_",CR/1000,"km_",sample,"_CR.rds")))
    
    rm(indo_df, df_list)
    CR <- CR + 20000
  }
}







#### Gather the fc2000 variables for each parcel_size and catchment_area combinations. ####

PS <- 3000  
sampleS <- c("IBS", "UML")
for(sample in sampleS){
  for(TT in c(2,4,6)){
    
    # For each Island, join columns of baseline forest extent pixelcounts
    df_list <- list()
    IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
    for(Island in IslandS){
      
      df_fc2000   <- readRDS(file.path(paste0("temp_data/processed_lu/fc2000_cs_",
                                              Island,"_30th_",
                                              PS/1000,"km_",
                                              TT,"h_",
                                              sample,"_CA.tif"))) 
      
      df_pfc2000 <- readRDS(file.path(paste0("temp_data/processed_lu/pfc2000_cs_",
                                             Island,"_total_",
                                             PS/1000,"km_",
                                             TT,"h_",
                                             sample,"_CA.tif"))) 
      
      
      df_pfc2000 <- dplyr::select(df_pfc2000, -lon, -lat)
      df_list[[match(Island, IslandS)]] <- inner_join(df_fc2000, df_pfc2000, by = "parcel_id")
    }
    
    # stack the three Islands together
    indo_df <- bind_rows(df_list)
    
    ### Add columns of converted pixel counts to hectares.
    pixel_area <- (27.8*27.6)/(1e4)
    # fc2000
    indo_df <- mutate(indo_df, fc2000_30th_ha = fc2000_30th_pixelcount*pixel_area) 
    # pfc2000
    indo_df <- mutate(indo_df, pfc2000_total_ha = pfc2000_total_pixelcount*pixel_area) 
    
    indo_df <- dplyr::select(indo_df, 
                             parcel_id,
                             fc2000_30th_ha,
                             pfc2000_total_ha, 
                             fc2000_30th_pixelcount,
                             pfc2000_total_pixelcount,
                             everything())
    
    
    ### Add categorical variables to distinguish parcels that had a positive forest extent in 2000. 
    indo_df$any_fc2000_30th <- indo_df$fc2000_30th_ha > 0
    indo_df$any_pfc2000_total <- indo_df$pfc2000_total_ha > 0
    
    
    saveRDS(indo_df, file.path(paste0("temp_data/processed_parcels/baseline_fc_cs_",PS/1000,"km_",TT,"h_",sample,"_CA.rds")))
    
    rm(indo_df, df_list)
  }
}





