### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# 
#       BUILDING PANEL DATAFRAMES OF LAND USE CHANGE FROM FOREST TO INDUSTRIAL PLANTATION (LUCFIP)
# 
#   Inputs:   * Island polygons 
#             ---> temp_data/processed_indonesia_spatial/island_sf 
# 
#             * GFC data downloaded, extracted, and thresholded in prepare_gfc.R 
#             ---> temp_data/processed_lu/gfc_data_Indonesia_30th.tif 
#             ---> temp_data/processed_lu/gfc_data_Indonesia_60th.tif 
#             ---> temp_data/processed_lu/gfc_data_Indonesia_90th.tif
#
#             * 2000 oil palm plantations (Austin et al. 2017) for Sumatra, Kalimantan and Papua - not already processed
#             ---> input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_2000_WGS1984.tif 
#             
#             * 2015 oil palm plantations (Austin et al. 2017) for Sumatra, Kalimantan and Papua - already processed in prepare_lucpfip.R
#             ---> temp_data/processed_lu/austin_ioppm_2015_Sumatra_aligned.tif
#             ---> temp_data/processed_lu/austin_ioppm_2015_Kalimantan_aligned.tif
#             ---> temp_data/processed_lu/austin_ioppm_2015_Papua_aligned.tif
#      
#             * Georeferenced mills (from georeferencing works)                                           
#             ---> IBS_UML_panel.dta
#
#   Main outputs:  panel dataframes of lucfip pixel-event count in parcels of a given size, from 2001 to 2018,  
#                 for the whole Indonesia (Sumatra, Kalimantan, Papua "row-binded"), 
#                 for 3 forest definitions (30, 60 and 90 percent canopy closure thresholds. 
#             
#             There is one such dataframe for each combination of parcel size (only 3x3km for now) and catchment radius (10, 30, 50km)
#             ---> lucfip_panel_3km_10CR.rds 
#             ---> lucfip_panel_3km_30CR.rds 
#             ---> lucfip_panel_3km_50CR.rds
# 
#   
#   Actions:  This script consists of mainly three functions.  
#             0. load needed packages; set working directory; set raster options; define the crs used throughout the script. 
#                             the chunksize and maxmemory raster options should be set accordingly with the machine used,  
#                             considering that this script executes parallel functions using parallel::detectCores() - 1
#                             It is recommended to not change default memory raster options. 
#
#             1. prepare_pixel_lucfip(island)
# 
#             2. aggregate_lucfip(island, parcel_size)
# 
#             3. to_panel_within_CR(island, parcel_size, catchment_radius)
#
#             Finally, functions are run and there outputs are merged across islands and primary forest types.   
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("data.table", "dplyr", "readstata13", "readxl", 
                   "raster", "rgdal", "sp", "sf",
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



### NEW FOLDERS USED IN THIS SCRIPT 

# annual layers will be stored in a subfolder just for the sake of tidyness in processed_lu
dir.create("temp_data/processed_lu/annual_maps")

# dataframes of parcels are stored here 
dir.create("temp_data/processed_parcels")

# raster temp files are stored there
dir.create("temp_data/raster_tmp")


### RASTER OPTIONS ### 
# Do not change chunksize, as it is not safe in combinatin with clusterR and machine specific. 
rasterOptions(timer = TRUE, 
              tmpdir = "temp_data/raster_tmp")

### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Applied to GFC data, this makes a resolution of 27.6 ; 27.8 meters.
# Applied to GFC data, the crs from Margono et al. 2014 
# (+proj=sinu +lon_0=140 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs )
# makes a resolution of 27.7 ; 27.8 meters. 


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# preliminary: read the island shapefile
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)


##### 1. PREPARE 30m PIXEL-LEVEL MAPS OF LUCFIP ##### 

prepare_pixel_lucfip <- function(island){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Define area of interest (AOI) ####
  
  aoi <- island_sf[island_sf$island_name == island,]
  aoi <- st_as_sfc(st_bbox(aoi))
  aoi_sp <- as(aoi, "Spatial")
  rm(aoi)
  
  # it is not projected # 
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Prepare loss layer from GFC data  ####
  th <- 30
  while(th < 100){
    # Import gfc data prepare (Downloaded, extracted, thresholded with gfcanalysis package in prepare_gfc.R)
    thed_gfc_data <- brick(file.path(paste0("temp_data/processed_lu/gfc_data_Indonesia_",th,"th.tif"))) 
    
    
    ### Select the loss layer (15 ad 40 are arbitrary, the max value of loss layer is an integer corresponding 
    # to the latest year after 2000 when loss is observed. We need a GFC version with this year being at least 2015. 
    # and we do not want to select a layer with percentage and values up to 100. 
    
    loss <- thed_gfc_data[[which(thed_gfc_data@data@max > 15 & thed_gfc_data@data@max < 40)]]
    
    rm(thed_gfc_data)
    
    ### Crop to extent of island 
    crop(loss, aoi_sp, 
         filename = file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_",th,"th.tif")),
         datatype = "INT1U",
         overwrite = TRUE)
    
    
    ### Project loss layer
    
    # This is necessary because we will need to make computations on this map within mills' catchment *areas*.
    # If one does not project maps, then catchment areas all have different areas while being defined with a common buffer.
    
    loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_",th,"th.tif")))
    
    beginCluster() # this uses by default detectCores() - 1 
    
    projectRaster(loss,
                  method = "ngb",
                  crs = indonesian_crs, 
                  filename = file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_",th,"th_prj.tif")), 
                  datatype = "INT1U",
                  overwrite = TRUE)
    
    endCluster()
    
    rm(loss) 
    removeTmpFiles(h=0)
    
    print(paste0("complete ", "temp_data/processed_lu/gfc_loss_",island,"_",th,"th_prj.tif"))
    
    th <- th + 30
  }
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Prepare plantation maps  (Austin et al. (2017))####
  
  ### Prepare 2000 industrial oil palm plantation maps, (2015 one has been prepared in prepare_lucpfip.R already). 
  ioppm2000 <- raster(file.path("input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_2000_WGS1984.tif")) 
  
  # As such, these raw data are lat-lon (WGS1984)
  # extent is Indonesia  
  # resolution is 0.002277, 0.002277
  # while gfc_data resolution is 0.00030, 0.00025
  # Cells have value 1 for oil palm plantation, NA else. 
  
  # first adjust ioppm2000 maps to the same extent as gfc_data, i.e. the island bounding box. 
  # intermediary step to reclassify NA to 0 (not time consuming)
  # then projectRaster to match crs res - disaggregate before is not necessary (yields the same result)
  
  
  ### Adjust to roughly the same extent as GFC loss, i.e. island bbox.  
  # First, crop from Indonesia wide to island extent (in the longitude mainly)
  # then, extend it (in the latitude mainly) because raw data do not cover northern part of Sumatra. 
  
  cropped_ioppm2000 <- raster::crop(ioppm2000, y = aoi_sp)
  
  extended_ioppm2000 <- raster::extend(cropped_ioppm2000, y = aoi_sp, value = NA) # NA is the default
  
  writeRaster(extended_ioppm2000,
              filename = file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,".tif")),
              datatype = "INT1U",
              overwrite = TRUE)
  
  
  rm(extended_ioppm2000, cropped_ioppm2000, ioppm2000)
  
  
  ### Reclassify NA into 0 
  ioppm2000 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,".tif")))
  
  raster::reclassify(ioppm2000, 
                     rcl = cbind(NA,0), 
                     filename = file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,"_reclassified.tif")), 
                     overwrite=TRUE, 
                     datatype = "INT1U")
  
  
  ### Align to GFC loss crs, resolution and exact extent.  
  
  # read the ioppm2000 map (with roughly the island extent and reclassified) 
  ioppm2000 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,"_reclassified.tif")))
  
  # Read the target GFC loss layer
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
  
  beginCluster() # this uses by default detectCores() - 1
  
  projectRaster(from = ioppm2000, to = loss,
                method = "ngb",
                filename = file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif")),
                datatype = "INT1U",
                overwrite = TRUE )
  endCluster()
  # ~4200s. 
  
  rm(ioppm2000, loss)
  removeTmpFiles(h=0)  
  
  print(paste0("complete ", "temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif"))
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Overlay forest loss and oil palm plantation maps ####
  
  # We want to keep forest loss pixels only within 2015 plantations in order to induce forest conversion to plantation,
  # BUT outside 2000 plantations, in order not to count plantation renewals as forest conversion to plantation.
  # po maps are binary with 1 meaning plantation in 2015 (or 2000 resp.))
  
  
  ioppm2000 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif")))
  ioppm2015 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,"_aligned.tif")))

  # primary forest
  pf <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif")))
  
  # overlay function
  overlay_maps <- function(rs){rs[[1]]*(1-rs[[2]])*rs[[3]]*is.na(rs[[4]])}
  # multiplies a cell of forest loss (rs[[1]]) by 0 (i.e. "removes" it) if it it is a plantation in 2000 (rs[[2]]) 
  # or if is not a plantation in 2015 (rs[[3]])
  # or if it is within primary forest (i.e. if pf is NA)
  
  ### For each threshold, overlay forest loss map with plantations maps in a clusterR setting 
  th <- 30
  while(th < 100){
    # call the loss layer for threshold th
    loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_",th,"th_prj.tif")))
    
    
    # stack loss with plantation maps (necessary for clusterR)
    rs <- stack(loss, ioppm2000, ioppm2015, pf)
    
    # run the computation in parallel with clusterR, as cells are processed one by one independently.
    beginCluster() # uses by default detectedCores() - 1
    clusterR(rs,
             fun = calc, # note we use calc but this is equivalent to using overlay 
             # (but more appropriate to the input being a stack)
             args = list(overlay_maps),
             filename = file.path(paste0("temp_data/processed_lu/lucfip_",island,"_",th,"th.tif")),
             datatype = "INT1U",
             overwrite = TRUE )
    endCluster()
    rm(loss)
    
    print(paste0("complete ", "temp_data/processed_lu/lucfip_",island,"_",th,"th.tif"))

    th <- th + 30
  }
  
  # ~ 4500 seconds / threshold
  rm(ioppm2000, ioppm2015, overlay_maps)
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Split lucfip to annual layers ####
  
  ### Function description
  # parallel_split has for input the single lucfip layer where each pixel has a value corresponding to the year when a lucfip event occured;
  # it outputs annual layers in each of which pixels are either 1 if a lucfip event occured that year, and 0 else.
  # the tasks are year specific and independent across years, therefore they are executed parallely over years.
  parallel_split <- function(th, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the threshold level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the threshold level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_split <- function(time){
      # define process
      process <- file.path(paste0("temp_data/processed_lu/lucfip_",island,"_",th,"th.tif"))
      # #set temp directory
      dir.create(paste0(process,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(process,"_Tmp")))
      # read in annual raster layer
      lucfip_prj <- raster(process)
      # split it into annual binary layers
      calc(lucfip_prj,
           fun = function(x){if_else(x == time, true = 1, false = 0)},
           filename = file.path(paste0("temp_data/processed_lu//annual_maps/lucfip_",island,"_",th,"th_", years[time],".tif")),
           datatype = "INT1U",
           overwrite = TRUE )
      # remove process temporary files
      unlink(file.path(paste0(process,"_Tmp")), recursive = TRUE)
    }
    
    ## register cluster
    registerDoParallel(cores = ncores)
    
    ## define foreach object.
    foreach(t = 1:length(years),
            # .combine combine the outputs as a mere character list (by default)
            .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
            .multicombine = TRUE,
            .export = c("island"),
            .packages = c("dplyr", "raster", "rgdal")
    ) %dopar%  annual_split(time = t)
  }
  
  ### Execute it for each forest definition
  th <- 30
  while(th < 100){
    parallel_split(th, detectCores() - 1) # ~500 seconds / annual layer
    
    removeTmpFiles(h=0)
    
    th <- th + 30
  }
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")  
  return(print(paste0("complete prepare_lucfip ", island)))
}




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 2. AGGREGATE THE PIXELS TO A GIVEN PARCEL SIZE. #####

aggregate_lucfip <- function(island, parcel_size){
  
  ### Function description
  # The function has for inputs annual layers of lucfip events at the pixel level.
  # It aggregates these pixels to a parcel size defined by parcel_size (in meters).
  # The aggregation operation is the sum of the pixel lucfip events.
  # Each annual aggregation is tasked in parallel.
  parallel_aggregate <- function(th, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the threshold level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the threshold level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_aggregate <- function(time){
      # Define which process (year and threshold) we are in:
      processname <- file.path(paste0("temp_data/processed_lu/annual_maps/lucfip_",island,"_",th,"th_", years[time],".tif"))
      #create process-specific temp directory
      dir.create(paste0(processname,"_Tmp"), showWarnings = FALSE)
      #set temp directory
      rasterOptions(tmpdir=file.path(paste0(processname,"_Tmp")))
      # read in the indonesia wide raster of lucfip at a given time and for a given threshold.
      annual_defo <- raster(processname)
      # aggregate it from the 30m cells to parcel_sizem cells with mean function.
      raster::aggregate(annual_defo, fact = c(parcel_size/res(annual_defo)[1], parcel_size/res(annual_defo)[2]),
                        expand = FALSE,
                        fun = sum,
                        na.rm = FALSE, # NA cells are in margins, see the NOTES part. If FALSE, aggregations at margins that use NA 
                        # are discarded because the sum would be spurious as it would count all NA as 0s while it is not necessary the case.
                        filename = paste0("temp_data/processed_lu/annual_maps/parcel_lucfip_",island,"_",parcel_size/1000,"km_",th,"th_",years[time],".tif"),
                        datatype = "INT4U", # because the sum may go up to ~ 10 000 with parcel_size = 3000,
                        # but to more than 65k with parcel_size = 10000 so INT4U will be necessary;
                        overwrite = TRUE)
      #removes entire temp directory without affecting other running processes (but there should be no temp file now)
      unlink(file.path(paste0(processname,"_Tmp")), recursive = TRUE)
      #unlink(file.path(tmpDir()), recursive = TRUE)
    }
    
    ## register cluster
    registerDoParallel(cores = ncores)
    
    ##  define foreach object.
    foreach(t = 1:length(years),
            # .combine combine the outputs as a mere character list (by default)
            .inorder = FALSE, # we don't care that the results be combine in the same order they were submitted
            .multicombine = TRUE,
            .export = c("island", "parcel_size"),
            .packages = c("raster", "rgdal")
    ) %dopar% annual_aggregate(time = t)
  }
  
  ### Execute the function to compute the RasterBrick object of 18 annual layers for each forest definition threshold
  th <- 30
  while(th < 100){
    # run the computation, that writes the layers 
    parallel_aggregate(th, detectCores() - 1)
    
    # brick the layers together and write the brick
    rasterlist <- list.files(path = "temp_data/processed_lu/annual_maps", 
                             pattern = paste0("parcel_lucfip_",island,"_",parcel_size/1000,"km_",th,"th_"), 
                             full.names = TRUE) %>% as.list()
    
    parcels_brick <- brick(rasterlist)
    
    writeRaster(parcels_brick,
                filename = file.path(paste0("temp_data/processed_lu/parcel_lucfip_",island,"_",parcel_size/1000,"km_",th,"th.tif")),
                datatype = "INT4U",
                overwrite = TRUE)
    
    rm(rasterlist, parcels_brick)
    removeTmpFiles(h=0)
    
    th <- th + 30
  }
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")    
  
  print(paste0("complete aggregate_lucfip ",island," ",parcel_size/1000, "km"))
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 3. CONVERT TO A DATAFRAME THE PARCELS WITHIN A GIVEN CATCHMENT RADIUS  #####

to_panel_within_CR <- function(island, parcel_size, catchment_radius){
  
  ### Function description
  # raster_to_df converts the raster bricks of annual layers of parcels to a panel dataframe.
  # This is executed for each threshold, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment radius.
  # 3. reshaping the values in these parcels to a long format panel dataframe
  raster_to_df <- function(th){
    
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## 1. Masking.
    # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
    
    ### IBS 
    
    ## 1. Masking.
    # Probably more efficient as the st_is_within does not need to be executed over all Indonesian cells but only those within the largest catchment_radius.
    
    # Make the mask
    mills <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
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
    
    # Mask
    parcels_brick_name <- paste0("parcel_lucfip_",island,"_",parcel_size/1000,"km_",th,"th")
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
    
    mask(x = parcels_brick, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(parcels_brick, total_ca_sp)
    
    
    ## 2. Selecting parcels within a given distance to a mill at least one year
    # (i.e. the parcel is present in the dataframe in all years even if it is within say 50km of a mill only since 2014)
    
    # Turn the masked raster to a sf dataframe
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))
    
    ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
    
    # the whole point of the below is to merge more safely the parcel data frames all together. 
    ibs_msk_df <- ibs_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
    ibs_msk_df <- st_as_sf(ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
    ibs_msk_df <- st_transform(ibs_msk_df, crs = 4326)
    ibs_msk_df$lon <- st_coordinates(ibs_msk_df)[,"X"]%>% round(6) # the rounding is bc otherwise there are very little differences in the decimals of the coordinates... 
    ibs_msk_df$lat <- st_coordinates(ibs_msk_df)[,"Y"]%>% round(6) 
    ibs_msk_df <- mutate(ibs_msk_df, lonlat = paste0(lon, lat))
    ibs_msk_df <- st_drop_geometry(ibs_msk_df)
    ibs_msk_df <- st_as_sf(ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
    
    # Remove here parcels that are not within the catchment area of a given size (defined by catchment radius)
    # coordinates of all mills (crs is indonesian crs, unit is meter)
    within <- st_is_within_distance(ibs_msk_df, mills_prj, dist = catchment_radius)
    ibs_msk_df <- ibs_msk_df %>% dplyr::filter(lengths(within) >0)
    ibs_msk_df <- ibs_msk_df %>% st_drop_geometry()
    
    rm(within, parcels_brick)
    
    
    ## 3. Reshaping to long format
    # # make parcel id
    # island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
    # ibs_msk_df$parcel_id <- paste0(island_id, c(1:nrow(ibs_msk_df))) %>% as.numeric()
    
    # vector of the names in the wide format of our time varying variables
    # the column names are the layer names in the parcels_brick + the layer index, separated by "." see examples in raster::as.data.frame
    varying_vars <- paste0(parcels_brick_name, "_IBS_masked.", seq(from = 1, to = 18))
    
    # reshape to long
    m.df <- stats::reshape(ibs_msk_df,
                           varying = varying_vars,
                           v.names = paste0("lucfip_pixelcount_",th,"th"),
                           sep = ".",
                           timevar = "year",
                           idvar = c("lonlat"), # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                           ids = "lonlat",
                           direction = "long",
                           new.row.names = seq(from = 1, to = nrow(ibs_msk_df)*length(years), by = 1))
    
    rm(varying_vars, ibs_msk_df)
    # replace the indices from the raster::as.data.frame with actual years.
    m.df <- mutate(m.df, year = years[year])
    
    m.df <- setorder(m.df, lonlat, year)
    saveRDS(m.df,
            file = file.path(paste0("temp_data/processed_parcels/lucfip_panel_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_IBS_CR_",
                                    th,"th.rds")))
    
    
    ### UML 
    
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
    
    # Mask
    parcels_brick_name <- paste0("parcel_lucfip_",island,"_",parcel_size/1000,"km_",th,"th")
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
    
    mask(x = parcels_brick, mask = total_ca_sp,
         filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")),
         datatype = "INT4U",
         overwrite = TRUE)
    
    rm(parcels_brick, total_ca_sp)
    
    
    ## 2. Selecting parcels within a given distance to a mill at least one year
    # (i.e. the parcel is present in the dataframe in all years even if it is within say 50km of a mill only since 2014)
    
    # Turn the masked raster to a sf dataframe
    parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_UML_masked.tif")))
    
    uml_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
    
    # the whole point of the below is to merge more safely the parcel data frames all together. 
    uml_msk_df <- uml_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
    uml_msk_df <- st_as_sf(uml_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
    uml_msk_df <- st_transform(uml_msk_df, crs = 4326)
    uml_msk_df$lon <- st_coordinates(uml_msk_df)[,"X"]%>% round(6) # the rounding is bc otherwise there are very little differences in the decimals of the coordinates... 
    uml_msk_df$lat <- st_coordinates(uml_msk_df)[,"Y"]%>% round(6) 
    uml_msk_df <- mutate(uml_msk_df, lonlat = paste0(lon, lat))
    uml_msk_df <- st_drop_geometry(uml_msk_df)
    uml_msk_df <- st_as_sf(uml_msk_df, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
    
    # Remove here parcels that are not within the catchment area of a given size (defined by catchment radius)
    # coordinates of all mills (crs is indonesian crs, unit is meter)
    within <- st_is_within_distance(uml_msk_df, mills_prj, dist = catchment_radius)
    uml_msk_df <- uml_msk_df %>% dplyr::filter(lengths(within) >0)
    uml_msk_df <- uml_msk_df %>% st_drop_geometry()
    
    rm(within, parcels_brick)
    
    
    ## 3. Reshaping to long format
    # # make parcel id
    # island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
    # uml_msk_df$parcel_id <- paste0(island_id, c(1:nrow(uml_msk_df))) %>% as.numeric()
    
    # vector of the names in the wide format of our time varying variables
    # the column names are the layer names in the parcels_brick + the layer index, separated by "." see examples in raster::as.data.frame
    varying_vars <- paste0(parcels_brick_name, "_UML_masked.", seq(from = 1, to = 18))
    
    # reshape to long
    m.df <- stats::reshape(uml_msk_df,
                           varying = varying_vars,
                           v.names = paste0("lucfip_pixelcount_",th,"th"),
                           sep = ".",
                           timevar = "year",
                           idvar = c("lonlat"), # don't put "lon" and "lat" in there, otherwise memory issue (see https://r.789695.n4.nabble.com/reshape-makes-R-run-out-of-memory-PR-14121-td955889.html)
                           ids = "lonlat",
                           direction = "long",
                           new.row.names = seq(from = 1, to = nrow(uml_msk_df)*length(years), by = 1))
    
    rm(varying_vars, uml_msk_df)
    # replace the indices from the raster::as.data.frame with actual years.
    m.df <- mutate(m.df, year = years[year])
    
    m.df <- setorder(m.df, lonlat, year)
    saveRDS(m.df,
            file = file.path(paste0("temp_data/processed_parcels/lucfip_panel_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_UML_CR_",
                                    th,"th.rds")))
    
  }
  
  ### Execute it
  th <- 30
  while(th < 100){
    
    raster_to_df(th)
  
    removeTmpFiles(h=0)
    
    th <- th + 30
  }  
  
  
      print(paste0("complete to_panel_within_CR ",island," ",parcel_size/1000,"km ",catchment_radius/1000,"CR"))

}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### EXECUTE FUNCTIONS AND MERGE THE OUTPUTS #####

#### Execute the functions ####
# Only if their outputs have not been already computed

### Prepare a 30m pixel map of lucfip for each Island
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/annual_maps/lucfip_",Island,"_total_2018.tif")))){
    
    prepare_pixel_lucfip(Island)
  }
}

### Aggregate this Island map to a chosen parcel size (3km, 6km and 9km for instance)
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/parcel_lucfip_",Island,"_",PS/1000,"km_total.tif")))){
    
    aggregate_lucfip(island = Island,
                      parcel_size = PS)
  }
}

### For that Island and for each aggregation factor, extract panels of parcels within different catchment area sizes 
# (radius of 10km, 30km and 50km)
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan", "Papua")
for(Island in IslandS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
    # on ne fait pas confiance à celui qui existe pour Sumatra car il vient d'une erreur dans aggregate_lucfip
    #if(!file.exists(file.path(paste0("temp_data/processed_parcels/lucfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"CR_total.rds")))){
    
    to_panel_within_CR(island = Island,
                       parcel_size = PS,
                       catchment_radius = CR)
    
    #} # only the function execution is conditioned to the file existance, not the loop incrementation
    
    CR <- CR + 20000
  }
}

#### Gather the lucfip variables for each parcel_size and catchment_radius combinations. ####
PS <- 3000  
sampleS <- c("IBS", "UML")
for(sample in sampleS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
    
    # For each Island, join columns of lucfip variable for different forest thresholds. 
    df_list <- list()
    IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
    for(Island in IslandS){
      
      df30 <- readRDS(paste0("temp_data/processed_parcels/lucfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_30th.rds"))
      df60 <- readRDS(paste0("temp_data/processed_parcels/lucfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_60th.rds"))
      df90 <- readRDS(paste0("temp_data/processed_parcels/lucfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_90th.rds"))
      
      df60 <- dplyr::select(df60, -lon, -lat, -idncrs_lon, -idncrs_lat)
      df <- inner_join(df30, df60, by = c("lonlat", "year"))

      df90 <- dplyr::select(df90, -lon, -lat, -idncrs_lon, -idncrs_lat)
      df_list[[match(Island, IslandS)]] <- inner_join(df, df90, by = c("lonlat", "year"))
      
      if(nrow(df_list[[match(Island, IslandS)]]) != nrow(df30)){stop("data frames do not all have the same set of grid cells")}
      rm(df, df30, df60, df90)
      
    }
    
    # stack the three Islands together
    indo_df <- bind_rows(df_list)
    
    ### Add columns of converted pixel counts to hectares.
    # pixel_area <- (27.8*27.6)/(1e4)
    # # 30th
    # indo_df <- mutate(indo_df, lucfip_ha_30th = lucfip_pixelcount_30th*pixel_area) 
    # # 60th
    # indo_df <- mutate(indo_df, lucfip_ha_60th = lucfip_pixelcount_60th*pixel_area) 
    # # 90th
    # indo_df <- mutate(indo_df, lucfip_ha_90th = lucfip_pixelcount_90th*pixel_area) 
    # 
    # indo_df <- dplyr::select(indo_df, lonlat, year, 
    #                          lucfip_ha_30th,
    #                          lucfip_ha_60th, 
    #                          lucfip_ha_90th,
    #                          lucfip_pixelcount_30th,
    #                          lucfip_pixelcount_60th, 
    #                          lucfip_pixelcount_90th,
    #                          everything())
    
    
    # ### Add categoric variable to distinguish parcels that experienced at least one lucpfip pixel event 
    # indo_df <- merge(indo_df,
    #                  ddply(indo_df, "lonlat", summarize, 
    #                        any_lucfip = sum(lucfip_ha_total, na.rm = TRUE) >0 ), 
    #                  by = "lonlat")
    
    
    saveRDS(indo_df, paste0("temp_data/processed_parcels/lucfip_panel_",PS/1000,"km_",CR/1000,"km_",sample,"_CR.rds"))
    
    rm(indo_df, df_list)
    
    CR <- CR + 20000
  }
}  
