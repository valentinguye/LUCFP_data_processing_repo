### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# 
#       BUILDING PANEL DATAFRAMES OF LAND USE CHANGE FROM PRIMARY FOREST TO INDUSTRIAL PLANTATION (LUCPFIP)
# 
#   Inputs:   * Island polygons 
#             ---> temp_data/processed_indonesia_spatial/island_sf 
#
#             * GFC data downloaded, extracted, and thresholded in prepare_gfc.R 
#             ---> temp_data/processed_lu/gfc_data_Indonesia_30th.tif
#
#             * 2015 oil palm plantations (Austin et al. 2017) for Sumatra, Kalimantan and Papua,
#             ---> oilpalm_2015_WGS1984.tif 
# 
#             * Primary forest extent in 2000 (Margono et al. 2014) 
#             ---> timeseq_change00_12.tif
#      
#             * Georeferenced mills (from georeferencing works)                                           
#             ---> IBS_UML_panel.dta
#
#   Main outputs:  panel dataframes of LUCPFIP pixel-event count in parcels of a given size, from 2001 to 2018,  
#                 for the whole Indonesia (Sumatra, Kalimantan, Papua "row-binded"), 
#                 for 3 forest definitions (30% canopy closure in intact, degraded, and intact or degraded (total) primary forest).
#             
#             There is one such dataframe for each combination of parcel size (only 3x3km for now) and catchment radius (10, 30, 50km)
#             ---> lucpfip_panel_3km_10CR.rds 
#             ---> lucpfip_panel_3km_30CR.rds 
#             ---> lucpfip_panel_3km_50CR.rds
# 
#   
#   Actions:  This script consists of mainly three functions.  
#             0. load needed packages; set working directory; set raster options; define the crs used throughout the script. 
#                             the chunksize and maxmemory raster options should be set accordingly with the machine used,  
#                             considering that this script executes parallel functions using parallel::detectCores() - 1
#                             It is recommended to not change default memory raster options. 
#
#             1. prepare_pixel_lucpfip(island)
# 
#             2. aggregate_lucpfip(island, parcel_size)
# 
#             3. to_panel_within_CR(island, parcel_size, catchment_radius)
#
#             Finally, functions are run and there outputs are merged across islands and primary forest types.   
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
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


# preliminary: read the island shapefile
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 1. PREPARE 30m PIXEL-LEVEL MAPS OF LUCPFIP ##### 
#island <- "Sumatra"

prepare_pixel_lucpfip_dynamics <- function(island){
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Define area of interest (AOI) ####
  aoi <- island_sf[island_sf$shape_des == island,]
  aoi <- st_as_sfc(st_bbox(aoi))
  aoi_sp <- as(aoi, "Spatial")
  rm(aoi)
  
  # it is not projected # 
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Prepare loss layer from GFC data  ####
  
  # Import gfc data prepared (Downloaded, extracted, thresholded with gfcanalysis package in prepare_gfc.R)
  thed_gfc_data <- brick(file.path(paste0("temp_data/processed_lu/gfc_data_Indonesia_30th.tif"))) 
  
  
  ### Select the loss layer (15 and 40 are arbitrary, the max value of loss layer is an integer corresponding 
  # to the latest year after 2000 when loss is observed. We need a GFC version with this year being at least 2015. 
  # and we do not want to select a layer with percentage and values up to 100. 
  
  loss <- thed_gfc_data[[which(thed_gfc_data@data@max > 15 & thed_gfc_data@data@max < 40)]]
  
  rm(thed_gfc_data)
  
  ### Crop to extent of island 
  crop(loss, aoi_sp, 
       filename = file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th.tif")),
       datatype = "INT1U",
       overwrite = TRUE)
  
  
  ### Project loss layer
  
  # This is necessary because we will need to make computations on this map within mills' catchment *areas*.
  # If one does not project maps, then catchment areas all have different areas while being defined with a common buffer.
  
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th.tif")))
  
  beginCluster() # this uses by default detectCores() - 1 
  
  projectRaster(loss,
                method = "ngb",
                crs = indonesian_crs, 
                filename = file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")), 
                datatype = "INT1U",
                overwrite = TRUE)
  
  endCluster()
  
  rm(loss) 
  removeTmpFiles(h=0)
  
  print(paste0("complete ", "temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif"))
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Make a layer of earliest detected industrial oil palm plantations ####
  
  # Read all layers of ioppm of interest
  ioppm_list <- list()
  decades <- c("2000", "2005", "2010", "2015")
  length(ioppm_list) <- length(decades)
  names(ioppm_list) <- decades
  
  for(decadal_year in decades){
    ioppm <- raster(file.path(paste0("input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/oilpalm_",decadal_year,"_WGS1984.tif")))
    
    # stabilizes the transformations implemented in this loop (otherwise 2005 turns NAs into 0s for instance...)
    # dataType(ioppm) # this changes to FLT4S with the operations below but it is not important. 
    # NAvalue(ioppm)
    # First, crop from Indonesia wide to island extent (in the longitude mainly)
    ioppm <- raster::crop(ioppm, y = aoi_sp)
    # then extend it (in the latitude mainly) because raw data do not cover northern part of Sumatra. 
    ioppm <- raster::extend(ioppm, y = aoi_sp, value = NA)
    
    # change binary occurrence of event into time of event
    ioppm_list[[decadal_year]] <- ioppm*(as.numeric(decadal_year)-2000+200) # (-1800 is a trick to enable dataType to pass under the bar of dataType from INT2 to INT1, without losing information)
  }
  
  rs <- raster::stack(ioppm_list[["2000"]], ioppm_list[["2005"]], ioppm_list[["2010"]], ioppm_list[["2015"]])
  
  # the min requires that cells that are never plantations are NA and not O (and hence the na.rm = TRUE) 
  earliest <- calc(rs, min, na.rm = TRUE)
  # dataType(earliest)
  # NAvalue(earliest)
  # getValues(earliest) %>% unique()  

  # lightening the dataType requires that there are no NAs
  earliest <- raster::reclassify(earliest, rcl = cbind(NA,0))
  # getValues(earliest) %>% unique()  
  # plot(earliest)
  
  rm(rs, ioppm_list)
  
  ### COMPUTE TIME LAPS BY OVERLAYING EARLIEST AND LOSS  
  
  ## PRELIMINARY ALIGNEMENT
  if(!file.exists(file.path(paste0("temp_data/processed_lu/austin_earliest_detected_",island,"_aligned.tif")))){
  # Read the target GFC loss layer
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
  
  # get the two rasters on the same projection
  beginCluster() # this uses by default detectCores() - 1
  projectRaster(from = earliest, to = loss,
                method = "ngb",
                filename = file.path(paste0("temp_data/processed_lu/austin_earliest_detected_",island,"_aligned.tif")),
                datatype = "INT1U", # INT2U will allow NA at the margins. 
                overwrite = TRUE)
  endCluster()
  print(paste0("completed austin_earliest_detected_",island,"_aligned.tif"))
  removeTmpFiles(h=0)  
  }
  
  ## OVERLAY FUNCTION
  # Compute the difference between the year when iopp is detected the earliest, and the year when forest loss occurred. 
  # We do not condition on being within primary forest as of now, in order to remain as general as possible 
  # this is not a loss of coding efficiency as an overlay will be required anyways afterwards to qualify the loss years either as immediate or long.  
  
  # NOTE ON NAs
  # NAs are present, at the margins, due to project operations. 
  # We do not have to reclassify them, as long as we use if_else in the make_time_laps function below. 
  # if_else(NA>0 & NA >0, true = 1, false = 0) does not throw an error. But something using if(NA){} would! 
  
  # Moreover, we would like to stick to INT1 dataType as long as possible to reduce file sizes, but "data type "INT1S" is not available in GDAL."
  # Therefore, we "upgrade" to INT2S. 
  # We arbitrarily chose a value of -51 for cells that do not satisfy the conditions of having loss and having plantations. 
  # This makes it easier (than setting to NA) to reclassify these cells to 0 later on. And of course we cannot set them to 0 because 0 is a value
  # that can be reached in the difference, and has a different meaning (immediate conversion)
  
  # loss rs [[1]]
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
  
  #earliest rs[[2]]
  earliest <- raster(file.path(paste0("temp_data/processed_lu/austin_earliest_detected_",island,"_aligned.tif")))
  
  rs <- raster::stack(loss, earliest)
  
  # condition on the pixel having a loss event once (rs[[1]]>0), within a plantation (rs[[2]]>0)
  # outside this condition, this is not an event we are interested in. 
  make_time_laps <- function(loss1, earliest2){if_else(loss1>0 & earliest2>0, 
                                                       true = earliest2 - loss1-200, # -200 is because values in loss are 1, 2,..., 18 while in earliest they are 200, 205, 210, 215
                                                       false = -51)} 
  
  
  beginCluster() # uses by default detectedCores() - 1
  clusterR(rs,
           fun = overlay, 
           args = list(fun = make_time_laps, unstack = TRUE),
           filename = file.path(paste0("temp_data/processed_lu/earliest_ioppm_to_loss_timelaps_",island,".tif")),
           datatype = "INT2S", 
           overwrite = TRUE)
  endCluster()
  print(paste0("completed earliest_ioppm_to_loss_timelaps_",island,".tif"))
  
  
  
  ### RECLASSIFY
  
  # Finally, split this into  different rasters, with reclassify
  
  # Recall that reclassify works with a reclassiying matrix
  # The two first columns of rclmat are intervals, values within which should be converted in the third column's value.
  # The right = TRUE (the default) means that intervals are open on the left and closed on the right. 
  # Here, with right = FALSE we set this the other way round, so the kind: [x;y[
  # include.lowest means that the lowest interval is closed on the left if right = TRUE, and 
  # that the highest interval is closed on the right if right = FALSE, so in our case here. 
  time_laps <- raster(file.path(paste0("temp_data/processed_lu/earliest_ioppm_to_loss_timelaps_",island,".tif")))
  
  beginCluster() # this uses by default detectCores() - 1
  
  # REPLACEMENT OF TREES WITHIN PLANTATIONS: loss that occurs within plantations. 
  m1 <- c(-51,-50,0,
          -18,0,1,
          0,5,0,
          5,14,0)
  rclmat1 <- matrix(m1, ncol = 3, byrow = TRUE)
  clusterR(time_laps,
           fun = reclassify,
           args = list(rcl = rclmat1, right = FALSE, include.lowest = TRUE),
           #export = "rclmat", # works with or without
           filename = file.path(paste0("temp_data/processed_lu/loss_in_iopp_",island,".tif")),
           datatype = "INT1U",
           overwrite = TRUE)
  rm(m1, rclmat1)
  
  # SHORT TIMELAPS BETWEEN LOSS AND EALIEST IOPP: loss that occurs the same year or up to 4 years prior to the earliest iopp detected
  m2 <- c(-51,-50,0, 
          -18,0,0,
          0,5,1,
          5,14,0)
  rclmat2 <- matrix(m2, ncol = 3, byrow = TRUE)
  clusterR(time_laps,
           fun = reclassify,
           args = list(rcl = rclmat2, right = FALSE, include.lowest = TRUE), 
           #export = "rclmat", # works with or without 
           filename = file.path(paste0("temp_data/processed_lu/loss_to_iopp_rapid_",island,".tif")), 
           datatype = "INT1U",
           overwrite = TRUE)
  rm(m2, rclmat2)
  
  # LONG TIMELAPS BETWEEN LOSS AND EARLIEST IOPP: loss that occurs between 5 to 14 years prior to the earliest iopp detected
  m3 <- c(-51,-50,0, 
          -18,0,0,
          0,5,0,
          5,14,1)
  rclmat3 <- matrix(m3, ncol = 3, byrow = TRUE)
  clusterR(time_laps,
           fun = reclassify,
           args = list(rcl = rclmat3, right = FALSE, include.lowest = TRUE), 
           #export = "rclmat", # works with or without 
           filename = file.path(paste0("temp_data/processed_lu/loss_to_iopp_slow_",island,".tif")), 
           datatype = "INT1U",
           overwrite = TRUE)
  rm(m3, rclmat3)
  
  endCluster()  
  
  print(paste0("completed loss_to_iopp_slow_",island,"_aligned.tif"))
  
  removeTmpFiles(h=0)  
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Overlay forest loss YEARS and different dynamics of oil palm conversion ####
  
  # GFC loss layer (rs[[1]]) of years of a loss event
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
  
  # selecting rasters for either rapid or slow LUCFP (rs[[2]])
  replacement <- raster(file.path(paste0("temp_data/processed_lu/loss_in_iopp_",island,".tif")))
  rapid <- raster(file.path(paste0("temp_data/processed_lu/loss_to_iopp_rapid_",island,".tif")))
  slow <- raster(file.path(paste0("temp_data/processed_lu/loss_to_iopp_slow_",island,".tif")))
  
  # primary forest (rs[[3]])
  pf <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif")))
  
  # stack is necessary for clusterR
  rs_replacement <- raster::stack(loss, replacement)
  rs_rapid <- raster::stack(loss, rapid, pf)
  rs_slow <- raster::stack(loss, slow, pf)
  
  overlay_replacement <- function(rs){rs[[1]]*rs[[2]]}
  overlay_lucpfp <- function(rs){rs[[1]]*rs[[2]]*(rs[[3]] != 0)} # (recall that pf is either 2 or 1 or 0 for degraded, intact or non primary forest resp. Here we want total primary forest)
  
  # note that using calc below is equivalent to using overlay but more appropriate to the input being a stack, 
  # which is necessary to pass several raster layers to the first argument of clusterR
  
  # run the computation in parallel with clusterR, as cells are processed one by one independently.
  beginCluster() # uses by default detectedCores() - 1
  
  # For replacement
  # clusterR(rs_replacement,
  #          fun = calc, #
  #          args = list(overlay_replacement),
  #          filename = file.path(paste0("temp_data/processed_lu/lucpfip_replace_",island,"_total.tif")),
  #          # the name is not very meaningful, but it's for coding purpose along this script, to have this outcome in the loops easily. 
  #          datatype = "INT1U",
  #          overwrite = TRUE )
  
  # For rapid conversion
  clusterR(rs_rapid,
           fun = calc, 
           args = list(overlay_lucpfp),
           filename = file.path(paste0("temp_data/processed_lu/lucpfip_rapid_",island,"_total.tif")),
           datatype = "INT1U",
           overwrite = TRUE )
  
  # For slow conversion
  clusterR(rs_slow,
           fun = calc, 
           args = list(overlay_lucpfp),
           filename = file.path(paste0("temp_data/processed_lu/lucpfip_slow_",island,"_total.tif")),
           datatype = "INT1U",
           overwrite = TRUE )
  
  endCluster()
  
  rm(loss, replacement, rapid, slow, pf, 
     rs_replacement, rs_rapid, rs_slow, overlay_replacement, overlay_lucpfp)
  removeTmpFiles(h=0) 
  
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Split LUCFP to annual layers ####
  
  ### Function description
  # parallel_split has for input a single lucpfip layer where each pixel has a value corresponding to the year when a lucpfip event occured;
  # it outputs annual layers in each of which pixels are either 1 if a lucpfip event occured that year, and 0 else.
  # the tasks are year specific and independent across years, therefore they are executed parallelly over years.
  parallel_split <- function(dyna, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the forest type level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the forest type level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_split <- function(time){
      # define process (island, dyna) we are in 
      process <- file.path(paste0("temp_data/processed_lu/lucpfip_",dyna,"_",island,"_total.tif"))
      
      # #set temp directory
      dir.create(paste0(process,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(process,"_Tmp")))
      
      # read in the input
      lucpfip_prj <- raster(process)
      
      # define output file name 
      output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",dyna,"_",island,"_total_", years[time],".tif"))
      
      # split it into annual binary layers
      calc(lucpfip_prj,
           fun = function(x){if_else(x == time, true = 1, false = 0)},
           filename = output_filename,
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
  
  ### Execute it for each lucpfip dynamic
  
  dynamicS <- c("rapid", "slow") # "replace",
  for(dyna in dynamicS){
    
    parallel_split(dyna = dyna, detectCores() - 1) # ~500 seconds / annual layer
    
    removeTmpFiles(h=0)
  }  
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")  
  
  return(print(paste0("complete prepare_lucpfip_dynamics ", island)))
  
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 2. AGGREGATE THE PIXELS TO A GIVEN PARCEL SIZE. #####

aggregate_lucpfip_dynamics <- function(island, parcel_size){
  
  ### Function description
  # The function has for inputs annual layers of lucfp events at the pixel level.
  # It aggregates these pixels to a parcel size defined by parcel_size (in meters).
  # The aggregation operation is the sum of the pixel lucfp events.
  # Each annual aggregation is tasked in parallel.
  parallel_aggregate <- function(dyna, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the dyna level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the dyna level.
    years <- seq(from = 2001, to = 2015, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_aggregate <- function(time){
      # Define which process (island, dyna, and year) we are in:
      processname <- file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",dyna,"_",island,"_total_", years[time],".tif"))
      
      #set temp directory
      dir.create(paste0(processname,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(processname,"_Tmp")))
      
      # read in the input.
      lucpfip_annual <- raster(processname)
      
      # define output file name
      output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total_",years[time],".tif"))
      
      # aggregate it from the ~30m cells to parcel_size cells with mean function.
      raster::aggregate(lucpfip_annual, fact = c(parcel_size/res(lucpfip_annual)[1], parcel_size/res(lucpfip_annual)[2]),
                        expand = FALSE,
                        fun = sum,
                        na.rm = FALSE, # NA cells are in margins, see the NOTES part. If FALSE, aggregations at margins that use NA 
                        # are discarded because the sum would be spurious as it would count all NA as 0s while it is not necessary the case.
                        # as there is no land on margins anyways, we can set na.rm = TRUE. It is safer, as it will count 
                        filename = output_filename,
                        datatype = "INT4U", # because the sum may go up to ~ 10 000 with parcel_size = 3000,
                        # but to more than 65k with parcel_size = 10000 so INT4U will be necessary;
                        overwrite = TRUE)
      #removes entire temp directory without affecting other running processes (but there should be no temp file now)
      unlink(file.path(paste0(processname,"_Tmp")), recursive = TRUE)
      #unlink(file.path(tmpDir()), recursive = TRUE)
      # return the path to this parcels file
      #return(output_filename)
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
  
  
  ### Execute the function to compute the RasterBrick object of 18 annual layers for each primary forest type
  
  dynamicS <- c("rapid", "slow") # "replace", 
  for(dyna in dynamicS){
    # run the computation, that writes the layers 
    parallel_aggregate(dyna = dyna, ncores = detectCores() - 1)
    
    # brick the layers together and write the brick
    rasterlist <- list.files(path = "temp_data/processed_lu/annual_maps", 
                             pattern = paste0("parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total_"), 
                             full.names = TRUE) %>% as.list()
    
    parcels_brick <- brick(rasterlist)
    
    writeRaster(parcels_brick,
                filename = file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total.tif")),
                datatype = "INT4U",
                overwrite = TRUE)
    
    rm(rasterlist, parcels_brick)
    removeTmpFiles(h=0)
  }
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")    
  
  print(paste0("complete aggregate_lucpfip_dynamics ",island," ",parcel_size/1000, "km"))
}
# 
# timelaps <- raster(file.path(paste0("temp_data/processed_lu/earliest_ioppm_to_loss_timelaps_",island,".tif")))
# timelaps_s <- sampleRandom(timelaps, 10000)
# unique(timelaps_s)
# dataType(timelaps)
# NAvalue(timelaps)
# plot(timelaps)
# 
# # binaries
# replacement <- raster(file.path(paste0("temp_data/processed_lu/loss_in_iopp_",island,".tif")))
# rapid <- raster(file.path(paste0("temp_data/processed_lu/loss_to_iopp_rapid_",island,".tif")))
# slow <- raster(file.path(paste0("temp_data/processed_lu/loss_to_iopp_slow_",island,".tif")))
# plot(slow)
# NAvalue(slow)
# slow_s <- sampleRandom(slow, size = 1e5)
# unique(slow_s)
# dataType(slow)

# 
# lucpfip <- raster(file.path(paste0("temp_data/processed_lu/annual_maps/parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total_2010.tif")))
# 
# lucpfip_v <- getValues(lucpfip)
# lucpfip_v[!is.na(lucpfip_v)] 
# 
# pixels <- raster(file.path(paste0("temp_data/processed_lu/annual_maps/lucpfip_",dyna,"_",island,"_total_2010.tif")))
# plot(pixels)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 3. CONVERT TO A DATAFRAME THE PARCELS WITHIN A GIVEN CATCHMENT RADIUS  #####

# Note on the use of mask and then st_is_within_distance:
# Parcels that are near the coast and hence consist of less ground than others are not discarded by the use 
# of island polygons when making the mask. This is just a particular fixed feature.  
# Mills are selected, as belonging to the island, and then the mask is made by extending large (60km) CR around them. 
# Then, parcels that remain (i.e. are not masked) are not discarded by st_is_within_distance,
# even if parts of their areas are in the sea, as long as their centroids are less than a given distance. 

to_panel_within_CR_dynamics <- function(island, parcel_size, catchment_radius){
  
  ### Function description
  # This repeats roughly twice the same actions, once for parcels within catchment radiuses of IBS mills only,
  # and once for all UML mills.
  
  # raster_to_df converts the raster bricks of annual layers of parcels to a panel dataframe.
  # This is executed for each dyna, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment radius.
  # 3. reshaping the values in these parcels to a long format panel dataframe
  raster_to_df <- function(dyna){
    
    years <- seq(from = 2001, to = 2015, by = 1)
    
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

    ## Mask
    parcels_brick_name <- paste0("parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total")
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
    # the layer index (1-18) indeed corresponds to the year (2001-2018) because, in aggregate_lucpfip, at the end,
    # rasterlist is ordered along 2001-2018 and this order is preserved when the layers are bricked.
    varying_vars <- paste0(parcels_brick_name, "_IBS_masked.", seq(from = 1, to = 15))

    # reshape to long
    m.df <- stats::reshape(ibs_msk_df,
                           varying = varying_vars,
                           v.names = paste0("lucpfip_",dyna,"_pixelcount"), # not precising that it is total primary forest in var name, as we do not compute intact and degraded but only total.
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
            file = file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                    dyna,"_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_IBS_CR_",
                                    "total.rds")))

  }
  
  ### Execute it
  dynamicS <- c("rapid", "slow") # "replace", 
  for(dyna in dynamicS){
    
    raster_to_df(dyna = dyna)
    
    removeTmpFiles(h=0)
  }
  
  
  print(paste0("complete to_panel_within_CR ",island," ",parcel_size/1000,"km ",catchment_radius/1000,"CR"))
  
}



# this one does not depend on catchment radius because we apply a unique 82km maximal CR.
to_panel_within_UML_CR_dynamics <- function(island, parcel_size){
  
  
  ### Function description
  
  # raster_to_df converts the raster bricks of annual layers of parcels to a panel dataframe.
  # This is executed for each dyna, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment radius.
  # 3. reshaping the values in these parcels to a long format panel dataframe
  raster_to_df <- function(dyna){
    
    years <- seq(from = 2001, to = 2015, by = 1)
    
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
    mills_ca <- st_buffer(mills_prj, dist = 82000) # 82km as in Heilmayr et al. 2020
    # work with squares rather than with circles
    for(i in 1:length(mills_ca)){
      mills_ca[i] <- st_as_sfc(st_bbox(mills_ca[i]))
    }
    total_ca <- st_union(st_geometry(mills_ca))
    # coerce to a SpatialPolygon
    total_ca_sp <- as(total_ca, "Spatial")
    # keep mills_prj we need it below
    rm(total_ca, mills_ca, mills)
    
    ## Mask
    parcels_brick_name <- paste0("parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total")
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
    # within <- st_is_within_distance(uml_msk_df, mills_prj, dist = catchment_radius)
    # uml_msk_df <- uml_msk_df %>% dplyr::filter(lengths(within) >0)
    uml_msk_df <- uml_msk_df %>% st_drop_geometry()
    
    rm(parcels_brick)
    
    
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
                           v.names = paste0("lucpfip_",dyna,"_pixelcount"),
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
            file = file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                    dyna,"_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    "82km_UML_CR_",
                                    "total.rds")))
  }
  
  ### Execute it
  dynamicS <- c("rapid", "slow") # "replace",
  for(dyna in dynamicS){
    
    raster_to_df(dyna = dyna)
    
    removeTmpFiles(h=0)
  }
  
  
  print(paste0("complete to_panel_within_UML_CR ",island," ",parcel_size/1000,"km CR"))

}



##### EXECUTE FUNCTIONS AND MERGE THE OUTPUTS #####

#### Execute the functions ####
# Only if their outputs have not been already computed

### Prepare a 30m pixel map of lucfp for each Island
IslandS <- c("Sumatra", "Kalimantan")# 
for(Island in IslandS){
  prepare_pixel_lucpfip_dynamics(Island)
}

### Aggregate this Island map to a chosen parcel size (3km, 6km and 9km for instance)
PS <- 3000  # this was also run with PS = 1000
IslandS <- c("Sumatra", "Kalimantan")#
for(Island in IslandS){
  aggregate_lucpfip_dynamics(island = Island,
                             parcel_size = PS)
}

### For that Island and for each aggregation factor, extract panels of parcels within different catchment area sizes 
# (radius of 10km, 30km and 50km)
PS <- 3000  # this was also run with PS = 1000
IslandS <- c("Sumatra","Kalimantan")#
for(Island in IslandS){
  CR <- 30000 # i.e. 30km radius
  while(CR < 60000){
    to_panel_within_CR_dynamics(island = Island,
                                parcel_size = PS,
                                catchment_radius = CR)
    
    CR <- CR + 20000
  }
}


### Transform to panel data but within the maximal 82km CR from UML 
PS <- 3000  # this was also run with PS = 1000
IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
for(Island in IslandS){
  to_panel_within_UML_CR_dynamics(island = Island,
                                    parcel_size = PS)
} 


#### Gather the lucfip variables for each parcel_size and catchment_radius combinations. ####
PS <- 3000  # this was also run with PS = 1000

### IBS
CR <- 30000 # i.e. 30km radius
while(CR < 60000){
  
  # For each Island, join columns of lucfip variable for different forest definitions 
  df_list <- list()
  IslandS <- c("Sumatra", "Kalimantan")
  for(Island in IslandS){
    
    # df_replace   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_replace_",Island,"_",PS/1000,"km_",CR/1000,"km_IBS_CR_total.rds")))
    df_rapid <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_rapid_",Island,"_",PS/1000,"km_",CR/1000,"km_IBS_CR_total.rds")))
    df_slow    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_slow_",Island,"_",PS/1000,"km_",CR/1000,"km_IBS_CR_total.rds")))
    
    df <- dplyr::select(df_rapid, -lon, -lat, -idncrs_lon, -idncrs_lat)
    df <- inner_join(df_replace, df, by = c("lonlat", "year"))
    
    df_slow <- dplyr::select(df_slow, -lon, -lat, -idncrs_lon, -idncrs_lat)
    df_list[[match(Island, IslandS)]] <- inner_join(df, df_slow, by = c("lonlat", "year"))

    # here, it's normal that df_replace do not have the same size as the two others. 
    if(nrow(df_list[[match(Island, IslandS)]]) != nrow(df_slow)){stop("data frames do not all have the same set of grid cells")}
    rm(df, df_replace, df_rapid, df_slow)
  }
  
  # stack the three Islands together
  indo_df <- bind_rows(df_list)
  
  
  # ### Add columns of converted pixel counts to hectares.
  # pixel_area <- (27.8*27.6)/(1e4)
  # # replace
  # indo_df <- mutate(indo_df, lucpfip_replace_ha = lucpfip_replace_pixelcount*pixel_area) 
  # # rapid
  # indo_df <- mutate(indo_df, lucpfip_rapid_ha = lucpfip_rapid_pixelcount*pixel_area) 
  # # slow
  # indo_df <- mutate(indo_df, lucpfip_slow_ha = lucpfip_slow_pixelcount*pixel_area) 
  # 
  # indo_df <- dplyr::select(indo_df, parcel_id, year, 
  #                          lucpfip_replace_ha,
  #                          lucpfip_rapid_ha, 
  #                          lucpfip_slow_ha,
  #                          lucpfip_replace_pixelcount,
  #                          lucpfip_rapid_pixelcount, 
  #                          lucpfip_slow_pixelcount,
  #                          everything())
  
  
  saveRDS(indo_df, file.path(paste0("temp_data/processed_parcels/lucpfip_panel_dynamics_",PS/1000,"km_",CR/1000,"km_IBS_CR.rds")))
  
  rm(indo_df, df_list)
  CR <- CR + 20000
}

### UML
# For each Island, join columns of lucfip variable for different forest definitions 
df_list <- list()
IslandS <- c("Sumatra", "Kalimantan")
for(Island in IslandS){
  
  # df_replace   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_replace_",Island,"_",PS/1000,"km_82km_UML_CR_total.rds")))
  df_rapid <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_rapid_",Island,"_",PS/1000,"km_82km_UML_CR_total.rds")))
  df_slow    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_slow_",Island,"_",PS/1000,"km_82km_UML_CR_total.rds")))
  
  df <- dplyr::select(df_rapid, -lon, -lat, -idncrs_lon, -idncrs_lat)
  df <- inner_join(df_replace, df, by = c("lonlat", "year"))
  
  df_slow <- dplyr::select(df_slow, -lon, -lat, -idncrs_lon, -idncrs_lat)
  df_list[[match(Island, IslandS)]] <- inner_join(df, df_slow, by = c("lonlat", "year"))
  
  # here, it's normal that df_replace do not have the same size as the two others. 
  if(nrow(df_list[[match(Island, IslandS)]]) != nrow(df_slow)){stop("data frames do not all have the same set of grid cells")}
  rm(df, df_replace, df_rapid, df_slow)
}

# stack the three Islands together
indo_df <- bind_rows(df_list)


# ### Add columns of converted pixel counts to hectares.
# pixel_area <- (27.8*27.6)/(1e4)
# # replace
# indo_df <- mutate(indo_df, lucpfip_replace_ha = lucpfip_replace_pixelcount*pixel_area) 
# # rapid
# indo_df <- mutate(indo_df, lucpfip_rapid_ha = lucpfip_rapid_pixelcount*pixel_area) 
# # slow
# indo_df <- mutate(indo_df, lucpfip_slow_ha = lucpfip_slow_pixelcount*pixel_area) 
# 
# indo_df <- dplyr::select(indo_df, parcel_id, year, 
#                          lucpfip_replace_ha,
#                          lucpfip_rapid_ha, 
#                          lucpfip_slow_ha,
#                          lucpfip_replace_pixelcount,
#                          lucpfip_rapid_pixelcount, 
#                          lucpfip_slow_pixelcount,
#                          everything())


saveRDS(indo_df, file.path(paste0("temp_data/processed_parcels/lucpfip_panel_dynamics_",PS/1000,"km_82km_UML_CR.rds")))

rm(indo_df, df_list)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 