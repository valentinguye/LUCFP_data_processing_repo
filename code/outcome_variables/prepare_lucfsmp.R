### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# 
#       BUILDING PANEL DATAFRAMES OF LAND USE CHANGE FROM PRIMARY FOREST TO INDUSTRIAL PLANTATION (lucfsmp)
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
#   Main outputs:  panel dataframes of lucfsmp pixel-event count in parcels of a given size, from 2001 to 2018,  
#                 for the whole Indonesia (Sumatra, Kalimantan, Papua "row-binded"), 
#                 for 1 forest definition: 30% canopy closure outside industrial plantations in 2000.
#             
#             There is one such dataframe for each combination of parcel size (only 3x3km for now) and catchment radius (10, 30, 50km)
#             ---> lucfsmp_panel_3km_30CR.rds 
#             ---> lucfsmp_panel_3km_50CR.rds
# 
#   
#   Actions:  This script consists of mainly three functions.  
#             0. load needed packages; set working directory; set raster options; define the crs used throughout the script. 
#                             the chunksize and maxmemory raster options should be set accordingly with the machine used,  
#                             considering that this script executes parallel functions using parallel::detectCores() - 1
#                             It is recommended to not change default memory raster options. 
#
#             1. prepare_pixel_lucfsmp(island)
# 
#             2. aggregate_lucfsmp(island, parcel_size)
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
                   "raster", "rgdal", "sp", "sf", "osrm",
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

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
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


### PRELIMINARY

# read in island shapefile
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)


#### Prepare polygons of Indonesian small and medium mosaics of oil palm plantations #### 
smop <- st_read(file.path("input_data/tree_plantations"))
smop <- smop[smop$country == "Indonesia",]
smop <- smop[smop$spec_simp == "Oil palm" | smop$spec_simp == "Oil palm mix",]
# check that we get the same result as in petersen et al. 2016 table 4: 2440 kha of Indonesian oil palm mix area 
sum(smop$area_ha[smop$spec_simp == "Oil palm mix"])/1000
# the first species is always oil palm
unique(smop$spec_1)
# but it's mixed with other species
unique(smop$spec_2) 

# restrict to small and medium sized plantations
# unique(smop$type[smop$type_text == "Mosaic of small-sized plantations"])
# unique(smop$type[smop$type_text == "Mosaic of medium-sized plantations"]) 

sop <- smop[smop$type == 3,] %>% st_geometry()
mop <- smop[smop$type == 2,] %>% st_geometry()

# dissolve by type
sop <- st_union(x = sop, by_feature = FALSE)
mop <- st_union(x = mop, by_feature = FALSE)

sop <- st_transform(sop, crs = indonesian_crs)
mop <- st_transform(mop, crs = indonesian_crs)

# make it Spatial class
sop_sp <- as(sop, "Spatial")
mop_sp <- as(mop, "Spatial")

rm(smop, sop, mop)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 1. PREPARE 30m PIXEL-LEVEL MAPS OF lucfsmp ##### 

prepare_pixel_lucfsmp <- function(island){
  
  #### Overlay forest loss, primary forest, and small and medium sized oil palm plantations ####
  
  ### First step: keep forest loss pixels only OUTSIDE INDUSTRIAL PLANTATIONS 
  
  # We do this only for 30th forest type: 
  # Thus, we define only one overlay function.
  # overlay function
  overlay_maps <- function(rs){rs[[1]]*(1 - rs[[2]])*(1 - rs[[3]])*(rs[[4]]==0)}
  # multiplies a cell of forest loss (rs[[1]]) by 0 (i.e. "removes" it) if it is a plantation in 2000 (rs[[2]])
  # or if it is an industrial plantation in 2015 (rs[[3]])
  # or if it is within primary forest (i.e. if rs[[4]] is either 1 or 2)
  
  ## Read necessary layers and stack them 
  
  # GFC loss layer (rs[[1]])
  loss <- raster(file.path(paste0("temp_data/processed_lu/gfc_loss_",island,"_30th_prj.tif")))
  
  # 2000 industrial plantations (rs[[2]])
  ioppm2000 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif")))
  
  # 2015 industrial plantations (rs[[3]])
  ioppm2015 <- raster(file.path(paste0("temp_data/processed_lu/austin_ioppm_2015_",island,"_aligned.tif")))
  
  # primary forest
  pf <- raster(file.path(paste0("temp_data/processed_lu/margono_primary_forest_",island,"_aligned.tif")))
  
  # stack is necessary for clusterR
  rs <- raster::stack(loss, ioppm2000, ioppm2015, pf)
  
  # note that using calc below is equivalent to using overlay but more appropriate to the input being a stack, 
  # which is necessary to pass several raster layers to the first argument of clusterR
  
  # run the computation in parallel with clusterR, as cells are processed one by one independently.
  beginCluster() # uses by default detectedCores() - 1
  # For 30th forest
  clusterR(rs,
           fun = calc, 
           args = list(overlay_maps),
           filename = file.path(paste0("temp_data/processed_lu/loss_out_30th_",island,".tif")),
           datatype = "INT1U",
           overwrite = TRUE )
  
  endCluster()
  
  rm(loss, ioppm2000, overlay_maps)
  removeTmpFiles(h=0)
  
  
  print(paste0("complete ", "temp_data/processed_lu/loss_out_30th_",island,".tif"))
  
  
  ### Second step: mask this raster of forest loss with polygons of small and medium sized plantations
  
  loss_30th <- raster(file.path(paste0("temp_data/processed_lu/loss_out_30th_",island,".tif")))
  
  # all cells that are not covered by the Spatial object are set to updatevalue
  
  # small sized plantations
  mask(loss_30th,
       sop_sp,
       updatevalue = 0, 
       filename = file.path(paste0("temp_data/processed_lu/lucfsp_",island,"_30th.tif")),
       datatype = "INT1U",
       overwrite = TRUE)
  
  # medium sized plantations
  mask(loss_30th,
       mop_sp,
       updatevalue = 0, 
       filename = file.path(paste0("temp_data/processed_lu/lucfmp_",island,"_30th.tif")),
       datatype = "INT1U",
       overwrite = TRUE)
  
  removeTmpFiles(h=0)
  
  
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
  #### Split LUCFP to annual layers ####
  
  ### Function description
  # parallel_split has for input a single lucfsmp layer where each pixel has a value corresponding to the year when a lucfsmp event occured;
  # it outputs annual layers in each of which pixels are either 1 if a lucfsmp event occured that year, and 0 else.
  # the tasks are year specific and independent across years, therefore they are executed parallely over years.
  parallel_split <- function(plant_type, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the plantation type level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the plantation type level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_split <- function(time){
      # define process (island, plant_type) we are in 
      process <- file.path(paste0("temp_data/processed_lu/lucf",plant_type,"p_",island,"_30th.tif"))
      
      # #set temp directory
      dir.create(paste0(process,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(process,"_Tmp")))
      
      # read in the input
      lucfsmp_prj <- raster(process)
      
      # define output file name 
      output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/lucf",plant_type,"p_",island,"_30th_", years[time],".tif"))
      
      # split it into annual binary layers
      calc(lucfsmp_prj,
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
  
  ### Execute it for each primary forest type
  
  plant_typeS <- c("s", "m")
  for(plant_type in plant_typeS){
    
    parallel_split(plant_type = plant_type, detectCores() - 1) # ~500 seconds / annual layer
    
    removeTmpFiles(h=0)
  }
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")  
  
  return(print(paste0("complete prepare_lucfsmp ", island)))
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 2. AGGREGATE THE PIXELS TO A GIVEN PARCEL SIZE. #####

aggregate_lucfsmp <- function(island, parcel_size){
  
  ### Function description
  # The function has for inputs annual layers of lucfp events at the pixel level.
  # It aggregates these pixels to a parcel size defined by parcel_size (in meters).
  # The aggregation operation is the sum of the pixel lucfp events.
  # Each annual aggregation is tasked in parallel.
  parallel_aggregate <- function(plant_type, ncores){
    
    ## sequence over which to execute the task.
    # We attribute the tasks to CPU "workers" at the annual level and not at the plant_type level.
    # Hence, if a worker is done with its annual task before the others it can move on to the next one and workers' labor is maximized wrt.
    # attributing tasks at the plant_type level.
    years <- seq(from = 2001, to = 2018, by = 1)
    
    ## read the input to the task
    # is done within each task because it is each time different here.
    
    ## define the task
    annual_aggregate <- function(time){
      # Define which process (island, plant_type, and year) we are in:
      processname <- file.path(paste0("temp_data/processed_lu/annual_maps/lucf",plant_type,"p_",island,"_30th_", years[time],".tif"))
      
      #set temp directory
      dir.create(paste0(processname,"_Tmp"), showWarnings = FALSE)
      rasterOptions(tmpdir=file.path(paste0(processname,"_Tmp")))
      
      # read in the input.
      lucfsmp_annual <- raster(processname)
      
      # define output file name
      output_filename <- file.path(paste0("temp_data/processed_lu/annual_maps/parcel_lucf",plant_type,"p_",island,"_",parcel_size/1000,"km_30th_",years[time],".tif"))
      
      # aggregate it from the ~30m cells to parcel_size cells with mean function.
      raster::aggregate(lucfsmp_annual, fact = c(parcel_size/res(lucfsmp_annual)[1], parcel_size/res(lucfsmp_annual)[2]),
                        expand = FALSE,
                        fun = sum,
                        na.rm = FALSE, # NA cells are in margins, see the NOTES part. If FALSE, aggregations at margins that use NA 
                        # are discarded because the sum would be spurious as it would count all NA as 0s while it is not necessary the case.
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
  
  
  ### Execute the function to compute the RasterBrick object of 18 annual layers for each plantation size
  
  plant_typeS <- c("s", "m")
  for(plant_type in plant_typeS){
    # run the computation, that writes the layers 
    parallel_aggregate(plant_type = plant_type, ncores = detectCores() - 1)
    
    # brick the layers together and write the brick
    rasterlist <- list.files(path = "temp_data/processed_lu/annual_maps", 
                             pattern = paste0("parcel_lucf",plant_type,"p_",island,"_",parcel_size/1000,"km_30th_"), 
                             full.names = TRUE) %>% as.list()
    
    parcels_brick <- brick(rasterlist)
    
    writeRaster(parcels_brick,
                filename = file.path(paste0("temp_data/processed_lu/parcel_lucf",plant_type,"p_",island,"_",parcel_size/1000,"km_30th.tif")),
                datatype = "INT4U",
                overwrite = TRUE)
    
    rm(rasterlist, parcels_brick)
    removeTmpFiles(h=0)
  }
  
  rasterOptions(tmpdir = "temp_data/raster_tmp")    
  
  print(paste0("complete aggregate_lucfsmp ",island," ",parcel_size/1000, "km"))
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 3. CONVERT TO A DATAFRAME THE PARCELS WITHIN A GIVEN CATCHMENT RADIUS  #####

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
  # This is executed for each plant_type, iteratively and not in parallel (not necessary because quite fast)
  # The tasks are:
  # 1. masking the brick of parcels of a given size (parcel_size) on a given island with the maximal CA of mills on that island;
  # 2. selecting only the parcels that are within a given catchment radius.
  # 3. reshaping the values in these parcels to a long format panel dataframe
  raster_to_df <- function(plant_type){
    
    years <- seq(from = 2001, to = 2018, by = 1)
    
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
    parcels_brick_name <- paste0("parcel_lucf",plant_type,"p_",island,"_",parcel_size/1000,"km_30th")
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
    # the layer index (1-18) indeed corresponds to the year (2001-2018) because, in aggregate_lucfsmp, at the end, 
    # rasterlist is ordered along 2001-2018 and this order is preserved when the layers are bricked. 
    varying_vars <- paste0(parcels_brick_name, "_IBS_masked.", seq(from = 1, to = 18))
    
    # reshape to long
    m.df <- stats::reshape(ibs_msk_df,
                           varying = varying_vars,
                           v.names = paste0("lucf",plant_type,"p_pixelcount_30th"),
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
            file = file.path(paste0("temp_data/processed_parcels/lucf",plant_type,"p_panel_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_IBS_CR_",
                                    "30th.rds")))
    
    
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
    
    ## Mask
    parcels_brick_name <- paste0("parcel_lucf",plant_type,"p_",island,"_",parcel_size/1000,"km_30th")
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
                           v.names = paste0("lucf",plant_type,"p_pixelcount_30th"),
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
            file = file.path(paste0("temp_data/processed_parcels/lucf",plant_type,"p_panel_",
                                    island,"_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_UML_CR_",
                                    "30th.rds")))
  }
  
  ### Execute it
  plant_typeS <- c("s", "m")
  for(plant_type in plant_typeS){
    
    raster_to_df(plant_type = plant_type)
    
    removeTmpFiles(h=0)
  }  
  
  
  print(paste0("complete to_panel_within_CR ",island," ",parcel_size/1000,"km ",catchment_radius/1000,"CR"))
  
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### EXECUTE FUNCTIONS AND MERGE THE OUTPUTS #####

#### Execute the functions ####
# Only if their outputs have not been already computed

### Prepare a 30m pixel map of lucfp for each Island
IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
for(Island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/annual_maps/lucfmp_",Island,"_30th_2018.tif")))){
    
    prepare_pixel_lucfsmp(Island)
  }
}

### Aggregate this Island map to a chosen parcel size (3km, 6km and 9km for instance)
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
for(Island in IslandS){
  if(!file.exists(file.path(paste0("temp_data/processed_lu/parcel_lucfmp_",Island,"_",PS/1000,"km_30th.tif")))){
    
    aggregate_lucfsmp(island = Island,
                       parcel_size = PS)
  }
}

### For that Island and for each aggregation factor, extract panels of parcels within different catchment area sizes 
# (radius of 10km, 30km and 50km)
PS <- 3000
IslandS <- c("Sumatra", "Kalimantan")#, "Papua"
for(Island in IslandS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
  
  to_panel_within_CR(island = Island,
                     parcel_size = PS,
                     catchment_radius = CR)
  
    CR <- CR + 20000
  }
}



#### Gather the lucfsmp variables for each parcel_size and catchment_radius combinations. ####
PS <- 3000  
sampleS <- c("IBS", "UML")
for(sample in sampleS){
  CR <- 10000 # i.e. 10km radius
  while(CR < 60000){
    
    # For each Island, join columns of plantation type variable 
    df_list <- list()
    IslandS <- c("Sumatra", "Kalimantan")# , "Papua"
    for(Island in IslandS){
      
      df_small   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucfsp_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_30th.rds")))
      df_medium <- readRDS(file.path(paste0("temp_data/processed_parcels/lucfmp_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_30th.rds")))
      
      df_medium <- dplyr::select(df_medium, -lon, -lat, -idncrs_lon, -idncrs_lat)
      df_list[[match(Island, IslandS)]] <- inner_join(df_small, df_medium, by = c("lonlat", "year"))
      
      if(nrow(df_list[[match(Island, IslandS)]]) != nrow(df_medium)){stop("data frames do not all have the same set of grid cells")}
      rm(df_small, df_medium)
    }
    
    # stack the three Islands together
    indo_df <- bind_rows(df_list)
    
    
    # ### Add columns of converted pixel counts to hectares.
    # pixel_area <- (27.8*27.6)/(1e4)
    # # small-sized plantations
    # indo_df <- mutate(indo_df, lucfsp_ha_30th = lucfsp_pixelcount_30th*pixel_area) 
    # # medium-sized plantations
    # indo_df <- mutate(indo_df, lucfmp_ha_30th = lucfmp_pixelcount_30th*pixel_area) 
    # 
    # indo_df <- dplyr::select(indo_df, lonlat, year, 
    #                          lucfsp_ha_30th,
    #                          lucfmp_ha_30th, 
    #                          lucfsp_pixelcount_30th,
    #                          lucfmp_pixelcount_30th,
    #                          everything())
    
    
    saveRDS(indo_df, file.path(paste0("temp_data/processed_parcels/lucfsmp_panel_",PS/1000,"km_",CR/1000,"km_",sample,"_CR.rds")))
    
    rm(indo_df, df_list)
    CR <- CR + 20000
  }
}

