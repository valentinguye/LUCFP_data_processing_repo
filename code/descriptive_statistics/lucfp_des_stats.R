### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#               LUCFP DESCRIPTIVE STATISTICS 
# 
# The aim here is to produce extensive descriptive statistics of the data on land use change from forest to industrial plantations (LUCFP)
#
# There are three parts in this script, one for each land use: 
#
#   I.  Forest cover in 2000
#       We produce six baseline (2000) forest cover maps: GFC forest cover at 30, 60 and 90% canopy closure threshold, outside of 2000 industrial oil palm plantations 
#                                                         GFC forest cover at 30% canopy closure threshold within intact, degraded, and total primary forest.  
# 
# 
#       Then, we extract the forest extent within different polygons (islands, mill catchment radiuses)
#       
#     Input:  - Unprojected GFC data for three different canopy closure thresholds, as outputted by prepare_gfc.R
#             ---> temp_data/processed_lu/gfc_data_Indonesia_30th.tif, gfc_data_Indonesia_60th.tif, gfc_data_Indonesia_90th.tif 
# 
#             - Projected industrial oil palm plantations in 2000 (Austin et al. 2017)
#             ---> temp_data/processed_lu/austin_ioppm_2000_",island,"_aligned.tif
# 
#             - Primary forest 
#             ---> margono_primary_forest_island_aligned.tif 
#   
#                      polygon extract direct
#
#
#   II. Accumulated LUC to industrial plantations    
#       Over each of six baseline forest types. 
#     
#       input:  - brick of annual maps of parcels of lucfip events 
#               ---> temp_data/processed_lu/parcel_lucfip_island_parcel_size/1000km_th.tif
#         
#               - lucPfip equivalent
#               ---> temp_data/processed_lu/parcel_lucpfip_island_parcel_size/1000km_pf_type.tif
# 
# 
#               sum over years and extract in polygons
#
# 
#   III.  Accumulated LUC to large, medium and small size oil palm plantations.  
# 
# 
# 
# 
#
# 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
# PB: on n'arrive pas à installed leaflet, kableExtra (qui n'est peut être pas nécessaire) et velox (qui n'est pas nécessaire)
neededPackages = c("dplyr", "data.table",  "stringr", "sjmisc",
                   "foreign", "readxl", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
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

# # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
troublePackages <- c("leaflet", "kableExtra", "leaflet.providers", "png")
# Attempt to load packages from user's default libraries.
lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 



### NEW FOLDERS USED IN THIS SCRIPT 

### PARCEL SIZE in METERS
parcel_size <- 3000

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

##### Prepare polygons of interest #####

# st_write(island_sf_prj, "temp_data/processed_indonesia_spatial/island_sf_prj", 
#          driver = "ESRI Shapefile", 
#          delete_dsn = TRUE, 
#          overwrite = TRUE)

### Prepare polygons of mill catchment radius ###
# We would like a multipolygon Simple feature collection with X features and 1 field
# The one field describes the feature (i.e. polygon).
# There is one feature for each CR for each mill group (IBS or UML) within each island (so 3*2*3 = 18 features). 

  ## Prepare UML mills
  uml <- read_xlsx(file.path("input_data/uml/mills_20200129.xlsx"))
  uml$latitude <- as.numeric(uml$latitude)
  uml$longitude <- as.numeric(uml$longitude)
  uml$lat <- uml$latitude
  uml$lon <- uml$longitude
  uml <- st_as_sf(uml,	coords	=	c("longitude",	"latitude"), crs = 4326)
  uml <- st_geometry(uml)
  uml_prj <- st_transform(uml, crs = indonesian_crs)
  
  ## Prepare IBS mills
  ibs <- read.dta13(file.path("temp_data/processed_mill_geolocalization/IBS_UML_panel.dta"))
  # keep only a cross section of those that are geolocalized mills
  ibs <- ibs[ibs$analysis_sample==1,]
  ibs <- ibs[!duplicated(ibs$firm_id),]
  
  ibs <- st_as_sf(ibs,	coords	=	c("lon",	"lat"), crs=4326)
  ibs <- st_geometry(ibs)
  ibs_prj <- st_transform(ibs, crs = indonesian_crs)
  


# Fill the lists with shapes and shape descriptions (names)
IslandS <- c("Sumatra", "Kalimantan", "Papua")
catchment_radius <- c(10, 30, 50)

for(island in IslandS){
  # Create empty list for shapes and shape descriptions
  shape_des <- list()
  mill_cr <- list()
  
  ## UML lists
  shape_des[[1]] <- list()
  mill_cr[[1]] <- list()
  
  for(CR in catchment_radius){
    ## create union of individual catchment areas.
    
    # buffer
    shape <- st_buffer(uml_prj, dist = CR*1000)
    
    # # work with squares rather than circles. 
    # shape <- sapply(shape, FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
    
    shape <- st_union(shape)
    
    # keep only the part of this total catchment area that is on our island of interest
    shape <- st_intersection(x = shape, y = island_sf_prj[island_sf_prj$shape_des == island,])
    
    # add name and geometry to lists
    shape_des[[1]][[match(CR, catchment_radius)]] <- paste0(island,"_uml_",CR,"km")
    mill_cr[[1]][[match(CR, catchment_radius)]] <-  shape %>% st_geometry()
  }  
  
  ## IBS lists
  shape_des[[2]] <- list()
  mill_cr[[2]] <- list()
  
  for(CR in catchment_radius){
    ## create union of individual catchment areas.
    # buffer
    shape <- st_buffer(ibs_prj, dist = CR*1000)
    
    # work with squares rather than circles. 
    shape <- sapply(shape, FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
    
    shape <- st_union(shape)
    
    # # keep only the part of this total catchment area that is on our island of interest
    # shape <- st_intersection(x = shape, y = island_sf_prj[island_sf_prj$shape_des == island,])

    # add name and geometry to lists
    shape_des[[2]][[match(CR, catchment_radius)]] <- paste0(island,"_ibs_",CR,"km")
    mill_cr[[2]][[match(CR, catchment_radius)]] <-  shape %>% st_geometry()
  }  
  
  # flatten the lists
  shape_des <- unlist(shape_des)
  mill_cru <- unlist(mill_cr, recursive = FALSE) 
  mill_cru <- unlist(mill_cru, recursive = FALSE)
  
  # make the sf object
  mill_cr_sf_prj <- st_set_geometry(data.frame(shape_des), st_as_sfc(mill_cru))
  st_crs(mill_cr_sf_prj) <- indonesian_crs
  
  # add the total island shape
  all_shapes <- rbind(mill_cr_sf_prj, island_sf_prj[island_sf_prj$shape_des==island,])
  
  # Save the shapes 
  dir.create(file.path(paste0("temp_data/processed_indonesia_spatial/mill_cr_prj_",island)))
  st_write(all_shapes, file.path(paste0("temp_data/processed_indonesia_spatial/mill_cr_prj_",island)), 
           driver = "ESRI Shapefile", 
           delete_dsn = TRUE, 
           overwrite = TRUE)
}



##### 1. FORET COVER IN 2000 ##### 
#### Prepare 30, 60, 90% forest cover outside industrial plantations in 2000 ####

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



#### Prepare intact, degraded and total primary forest cover in 2000 ####
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









#### Extract pixel sums in polygons #### 
# We have 7 polygons / island (3 CR sizes for 2 groups of mills (IBS and total), plus the total island polygon). 
# This makes 21 polygons in which to extract from 6 layers of 2000 forest cover.
# For efficiency reasons, we do this in GEE, because all polygons are pretty large and even a single extract 
# operation is long in R (and there 126 of them).  
# - GEE code for extraction on forest of different canopy densities outside indsutrial plantations in 2000
# https://code.earthengine.google.com/?scriptPath=users%2Fvalentinguye%2Frepo1%3Afc2000_outside_plantations_extractions
# - GEE code for extraction on different primary forest types
# https://code.earthengine.google.com/?scriptPath=users%2Fvalentinguye%2Frepo1%3Aprimary_forest_cover_2000





#### Make an ordered data frame that is easily tabulized afterwards. ####

# We proceed by blocks of island x forest definition. 

# FOR EACH ISLAND X FOREST DEFINITION, EXTRACTION IN GEE OUTPUTS DATAFRAMES WITH THREE COLUMNS AND 7 LINES 
# Each forest definition is an element of this list
IslandS <- c("Sumatra", "Kalimantan", "Papua")

pfc2000_list <- list()

for(island in IslandS){
  pc_intact <- st_read(file.path(paste0("input_data/GEE_outputs/pixelcount_",island,"_intact.geojson"))) %>% st_drop_geometry()  
  pc_degraded <- st_read(file.path(paste0("input_data/GEE_outputs/pixelcount_",island,"_degraded.geojson"))) %>% st_drop_geometry()
  pc_total <- st_read(file.path(paste0("input_data/GEE_outputs/pixelcount_",island,"_total.geojson"))) %>% st_drop_geometry()
  
  # change column names from sum to its respective forest definition
  names(pc_intact)[names(pc_intact)=="sum"] <- "intact"
  names(pc_degraded)[names(pc_degraded)=="sum"] <- "degraded"
  names(pc_total)[names(pc_total)=="sum"] <- "total"
  
  # coerce shape_des column to character (from vector) 
  pc_intact$shape_des <- as.character(pc_intact$shape_des)
  pc_degraded$shape_des <- as.character(pc_degraded$shape_des)
  pc_total$shape_des <- as.character(pc_total$shape_des)
  
  # useless pixel id column
  pc_intact <- dplyr::select(pc_intact, -id)
  pc_degraded <- dplyr::select(pc_degraded, -id)
  pc_total <- dplyr::select(pc_total, -id)
  
  # merge forest definitions
  pfc2000_tmp <- merge(pc_intact, pc_degraded, by = "shape_des")
  pfc2000_list[[match(island, IslandS)]] <- merge(pfc2000_tmp, pc_total, by = "shape_des")
}
pfc2000 <- bind_rows(pfc2000_list)

# now repeat for non-ip 30, 60 and 90 forest cover 2000 and cbind 
# something like : (tc : tree cover)
tc2000_list <- list()

for(island in IslandS){
  tc_30th <- st_read(file.path(paste0("input_data/GEE_outputs/pixelcount_",island,"_30th.geojson"))) %>% st_drop_geometry() 
  tc_60th <- st_read(file.path(paste0("input_data/GEE_outputs/pixelcount_",island,"_60th.geojson"))) %>% st_drop_geometry()
  tc_90th <- st_read(file.path(paste0("input_data/GEE_outputs/pixelcount_",island,"_90th.geojson"))) %>% st_drop_geometry()
  
  # change column names from sum to its respective forest definition
  names(tc_30th)[names(tc_30th)=="sum"] <- "30%"
  names(tc_60th)[names(tc_60th)=="sum"] <- "60%"
  names(tc_90th)[names(tc_90th)=="sum"] <- "90%"
  
  # coerce shape_des column to character (from vector) 
  tc_30th$shape_des <- as.character(tc_30th$shape_des)
  tc_60th$shape_des <- as.character(tc_60th$shape_des)
  tc_90th$shape_des <- as.character(tc_90th$shape_des)
  
  # useless pixel id column
  tc_30th <- dplyr::select(tc_30th, -id)
  tc_60th <- dplyr::select(tc_60th, -id)
  tc_90th <- dplyr::select(tc_90th, -id)
  
  # merge forest definitions (put 90% on the left)
  tc2000_tmp <- merge(tc_90th, tc_60th, by = "shape_des")
  tc2000_list[[match(island, IslandS)]] <- merge(tc2000_tmp, tc_30th, by = "shape_des")
}
tc2000 <- bind_rows(tc2000_list)

# and finally 
fc2000 <- merge(pfc2000, tc2000, by = "shape_des")

fc2000_torestore <- fc2000
#fc2000 <- fc2000_torestore

# then prepare row names, matrix etc. 
row.names(fc2000) <- fc2000$shape_des
fc2000 <- fc2000[c("Sumatra_ibs_10km",
                         "Sumatra_ibs_30km",
                         "Sumatra_ibs_50km",
                         "Sumatra_uml_10km",
                         "Sumatra_uml_30km",
                         "Sumatra_uml_50km",
                         "Sumatra",
                         "Kalimantan_ibs_10km",
                         "Kalimantan_ibs_30km",
                         "Kalimantan_ibs_50km",
                         "Kalimantan_uml_10km",
                         "Kalimantan_uml_30km",
                         "Kalimantan_uml_50km",
                         "Kalimantan",
                         "Papua_ibs_10km",
                         "Papua_ibs_30km",
                         "Papua_ibs_50km",
                         "Papua_uml_10km",
                         "Papua_uml_30km",
                         "Papua_uml_50km",
                         "Papua"), ]

# remove shape_des column 
fc2000 <- dplyr::select(fc2000, -shape_des)

# simplifying row names makes them not unique, which is not accepted in a data frame but is in a matrix. 
fc2000 <- as.matrix(fc2000)

# Fill the row of total over islands
fc2000 <- rbind(fc2000, rep(-1, ncol(fc2000)))
row.names(fc2000)[nrow(fc2000)] <- "Total three islands"
for(col in c(1:6)){
  fc2000["Total three islands", col] <- sum(fc2000["Sumatra", col], 
                                            fc2000["Kalimantan", col],
                                            fc2000["Papua", col])
}

# convert the sum value from number of pixels of resolution 27.8 by 27.6 meters (i.e. ~767.29 square meters) to million hectares 
pixel_area <- (27.8*27.6)/(1e10) # 1e10 is the convertion factor between a square meter and a MILLION hectare.  

fc2000 <- fc2000*pixel_area
fc2000 <- fc2000 %>% round(digits = 2)

row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                 pattern = "Sumatra_", 
                                 replacement = "")
row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                 pattern = "Sumatra", 
                                 replacement = "Total island")
row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                 pattern = "Kalimantan_", 
                                 replacement = "")
row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                 pattern = "Kalimantan", 
                                 replacement = "Total island")
row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                 pattern = "Papua_", 
                                 replacement = "")
row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                 pattern = "Papua", 
                                 replacement = "Total island")

for(KM in c(10, 30, 50)){
  row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                   pattern = paste0("ibs_",KM,"km"), 
                                   replacement = paste0("Within ",KM,"km of a sample mill"))
  
  
  row.names(fc2000) <- str_replace(string = row.names(fc2000), 
                                   pattern = paste0("uml_",KM,"km"), 
                                   replacement = paste0("Within ",KM,"km of an UML mill"))
}



#### Print the LateX table code ####
options(knitr.table.format = "latex") 

LU_stat_des <- fc2000
colnames(LU_stat_des) <- NULL

kable(LU_stat_des, booktabs = T, align = "c", 
      caption = "Forest extents in 2000, Mha") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "intact" = 1, "degraded" = 1, "total" = 1, 
                     "90%" = 1, "60%" = 1, "30%" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "30% pixel tree cover threshold,\nwithin primary forest" = 3, 
                     "Different pixel tree cover thresholds,\noutside industrial oil palm plantations" = 3), 
                   bold = T, 
                   align = "c") %>% 
  pack_rows("Sumatra", 1, 7)  %>%
  pack_rows("Kalimantan", 8, 14)  %>%
  pack_rows("Papua", 15, 21)  %>%
  column_spec(column = c(2:7),
              width = "5em") %>% 
  footnote(general = c("Sample mills are the 470 IBS palm oil mills we matched with the UML database."),
           threeparttable = TRUE, 
           escape = TRUE) 





##### 2. ACCUMULATED LUCFP #####

# ibs <- mills
# ibs <- ibs[!is.na(ibs$lat),]
# ibs <- st_as_sf(ibs, coords = c("lon", "lat"), crs = 4326)
# ibs[ibs$island_name == "Papua" & ibs$uml_matched_sample==1, "geometry"] %>% plot()

# Make the template table

accu_lucfp <- matrix(-1, nrow = 26, ncol = 6)
row.names(accu_lucfp) <- c("Sumatra_IBS_10km",
                           "Sumatra_IBS_30km",
                           "Sumatra_IBS_50km",
                           "Sumatra_UML_10km",
                           "Sumatra_UML_30km",
                           "Sumatra_UML_50km",
                           "Sumatra",
                           "",
                           "Kalimantan_IBS_10km",
                           "Kalimantan_IBS_30km",
                           "Kalimantan_IBS_50km",
                           "Kalimantan_UML_10km",
                           "Kalimantan_UML_30km",
                           "Kalimantan_UML_50km",
                           "Kalimantan",
                           "",
                           "Papua_IBS_10km",
                           "Papua_IBS_30km",
                           "Papua_IBS_50km",
                           "Papua_UML_10km",
                           "Papua_UML_30km",
                           "Papua_UML_50km",
                           "Papua", 
                           "", 
                           "Total three islands", 
                           "")
colnames(accu_lucfp) <- c("intact", "degraded", "total", "90%", "60%", "30%")



IslandS <- c("Sumatra", "Kalimantan", "Papua")
PS <- 3000
catchment_radiuseS <- c(10, 30, 50)
sampleS <- c("IBS", "UML")
pf_typeS <- c("intact", "degraded", "total")
thresholdS <- c(30, 60, 90)

for(island in IslandS){
  
  ## Extract total in different catchment radii using panel dataframes 
  # An initial way to do it was to use data frames from prepare_lucfip.R and prepare_lucpfip.R scripts
  # which feature all parcels that are within the given catchment radius *at least one year*. 
  # Therefore, aggregation using these dataframes includes lucfp in places and years where no mill is already established.
  # This is of course interesting, but it's not what we want to call "lucfp within at least one mill's catchment radius". 
  #  ---> therefore, this first method is commented out and we now compute accumulated LUCFP under the condition 
  # that at least one mill is reachable. This implies to use a different dataframe, that features n_reachable vars.  
  
  # this is a bit longer because it needs to load several times the heavier dataframe 
  
  # LUCPFIP
  for(sample in sampleS){
    for(CR in catchment_radiuseS){
      for(pf_type in pf_typeS){
        
        # read panel data frames of parcels prepared in prepare_lucpfip.R 
        # df <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
        #                                island,"_",PS/1000,"km_",CR,"km","_",sample,"_CR_",pf_type,".rds")))
        # accu_lucfp[paste0(island,"_",sample,"_",CR,"km"), pf_type] <- sum(df[df$year <= 2015, 5], na.rm = TRUE)
        # rm(df)
        
        # read in all parcels (panel from 2001-2015) within a given catchment radius at least one year 
        parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                            PS/1000,"km_",
                                            CR,"CR.rds")))
        # restrict to island
        parcels <- parcels[parcels$island == island,]
        
        # sum pixelcounts over parcels and years, under the condition that n_reachable is not null
        if(sample == "IBS"){
         accu_lucfp[paste0(island,"_",sample,"_",CR,"km"), pf_type] <- sum(parcels[parcels$n_reachable_ibs > 0, paste0("lucpfip_pixelcount_",pf_type)])
        }
        if(sample == "UML"){
         accu_lucfp[paste0(island,"_",sample,"_",CR,"km"), pf_type] <- sum(parcels[parcels$n_reachable_uml > 0, paste0("lucpfip_pixelcount_",pf_type)])
        }

      }
    }
  }
  
  # LUCFIP
  for(sample in sampleS){
    for(CR in catchment_radiuseS){
      for(th in thresholdS){
        
        # read panel data frames of parcels prepared in prepare_lucfip.R 
        # df <- readRDS(file.path(paste0("temp_data/processed_parcels/lucfip_panel_",
        #                                island,"_",PS/1000,"km_",CR,"km","_",sample,"_CR_",th,"th.rds")))
        # # these dataframes are all shaped in the same way, that the fifth column is always the pixelcount we are interested in summing here.
        # accu_lucfp[paste0(island,"_",sample,"_",CR,"km"), paste0(th,"%")] <- sum(df[df$year <= 2015, 5], na.rm = TRUE)
        # rm(df)
        
        # read in all parcels (panel from 2001-2015) within a given catchment radius at least one year 
        parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                            PS/1000,"km_",
                                            CR,"CR.rds")))
        # restrict to island
        parcels <- parcels[parcels$island == island,]
        
        # sum pixelcounts over parcels and years, under the condition that n_reachable is not null
        if(sample == "IBS"){
          accu_lucfp[paste0(island,"_",sample,"_",CR,"km"), paste0(th,"%")] <- sum(parcels[parcels$n_reachable_ibs > 0, paste0("lucfip_pixelcount_",th,"th")])
        }
        if(sample == "UML"){
          accu_lucfp[paste0(island,"_",sample,"_",CR,"km"), paste0(th,"%")] <- sum(parcels[parcels$n_reachable_uml > 0, paste0("lucfip_pixelcount_",th,"th")])
        }
        
      }
    }
  }
  
  ## Extract the island total using raster layers of parcels 
  
  # LUCPFIP
  for(pf_type in pf_typeS){
    brick_lucpfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",
                                          island,"_",PS/1000,"km_",pf_type,".tif")))
    
    # select 2001-2015 layers 
    layer_names <- paste0("parcel_lucpfip_",island,"_",PS/1000,"km_",pf_type,".",c(1:15))
                                          
    brick_lucpfip <- raster::subset(brick_lucpfip, layer_names)
    
    # Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
    r_accu_lucfp <- calc(brick_lucpfip, fun = sum, na.rm = TRUE)
    
    accu_lucfp[island, pf_type] <- raster::extract(x = r_accu_lucfp, 
                                                   y = island_sf_prj[island_sf_prj$shape_des==island,] %>% st_geometry() %>% as("Spatial"), 
                                                   fun = sum, 
                                                   na.rm = TRUE) 
  }

  # LUCFIP
  for(th in thresholdS){
    brick_lucfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucfip_",
                                            island,"_",PS/1000,"km_",th,"th.tif")))
    
    # select 2001-2015 layers 
    layer_names <- paste0("parcel_lucfip_",island,"_",PS/1000,"km_",th,"th.",c(1:15))
    
    brick_lucfip <- raster::subset(brick_lucfip, layer_names)
    
    # Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
    r_accu_lucfp <- calc(brick_lucfip, fun = sum, na.rm = TRUE)
    
    accu_lucfp[island, paste0(th,"%")] <- raster::extract(x = r_accu_lucfp, 
                                                   y = island_sf_prj[island_sf_prj$shape_des==island,] %>% st_geometry() %>% as("Spatial"), 
                                                   fun = sum, 
                                                   na.rm = TRUE) 
  }  
}


accu_lucfp_torestore <- accu_lucfp
# accu_lucfp <- accu_lucfp_torestore

### ### ### ### Alternative method ### ### ### ### 


# mais est ce qu'on a pas un problème car 
### ### ### ### ### ### ### ### ### 

# Fill the row of total over islands
for(col in c(1:6)){
  accu_lucfp["Total three islands", col] <- sum(accu_lucfp["Sumatra", col], 
                                                accu_lucfp["Kalimantan", col],
                                                accu_lucfp["Papua", col])
}

# convert the sum value from number of pixels of resolution 27.8 by 27.6 meters (i.e. ~767.29 square meters) to thousand hectares 
pixel_area <- (27.8*27.6)/(1e7) # 1e7 is the convertion factor between a square meter and a THOUSAND hectare.  

accu_lucfp <- accu_lucfp*pixel_area
accu_lucfp <- accu_lucfp %>% round(digits = 2)

# simplifying row names makes them not unique, which is not accepted in a data frame but is in a matrix. 
row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                pattern = "Sumatra_", 
                                replacement = "")
row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                pattern = "Sumatra", 
                                replacement = "Total island")
row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                pattern = "Kalimantan_", 
                                replacement = "")
row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                pattern = "Kalimantan", 
                                replacement = "Total island")
row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                pattern = "Papua_", 
                                replacement = "")
row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                pattern = "Papua", 
                                replacement = "Total island")

for(KM in c(10, 30, 50)){
  row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                  pattern = paste0("IBS_",KM,"km"), 
                                  replacement = paste0("Within ",KM,"km of a sample mill"))
  
  
  row.names(accu_lucfp) <- str_replace(string = row.names(accu_lucfp), 
                                  pattern = paste0("UML_",KM,"km"), 
                                  replacement = paste0("Within ",KM,"km of an UML mill"))
}


# Finally, add figures from Austin et al. 2017 (SI)
accu_lucfp[8,"intact"] <- "(35)"
accu_lucfp[16,"intact"] <- "(14)"
accu_lucfp[24,"intact"] <- "(36)"
accu_lucfp[26,"intact"] <- "(85)"

accu_lucfp[8,"degraded"] <- "(413)"
accu_lucfp[16,"degraded"] <- "(1007)"
accu_lucfp[24,"degraded"] <- "(55)"
accu_lucfp[26,"degraded"] <- "(1473)"

accu_lucfp[8,"total"] <- "(448)"
accu_lucfp[16,"total"] <- "(1021)"
accu_lucfp[24,"total"] <- "(91)"
accu_lucfp[26,"total"] <- "(1558)"

accu_lucfp[c(8,16,24,26), c("90%", "60%", "30%")] <- ""

#### Print the LateX table code ####
options(knitr.table.format = "latex") 

LU_stat_des <- accu_lucfp
colnames(LU_stat_des) <- NULL

kable(LU_stat_des, booktabs = T, align = "c", 
      caption = "Land use change from different forest definitions to industrial oil palm plantations, accumulated over 2001-2015, in kha.") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "intact" = 1, "degraded" = 1, "total" = 1, 
                     "90%" = 1, "60%" = 1, "30%" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "From 30% pixel tree cover threshold,\nwithin primary forest in 2000" = 3, 
                     "From different pixel tree cover thresholds,\noutside industrial plantations in 2000" = 3), 
                   bold = T, 
                   align = "c") %>% 
  # add_header_above(c(" " = 1,
  #                    "LUC to industrial plantations \n accumulated 2001-2015, kha" = 6)) %>% 
  pack_rows("Sumatra", 1, 8)  %>%
  pack_rows("Kalimantan", 9, 16)  %>%
  pack_rows("Papua", 17, 24)  %>%
  column_spec(column = c(2:7),
              width = "5em") %>% 
  footnote(general = c(""),
           threeparttable = TRUE, 
           escape = TRUE) 





##### Make figures #####

# IBS
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
ibs <- ibs[ibs$analysis_sample == 1,]
ibs <- ibs[!duplicated(ibs$firm_id),]
length(unique(ibs$firm_id))

ibs <- st_as_sf(ibs,	coords	=	c("lon",	"lat"), crs=4326)
ibs <- st_geometry(ibs)
ibs_prj <- st_transform(ibs, crs = indonesian_crs)
ibs30 <- st_buffer(ibs_prj, dist = 30000)
ibs30 <- st_union(ibs30)

#ibs30 <- st_intersection(x = ibs30, y = island_sf_prj)

# keep only the part of this total catchment area that is on our island of interest
#ibs30 <- st_intersection(x = ibs30, y = island_sf_prj[island_sf_prj$shape_des == island,])
# un-project to leaflet lon lat 
ibs30_lonlat <- st_transform(ibs30, crs = 4326)
# (For raster we let leaflet do it) 

# UML
uml <- read_xlsx(file.path("input_data/uml/mills_20200129.xlsx"))
uml$latitude <- as.numeric(uml$latitude)
uml$longitude <- as.numeric(uml$longitude)
uml$lat <- uml$latitude
uml$lon <- uml$longitude
uml <- st_as_sf(uml,	coords	=	c("longitude",	"latitude"), crs = 4326)
uml <- st_geometry(uml)
uml_prj <- st_transform(uml, crs = indonesian_crs)

uml30 <- st_buffer(uml_prj, dist = 30000)
uml30 <- st_union(uml30)
#uml30 <- st_intersection(x = uml30, y = island_sf_prj)
uml30_lonlat <- st_transform(uml30, crs = 4326)


#### MAP LUCFIP ####
prepare_accu_lucfp <- function(island){
  parcel_size <- 3000
  th <- 30
  brick_lucfp <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucfip_",island,"_",parcel_size/1000,"km_",th,"th.tif")))
  
  # remove layers of years aftr 2015
  brick_lucfp <- raster::subset(brick_lucfp, c(1:15))
  # Add up annual aggregated LUCFP
  accu_lucfp <- calc(brick_lucfp, fun = sum, na.rm = TRUE)
  # convert the sum value from number of pixels of resolution 27.7m (i.e; 767.29 square meters) to million hectare 
  pixel_area <- (27.8*27.6)/(1e4) # 1e4 is the convertion factor between a square meter and ONE hectare.  
  accu_lucfp <- calc(accu_lucfp, fun = function(x){x*pixel_area})
  # turn 0 (which are most parcels in the raster) to NA for transparence
  accu_lucfp <- reclassify(accu_lucfp, rcl = cbind(0,NA))
  
  return(accu_lucfp)
}
# plot(accu_lucfp)
# ibs30 %>% st_geometry() %>% plot(add = TRUE)
# island_sf_prj[island_sf_prj$shape_des == "Sumatra", "geometry"] %>% plot(add = T)

accu_lucfp_suma <- prepare_accu_lucfp("Sumatra")
accu_lucfp_kali <- prepare_accu_lucfp("Kalimantan")
accu_lucfp_papu <- prepare_accu_lucfp("Papua")

# settings for lucfp legend 
bins <- seq(from = 0, to = 900, by = 300)
cb <- colorBin("plasma", 
                domain = bins, 
                bins = bins, 
                na.color = "transparent")
# "viridis", "magma", "inferno", or "plasma".
# cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
#                     na.color = "transparent")

# this does not do anything: 
# bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
# bin_labels <- paste0("<div style='display: inline-block;height: ", 
#                      size, "px;margin-top: 4px;line-height: ", 
#                      size, "px;'>", bin_labels, "</div>")

# settings for catchment radius legend
color <- "transparent"
label <- "30km catchment <br/> radius of IBS mills"
shape <- "circle"
border <- "red"
size <- 10

shape <- gsub("circle", "50%", shape)
legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)

# MAP
ibs30_lonlat%>% 
  leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  addPolygons(opacity = 0.5, color = "red", weight = 2, fill = FALSE) %>%
  addRasterImage(accu_lucfp_suma, project = TRUE, colors = cb) %>% 
  addRasterImage(accu_lucfp_kali, project = TRUE, colors = cb) %>% 
  addRasterImage(accu_lucfp_papu, project = TRUE, colors = cb) %>% 
  addLegend(pal = cb,  values = cb, opacity = 0.7,
            labFormat = labelFormat(suffix = "ha"),
            labels = bin_labels, # this lign does not do anything
            title = "Total 2001-2015 LUCFP <br/> within 900ha parcels",
            position = "topright") %>% 
  addLegend(colors = legend_colors, labels = label) 
            

#### MAP LUCPFIP ####

prepare_accu_lucpfip <- function(island){
  parcel_size <- 3000
  pf_type <- "total"
  brick_lucpfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,".tif")))
  
  # remove layers of years aftr 2015
  brick_lucpfip <- raster::subset(brick_lucpfip, c(1:15))
  # Add up annual aggregated lucpfip
  accu_lucpfip <- calc(brick_lucpfip, fun = sum, na.rm = TRUE)
  # convert the sum value from number of pixels of resolution 27.7m (i.e; 767.29 square meters) to million hectare 
  pixel_area <- (27.8*27.6)/(1e4) # 1e4 is the convertion factor between a square meter and ONE hectare.  
  accu_lucpfip <- calc(accu_lucpfip, fun = function(x){x*pixel_area})
  # turn 0 (which are most parcels in the raster) to NA for transparence
  accu_lucpfip <- reclassify(accu_lucpfip, rcl = cbind(0,NA))
  
  return(accu_lucpfip)
}
# plot(accu_lucpfip)
# ibs30 %>% st_geometry() %>% plot(add = TRUE)
# island_sf_prj[island_sf_prj$shape_des == "Sumatra", "geometry"] %>% plot(add = T)

accu_lucpfip_suma <- prepare_accu_lucpfip("Sumatra")
accu_lucpfip_kali <- prepare_accu_lucpfip("Kalimantan")
accu_lucpfip_papu <- prepare_accu_lucpfip("Papua")

max_lucpfip <- max(accu_lucpfip_suma@data@max, 
                   accu_lucpfip_kali@data@max, 
                   accu_lucpfip_papu@data@max)

# min is not 0 because we've turned 0s into NAs for transparence. 
min_lucpfip <- min(accu_lucpfip_suma@data@min, 
                   accu_lucpfip_kali@data@min, 
                   accu_lucpfip_papu@data@min)

dom <- c(getValues(accu_lucpfip_suma), getValues(accu_lucpfip_kali), getValues(accu_lucpfip_papu))

# settings for lucpfip legend 
# bins <- seq(from = 0, to = 900, by = 300)
cb <- colorNumeric(palette = "plasma", 
                   domain = dom, 
                   #bins = bins, 
                   na.color = "transparent")
# "viridis", "magma", "inferno", or "plasma".
# cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucpfip_suma),
#                     na.color = "transparent")

# this does not do anything: 
# bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
# bin_labels <- paste0("<div style='display: inline-block;height: ", 
#                      size, "px;margin-top: 4px;line-height: ", 
#                      size, "px;'>", bin_labels, "</div>")

# settings for catchment radius legend - IBS
color <- "transparent"
label <- "30km catchment radius - IBS geo-localized mills"
shape <- "circle"
border <- "black"
size <- 10

shape <- gsub("circle", "50%", shape)
legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)

# settings for catchment radius legend - UML
color <- "transparent"
label_uml <- "30km catchment radius - UML mills"
shape <- "circle"
border <- "red"
size <- 10

shape <- gsub("circle", "50%", shape)
legend_colors_uml <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)

# MAP
ibs30_lonlat%>% 
  leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  addPolygons(opacity = 0.5, color = "red", weight = 2, fill = FALSE, 
              data = uml30_lonlat) %>%
  addPolygons(opacity = 1, color = "black", weight = 2, fill = FALSE) %>%
  addRasterImage(accu_lucpfip_suma, project = TRUE, colors = cb) %>% 
  addRasterImage(accu_lucpfip_kali, project = TRUE, colors = cb) %>% 
  addRasterImage(accu_lucpfip_papu, project = TRUE, colors = cb) %>% 
  addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
  addLegend(colors = legend_colors_uml, labels = label_uml, position = "bottomright") %>% 
  addLegend(pal = cb,  values = dom, bins = 5, opacity = 0.7,
            labFormat = labelFormat(suffix = " ha"),
            #labels = bin_labels, # this lign does not do anything
            title = "LUC from primary forest <br/> 
            to industrial plantations <br/> 
            within 900ha parcels <br/> 
            2001-2015 accumulated",
            position = "bottomright") 
  




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



















# 
# make_stats_fc2000 <- function(island, threshold){
# 
#   ### Extract forest cover 2000 pixels in different polygons
#   areas_fc2000 <- matrix(nrow = 3, 
#                          ncol = 3, 
#                          dimnames = list(c("10km_catchment_radius", "30km_catchment_radius", "50km_catchment_radius"),
#                                          c("area_in_island", "area_in_uml_cr", "area_in_ibs_cr")))
# 
#   thed_gfc_data <- brick(paste0(here("build/input/outcome_variables/gfc_data_"),island,"_",threshold,"th.tif"))
#   # select the forestcover layer (band 1)
#   fc2000 <- thed_gfc_data[[1]] 
#   # remove useless other stack of gfc layers
#   #rm(thed_gfc_data)  
#   
#   
#   ## compute area of forest cover on the island
#   # areas_fc2000[,"area_in_island"] <- raster::extract(fc2000, 
#   #                                                    island_sf[island_sf$shape_des == island,"geometry"],
#   #                                                    fun = sum, 
#   #                                                    na.rm = TRUE, 
#   #                                                    progress = "text")
#   
#   aoi_island <- island_sf[island_sf$shape_des == island,"geometry"] %>% as("Spatial")
#   #aoi_island$label <- island
#   
#   rasterOptions(progress = "text")
#   
#   # 60k seconds and a 9Gb tmp file (with chunksize 2e+9). 
#   stat_list <- gfc_stats(aoi = aoi_island, 
#                          thed_gfc_data, 
#                          dataset = "GFC-2018-v1.6")
#   saveRDS(stat_list, 
#           paste0(here("analysis/output/gfc_stats_output_"),island,"_",threshold,"th.Rdata"))
#   
#   print("area_in_island done")
#   print(Sys.time())
#   
#   CR <- 10
#   while(CR < 60){  
#     ## compute area of forest cover within catchment radius of all UML mills
#     areas_fc2000[paste0(CR,"km_catchment_radius"), "area_in_uml_cr"] <- raster::extract(fc2000,
#                                                                                         uml_cr[[CR]], 
#                                                                                         fun = sum, 
#                                                                                         na.rm = TRUE, 
#                                                                                         progress = "text")
#     
#     print(paste0(CR, "CR UML done"))
#     print(Sys.time())
#     
#     ## compute area of forest cover within catchment radius of IBS_UML mills only
#     areas_fc2000[paste0(CR,"km_catchment_radius"), "area_in_ibs_cr"] <- raster::extract(fc2000,
#                                                                                         ibs_cr[[CR]], 
#                                                                                         fun = sum, 
#                                                                                         na.rm = TRUE, 
#                                                                                         progress = "text")
#     print(paste0(CR, "CR IBS done"))
#     print(Sys.time())
#     
#     CR <- CR + 20
#     }
#   
#   rm(fc_2000, uml_cr, ibs_cr)
#   # at this point, each value in areas_fc2000 is the count of pixels covered with forest. 
#   # This value is converted to an area by approximating each pixel's area with a unique value (in m2): 27.7
#   # (because once projected to indonesian_crs, the gfc raster has resolution(27.8 ; 27.6))
#   areas_fc2000 <- areas_fc2000*27.7
#   
#   saveRDS(areas_fc2000, paste0(here("analysis/output/areas_fc2000_"),island, "_",threshold,"th.Rdata"))
#   
#   return(areas_fc2000)
# }
# 
# tic()
# make_stats_fc2000(island = "Sumatra", threshold = 30)
# toc()







### VISUALIZATION

# provinces <- cbind(provinces, st_coordinates(st_centroid(provinces)))
# ggplot(data = provinces) +
#   geom_sf() +
#   geom_sf(data = provinces, fill = NA) + 
#   geom_text(data = provinces, aes(X, Y, label = NAME_1))

## Map building

# info you want the map to display
#provinces$popup <- provinces$NAME_1
# MAP
# provinces %>% 
#   leaflet() %>% 
#   addTiles()%>%
#   addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
#   addPolygons(opacity = 0.5, color = "red", weight = 2, fill = FALSE, popup = ~provinces$NAME_1)#%>%
#   #addAwesomeMarkers(data = provinces, icon = icons, popup = ~provinces$NAME_1)
# 


## What we want is to get an integer of the area of grid cells = 1 in this fc2000 
## Two options, 
# 1. either compute cell specific area (but for now area(fc2000) returns a raster with only 0)
# 2. or approximate with a unique cell area. 

# 1. each cell has its area for value
# area(fc2000,
#      filename = "cell_area.tif", 
#      datatype = "INT2U", 
#      overwrite = TRUE)
# cell_area <- raster("cell_area.tif")
# # each cell has its area for value if it is covered with forest
# # cell_area could also be just replace by the average cell_area in the region (~30m2)
# rs <- stack(fc2000, cell_area)
# fc2000_area <- calc(fc2000, 
#                     fun = function(x){x*cell_area})


