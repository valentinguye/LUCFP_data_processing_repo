### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# This script adds variables to the parcel panel data frame. 
# In particular, it compute the number of UML mills each parcel can reach within 10, 30 and 50km. 
# And it adds geographic variables (island and districts). 
#
#
#   Inputs: UML most complete version
#           --> UML_valentin_imputed_est_year.dta
# 
#           parcel panel from previous step (i.e. making weighted averages)
#           --> pattern: wa_panel_parcels_ ; for each parcel_size and catchment_radius combination
#         
#           island polygons
#           --> temp_data/processed_indonesia_spatial/island_sf
#   
#           Province polygons
#           --> input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
#
#           district polygons and names (prepare in prepare_crosswalks.do)
# 
#           baseline forest extents variables in cross-section (prepare in prepare_2000_forest_extents.R)
#
#
#   Outputs: parcel panel with new columns: the parcel and time varying numbers of UML mills reachable within 10, 30 and 50km. 
#                                           the island, the district, and the pixelcounts and areas (ha) of 2000 forest extents
#                                           (for 30% tree canopy density outside indsutrial plantations and total primary forest. 
#           --> pattern temp_data/processed_parcels/parcels_panel_reachable_uml_
#                       temp_data/processed_parcels/parcels_panel_geovars_
#                       temp_data/processed_parcels/parcels_panel_final_ 
#                       ; for each parcel_size and catchment_radius combination
# 
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("plyr", "dplyr", "readstata13", "foreign", "sjmisc",
                   "rgdal", "sf", "nngeo",
                   "DataCombine", 
                   "parallel")
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

# 3. If the troubling packages could not be loaded ("there is no package called ...) 
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

### GRID CELL SIZE
parcel_size <- 3000

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


# The sample coverage has to be computed as the number of reachable IBS-UML matched sample
# related to the number of reachable UML mills, with the all the former being included in the latter. 
# We cannot compute a ratio of all mills in the analysis sample - i.e. also including IBS not matched with 
# UML but with a desa centroid - because we would not know whether a mill from the latter group 
# is an additional mill that is not in UML or if it is in UML but it was not matched. 

# read the sample panel of IBS geolocalized mills
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
ibsuml <- ibs[!duplicated(ibs$firm_id),]
# keep only those that are matched with UML. 
ibsuml <- ibsuml[ibsuml$uml_matched_sample==1,]

# make it a cross section
ibsuml <- ibsuml[!duplicated(ibsuml$firm_id),]
ibsuml <- ibsuml[!is.na(ibsuml$lat),]
ibsuml <- st_as_sf(ibsuml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
ibsuml <- st_transform(ibsuml, crs = indonesian_crs)

# read the most complete version of UML we have. 
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
uml <- uml[!is.na(uml$lat),]
uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
uml <- st_transform(uml, crs = indonesian_crs)


# so we do not do that: 
# # all mills are IBS-UML matched + IBS not matched with UML but with desa centroid)
# all_mills <- bind_rows(uml[,c("trase_code", "year","lon","lat")], 
#                        ibs[ibs$analysis_sample==1 & ibs$uml_matched_sample == 0,c("firm_id", "year", "lon", "lat")])
# 
# ibs[ibs$analysis_sample==1 & ibs$uml_matched_sample == 0,"firm_id"] %>% unique() %>% length()
# length(unique(uml$trase_code))
# nrow(unique(dplyr::select(all_mills, -year)))
# break it down to cross sections just like for IBS above. 
# class(all_mills$year)
# all_mills_cs <- lapply(years, FUN = function(x) all_mills[all_mills$year == x,]) 
# all_mills_cs <- lapply(all_mills_cs, FUN = st_as_sf, coords =  )
# all_mills_cs <- lapply(all_mills_cs, FUN = st_transform, crs = indonesian_crs)
# all_mills_cs <- lapply(all_mills_cs, FUN = st_geometry)
### ### ###



#catchment_radius <- 3e4
# This takes ~1h 
#### ADD N REACHABLE UML AND SAMPLE COVERAGE ####
make_n_reachable_uml <- function(parcel_size, catchment_radius){
  
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"km_IBS_CR.rds")))
  
  
  # make a spatial cross section of it (parcels' coordinates are constant over time)
  parcels_centro <- parcels[!duplicated(parcels$lonlat), c("lonlat", "idncrs_lat", "idncrs_lon")]
  # (lon lat are already expressed in indonesian crs)
  parcels_centro <- st_as_sf(parcels_centro, coords = c("idncrs_lon", "idncrs_lat"), remove = T, crs = indonesian_crs)
  
  parcels$newv_uml <- rep(0, nrow(parcels))
  parcels$newv_ibsuml <- rep(0, nrow(parcels))
  
  for(t in 1:length(years)){
    
    # UML
    # This is not a panel, so the information on presence or not a given year is whether 
    # the establishment year is anterior. We impute NA establishment year to be older than 1998. 
    present_uml <- uml[uml$est_year_imp <= years[t] | is.na(uml$est_year_imp),]
    
    annual_reachable_uml <- st_is_within_distance(parcels_centro, present_uml, dist = catchment_radius)
    parcels[parcels$year == years[t], "newv_uml"] <- lengths(annual_reachable_uml)
    
    # IBS-UML
    present_ibsuml <- ibsuml[ibsuml$est_year_imp <= years[t] | is.na(ibsuml$est_year_imp),]
    
    annual_reachable_ibsuml <- st_is_within_distance(parcels_centro, present_ibsuml, dist = catchment_radius)
    parcels[parcels$year == years[t], "newv_ibsuml"] <- lengths(annual_reachable_ibsuml)
    
    # Note that n_reachable_ibs was already computed in wa_at_parcels.R, 
    # but this includes ibs that are not matched with uml, which we do not want to count here. 
    # We used the year variable from the IBS panel to determine whether a firm was present in a given year. 
    # This means that when a firm has a yearly record missing, although we know it was there 
    # that year bc we observe an older and an earlier records, we do not count it as reachable 
    # by the parcel, because no IBS information would be usable that year in such a case.
    # But when the mill has a record line in IBS this year, but some or all information is missing we still 
    # count the mill as one more being reachable, altough we use no info from it.
    
  }
  
  
  # IBS/all mills (IBS-UML matched + IBS not matched with UML but with desa centroid) -> sample coverage
  # any reachable ibs is also counted as a reachable uml : YES if ibs is defined as ibs[ibs$uml_matched_sample==1,]
  nrow(parcels[parcels$newv_uml< parcels$newv_ibsuml,])
  
  # we give to each parcel a ratio that informs on the share of the total influence (from all possible mills known)
  # that is catched by our sample of analysis, i.e. geo-localized palm oil mills. 
  parcels$ratio <- rep(0, nrow(parcels))  
  parcels[parcels$newv_uml != 0, "ratio"] <- 100*(parcels[parcels$newv_uml != 0, "newv_ibsuml"]/parcels[parcels$newv_uml != 0,"newv_uml"])
  
  colnames(parcels)[colnames(parcels) == "newv_ibsuml"] <- paste0("n_reachable_ibsuml")
  colnames(parcels)[colnames(parcels) == "newv_uml"] <- paste0("n_reachable_uml")
  colnames(parcels)[colnames(parcels) == "ratio"] <- paste0("sample_coverage")
  
  
  # lag n reachable 
  parcels <- dplyr::arrange(parcels, lonlat, year)
  parcels <- DataCombine::slide(parcels,
                                Var = "n_reachable_uml", 
                                TimeVar = "year",
                                GroupVar = "lonlat",
                                NewVar = paste0("n_reachable_uml_lag1"),
                                slideBy = -1, 
                                keepInvalid = TRUE)
  parcels <- dplyr::arrange(parcels, lonlat, year)
  
  
  ### MAKE ID FOR THE SET OF REACHABLE IBS MILLS
  
  ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))  
  # keep only geolocalized mills
  ibs <- ibs[ibs$analysis_sample == 1,]
  ibs_cs <- ibs[!duplicated(ibs$firm_id),]
  ibs_cs <- st_as_sf(ibs_cs, coords = c("lon", "lat"), crs = 4326)
  ibs_cs <- st_transform(ibs_cs, crs = indonesian_crs)
  
  
  # make spatial
  # parcels_cs <- parcels[!duplicated(parcels$lonlat),c("lonlat", "year", "idncrs_lat", "idncrs_lon", "lon", "lat")]
  # parcels_cs <- st_as_sf(parcels_cs, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  # row.names(parcels_cs) <- parcels_cs$lonlat
  
  # define the set of mills within a distance. 
  sgbp <- st_is_within_distance(parcels_centro, ibs_cs, dist = catchment_radius, sparse = TRUE)
  
  usgbp <- unique(sgbp)
  
  parcels_centro$reachable <- rep(NA,nrow(parcels_centro))
  for(i in 1:length(usgbp)){
    parcels_centro$reachable[sgbp %in% usgbp[i]] <- i
  }
  
  
  # # Define the nearest mill 
  # nearest_mill_idx <- st_nearest_feature(parcels_centro, ibs_cs)
  # 
  # parcels_centro$nearest_firm_id <- ibs_cs$firm_id[nearest_mill_idx]
  
  parcels <- left_join(parcels, parcels_centro[,c("lonlat", "reachable")], by = "lonlat")#, "nearest_firm_id"
  
  parcels <- st_drop_geometry(parcels)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_reachable_uml_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
}


catchment_radius <- 10000
while(catchment_radius < 60000){
  
  make_n_reachable_uml(parcel_size, catchment_radius) 
  
  catchment_radius <- catchment_radius + 20000
}

rm(ibs, ibsuml)

#### ADD GEOGRAPHIC VARIABLES AND THEIR TRENDS ####

# read geographic shapefiles
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

# make the operation faster by using island bbox (other wise the island polygons make 
# the computation very long)
# (and this also includes parcel centroids in the sea)
island_sf_prj_bbox <- sapply(island_sf_prj$geometry, function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)

# province
province_sf <- st_read(file.path("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp"))
province_sf <- dplyr::select(province_sf, NAME_1)
province_sf_prj <- st_transform(province_sf, crs = indonesian_crs)

# district
district_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/district2000_sf"))
district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)

#sub-district
subdistrict <- st_read(file.path("input_data/indonesia_spatial/podes_bps2014"))
subdistrict_prj <- st_transform(subdistrict, crs = indonesian_crs)

catchment_radiuseS <- c(1e4, 3e4, 5e4)#
for(catchment_radius in catchment_radiuseS){
  
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"km_IBS_CR.rds")))
  
  parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  
  # Work with a cross section for province and district attribution
  parcels_cs <- parcels[!duplicated(parcels$lonlat),]
  
  ### ISLAND variable
  
  parcels$island <- rep("", nrow(parcels))
  
  sgbp <- st_within(parcels$geometry, island_sf_prj_bbox)
  # the bboxes of Sumatra and Kalimantan intersect a bit, so we check that no parcel falls 
  # in the intersection, this is the case for catchment_radius = 50km, for 2538 parcels, 
  # these parcel centroids belong to Kalimantan (after visual check). 
  # intersect <- st_intersection(island_sf_prj_bbox[1], island_sf_prj_bbox[3])
  # plot(island_sf_prj_bbox[[1]])
  # plot(island_sf_prj, add = TRUE)
  # plot(parcels$geometry[parcels$island==4], col = "red", add = TRUE)
  
  sgbp[lengths(sgbp)==2] <- 3
  
  # island_sf_prj features are in this order : 1 Sumatra; 2 Papua; 3 Kalimantan
  unique(unlist(sgbp))
  parcels$island <- unlist(sgbp)
  
  parcels$island <- replace(parcels$island, parcels$island == 1, "Sumatra")
  unique(parcels$island)
  parcels$island <- replace(parcels$island, parcels$island == 2, "Papua")
  parcels$island <- replace(parcels$island, parcels$island == 3, "Kalimantan")
  
  
  ### PROVINCE variable
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_prov_idx <- st_nearest_feature(parcels_cs, province_sf_prj)
  
  parcels_cs$province <- province_sf_prj$NAME_1[nearest_prov_idx]
  
  ### DISTRICT variable
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_dstr_idx <- st_nearest_feature(parcels_cs, district_sf_prj)
  
  parcels_cs$district <- district_sf_prj$name_[nearest_dstr_idx]
  # (4 parcels are closest to district with no name (NA) )
  
  
  
  ### SUB-DISTRICT AND VILLAGE   VARIABLES
  # use st_join, don't know why I did not use it above...
  parcels_cs <- st_join(parcels_cs, 
                        subdistrict_prj[,c("KECAMATAN", "DESA")], 
                        join = st_nearest_feature)
  
  
  ### NEIGHBOR VARIABLE
  # create a grouping variable at the cross section (9 is to recall the the group id includes the 8 neighbors + the central grid cell.)
  
  # identify neighbors
  # this definition of neighbors includes the 8 closest, surrounding, grid cells.
  parcels_buf <- st_buffer(parcels_cs, dist = parcel_size - 10)
  row.names(parcels_buf) <- parcels_buf$lonlat
  sgbp <- st_intersects(parcels_buf)
  
  neighbors <- list()
  length(neighbors) <- nrow(parcels_cs)
  parcels_cs$neighbors <- neighbors
  parcels_cs$neighbors <- lapply(1:nrow(parcels_cs), FUN = function(i){parcels_cs$lonlat[sgbp[[i]]]}) 
  
  rm(parcels_buf)
  
  ### MERGE WITH PANEL  
  parcels <- st_drop_geometry(parcels)
  parcels_cs <- st_drop_geometry(parcels_cs)
  
  parcels <- left_join(parcels, parcels_cs[,c("lonlat", "province", "district", "KECAMATAN", "DESA", "neighbors")], by = "lonlat")
  
  names(parcels)[names(parcels)=="KECAMATAN"] <- "subdistrict"
  names(parcels)[names(parcels)=="DESA"] <- "village"
  
  rm(parcels_cs)
  
  
  # SPATIAL TRENDS VARIABLES
  parcels$island_year <- paste0(parcels$island,"_",parcels$year)
  parcels$province_year <- paste0(parcels$province,"_",parcels$year)
  parcels$district_year <- paste0(parcels$district,"_",parcels$year)
  parcels$subdistrict_year <- paste0(parcels$subdistrict,"_",parcels$year)
  parcels$village_year <- paste0(parcels$village,"_",parcels$year)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
}
rm(district_sf, district_sf_prj, province_sf, province_sf_prj, island_sf, island_sf_prj, island_sf_prj_bbox)
# parcels_list[[match(catchment_radius, catchment_radiuseS)]] <- parcels



#### IDENTIFY "EXTENSIVE MARGIN" PARCELS ####
# In theory, mills establish only if a sufficient supply base is available. In practice, this is insured by companies developing plantations 
# in concomitance with opening a mill. In such an integrated development, there is no reason that mills and plantations are far appart. 
# Hence, we identify these plantations as being very close to mills. 
# Concretely, we identify all grid cells that are within

# read the most complete version of UML we have. 
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
uml <- uml[!is.na(uml$lat),]
uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
uml <- st_transform(uml, crs = indonesian_crs)

catchment_radiuseS <- c(1e4, 3e4, 5e4)#
for(catchment_radius in catchment_radiuseS){
  
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"km_IBS_CR.rds")))
  
  # make a spatial cross section of it (parcels' coordinates are constant over time)
  parcels_grid <- parcels[!duplicated(parcels$lonlat), c("lonlat", "idncrs_lat", "idncrs_lon")]
  # (lon lat are already expressed in indonesian crs)
  
  parcels_grid <- st_as_sf(parcels_grid, coords = c("idncrs_lon", "idncrs_lat"), remove = TRUE, crs = indonesian_crs)
  
  parcels_grid <- st_buffer(parcels_grid, dist = 1500)
  
  nn <- st_nn(uml, parcels_grid, sparse = TRUE, k = 4, maxdist = 3000, parallel = detectCores() - 1)
  
  length(nn[lengths(nn)>0]) # 952 mills have parcels from our sample that are within 3km and are among the 4 nearest ones.
  # the other ~200 mills are UML mills that are not in our area of interest and hence have no parcel within 3km. 
  unique(lengths(nn))
  
  # A parcel may be in the 4 nearest parcels from 2 mills. For now, we do not care to diferentiate which mill it should get the attributes of. 
  extensive_idx <- unique(unlist(nn)) # length is 3487
  length(unlist(nn)) # 3763 so not so many more
  
  # make the indicator variable
  parcels_grid$extensive <- rep(FALSE,nrow(parcels_grid))
  parcels_grid[extensive_idx,"extensive"] <-  TRUE
  
  
  ## make a more restrictive indicator variable, based on time 
  
  # for this exercise, we consider that UML mills with unknown establishement date were not established after 2000
  # for simplicity, we set the year to 2000 as we do not need to be more accurate for this exercise
  uml[is.na(uml$est_year_imp),"est_year_imp"] <- 2000
  
  parcels_grid$mill_est_year_imp <- rep(NA,nrow(parcels_grid))
  for(idx in extensive_idx){
    # identify the mill(s) that each extensive parcel is closer to
    mills <- uml[sapply(nn, FUN = function(nn_elm){idx %in% nn_elm}),]
    # extract the earliest establishment year 
    parcels_grid$mill_est_year_imp[idx] <- mills[mills$est_year_imp == min(mills$est_year_imp),]$est_year_imp
  }
  
  parcels_grid <- st_drop_geometry(parcels_grid)
  
  parcels <- left_join(parcels[,c("lonlat","year")], parcels_grid, by = "lonlat")
  
  # switch the extensive indicator from TRUE to FALSE for the annual records that are posterior to the close mill establishment.
  # in other words: let deforestation occurring very close to mills but after their establishment be counted as intensive margin. 
  parcels$extensive_restr <- parcels$extensive
  parcels$extensive_restr[parcels$year > parcels$mill_est_year_imp] <- FALSE
  summary(parcels$extensive_restr)
  summary(parcels$extensive)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_extmar_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
}


#### TIME DYNAMICS VARIABLES ####

### TIME SERIES 
ts <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
ts <- dplyr::select(ts, year, spread, dom_blwn_cpo_y)
# taxeffectiverate,
# ref_int_cpo_price,
# cif_rtdm_cpo,
# dom_blwn_cpo,
# fob_blwn_cpo,
# spread_int_dom_paspi,
# rho,
# dom_blwn_pko,
# cif_rtdm_pko,
# spread1, spread2, spread3, spread4, spread5, spread6)
# we only need the time series
ts <- ts[!duplicated(ts$year),]

catchment_radiuseS <- c(1e4, 3e4, 5e4)# 
for(catchment_radius in catchment_radiuseS){ 
  
  
  ### MERGE WEIGHTED AVERAGE AND NEAREST MILL VARIABLES
  
  wavars <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"CR.rds")))
  
  nmvars <- readRDS(file.path(paste0("temp_data/processed_parcels/nm_panel_parcels_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"CR.rds")))
  
  names(nmvars) %in% names(wavars)
  
  parcels <- inner_join(wavars, nmvars, by = c("lonlat","year"))
  
  rm(wavars, nmvars)
  
  ### EXPORT SHARES FROM PCT TO FRACTION   
  parcels$wa_prex_cpo_imp1 <- parcels$wa_prex_cpo_imp1/100 
  parcels$wa_prex_cpo_imp2 <- parcels$wa_prex_cpo_imp2/100 
  parcels$wa_lag1_prex_cpo_imp1 <- parcels$wa_lag1_prex_cpo_imp1/100 
  parcels$wa_lag1_prex_cpo_imp2 <- parcels$wa_lag1_prex_cpo_imp2/100 
  parcels$wa_avg_prex_cpo_imp1 <- parcels$wa_avg_prex_cpo_imp1/100 
  parcels$wa_avg_prex_cpo_imp2 <- parcels$wa_avg_prex_cpo_imp2/100 
  
  
  ### SHORT LAGS OF OTHER VARIABLES THAN PRICES 
  
  variables <- c("wa_prex_cpo_imp1","wa_prex_cpo_imp2",       
                 "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp",  "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",     
                 #"wa_concentration_10",     "wa_concentration_30", "wa_concentration_50",
                 "prex_cpo_imp1", "prex_cpo_imp2",
                 "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp",
                 "concentration_30", "concentration_50")
  
  for(voi in variables){
    ## lags
    parcels <- dplyr::arrange(parcels, lonlat, year)
    parcels <- DataCombine::slide(parcels,
                                  Var = voi, 
                                  TimeVar = "year",
                                  GroupVar = "lonlat",
                                  NewVar = paste0(voi,"_lag1"),
                                  slideBy = -1, 
                                  keepInvalid = TRUE)
    # parcels <- dplyr::arrange(parcels, lonlat, year)
  }
  
  #parcels1 <- parcels
  
  
  ### Operations relating contemporaneous to past information - on prices only
  variables <- c("wa_ffb_price_imp1", "wa_ffb_price_imp2", 
                 "wa_cpo_price_imp1", "wa_cpo_price_imp2",
                 "ffb_price_imp1", "ffb_price_imp2",
                 "cpo_price_imp1", "cpo_price_imp2", 
                 "wa_prex_cpo_imp1", "wa_prex_cpo_imp2") # those necessary for controls in IV lagged strategy
                #,"wa_pko_price_imp1", "wa_pko_price_imp2"
  
  for(voi in variables){
    
    ## short to long lags
    for(lag in c(1:5)){
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = voi, 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_lag",lag),
                                    slideBy = -lag, 
                                    keepInvalid = TRUE)
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      
    }
    # ## leads                               
    # for(lag in c(1:5)){
    #   parcels <- dplyr::arrange(parcels, lonlat, year)
    #   parcels <- DataCombine::slide(parcels,
    #                                 Var = voi, 
    #                                 TimeVar = "year",
    #                                 GroupVar = "lonlat",
    #                                 NewVar = paste0(voi,"_lead",lag),
    #                                 slideBy = lag, 
    #                                 keepInvalid = TRUE) 
    #   parcels <- dplyr::arrange(parcels, lonlat, year)
    # } 
    
    for(py in c(2,3,4)){
      ## Past-year averages (2, 3 and 4 years) - LONG RUN MEASURE - 
      parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
      parcels[is.nan(parcels$newv),"newv"] <- NA
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_",py,"pya")
      
      # lag the past year average (useful if we use those per se. as measures of LR price signal)
      # note that 3pya_lag1 is different from 4pya. 
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_",py,"pya"), 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_",py,"pya_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      
      
      # ## and absolute deviation - SHORT RUN MEASURE -
      # parcels <- mutate(parcels,
      #                   !!as.symbol(paste0(voi,"_dev_",py,"pya")) := !!as.symbol(paste0(voi)) - 
      #                     !!as.symbol(paste0(voi,"_",py,"pya")))
      # # # and relative deviation
      # # parcels <- mutate(parcels,
      # #                   !!as.symbol(paste0(voi,"_rdev_",py,"pya")) := (!!as.symbol(paste0(voi)) - 
      # #                                                               !!as.symbol(paste0(voi,"_",py,"pya"))) /
      # #  
      # 
      # # Lag these deviations by one year
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      # parcels <- DataCombine::slide(parcels,
      #                               Var = paste0(voi,"_dev_",py,"pya"), 
      #                               TimeVar = "year",
      #                               GroupVar = "lonlat",
      #                               NewVar = paste0(voi,"_dev_",py,"pya_lag1"),
      #                               slideBy = -1, 
      #                               keepInvalid = TRUE)  
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      
      ## and mean of contemporaneous and pya - OVERALL MEASURE - 
      
      # note that we add voi column (not lagged) in the row mean
      parcels$newv <- rowMeans(x = parcels[,c(voi, paste0(voi,"_lag",c(1:py)))], na.rm = FALSE)
      parcels[is.nan(parcels$newv),"newv"] <- NA
      # note that we name it ya (year average) and not past year average (pya). It the average of past years AND
      # contemporaneous obs..
      # When e.g. the looping variable py = 2, then pya are computed as averages of t-1 and t-2 values and
      # ya are computed as averages of t, t-1 and t-2 values. 
      # Coherently, the names of ya variables have _3ya_ (the py+1) to reflect the average being made over 3 years.    
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_",py+1,"ya")    
      
      # and lag it
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_",py+1,"ya"), 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_",py+1,"ya_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      
      
      ## Repeat for standard deviation
      # need have input parcels as a rowwise dataframe here (see https://dplyr.tidyverse.org/articles/rowwise.html)
      newv <- rowwise(parcels, c("lonlat", "year")) %>% 
        summarise(newv = sd(c_across(contains(paste0(voi,"_lag",c(1:py)))), na.rm = FALSE)) %>% 
        as.data.frame()
      parcels <- left_join(parcels, newv, by = c("lonlat", "year"))
      parcels[is.nan(parcels$newv),"newv"] <- NA
      
      # note y*v* for variation in the new name
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_",py+1,"yv")    
      
      # and lag it
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_",py+1,"yv"), 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_",py+1,"yv_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
      
    }
    
    
    
    
    ### Contemporaneous and past-year averaged YEAR-ON-YEAR GROWTH - DEPRECATED, USED NOWHERE AND THUS COMMENTED OUT
    
    # ## contemporaneous yoyg - SHORT RUN MEASURE - (invalid for at least the first record of each lonlat)
    # parcels <- mutate(parcels,
    #                   !!as.symbol(paste0(voi,"_yoyg")) := 100*(!!as.symbol(paste0(voi)) - 
    #                                                              !!as.symbol(paste0(voi,"_lag1"))) /
    #                     !!as.symbol(paste0(voi,"_lag1")))
    # 
    # ## Lagged yoyg 
    # # (the first lag is invalid for at least two first records of each lonlat;    
    # # the fourth lag is invalid for at least 5 first records)
    # for(lag in c(1:4)){
    #   parcels <- dplyr::arrange(parcels, lonlat, year)
    #   parcels <- DataCombine::slide(parcels,
    #                                 Var = paste0(voi,"_yoyg"), 
    #                                 TimeVar = "year",
    #                                 GroupVar = "lonlat",
    #                                 NewVar = paste0(voi,"_yoyg_lag",lag),
    #                                 slideBy = -lag, 
    #                                 keepInvalid = TRUE)  
    #   parcels <- dplyr::arrange(parcels, lonlat, year)
    # }
    
    
    ## past-years averaged yoyg - LONG RUN MEASURE - DEPRECATED, USED NOWHERE AND THUS COMMENTED OUT WITHIN THE LOOP
    for(py in c(2,3,4)){
      # parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_yoyg_lag",c(1:py))], na.rm = FALSE)
      # # treat NaNs that arise from means over only NAs when na.rm = T 
      # parcels[is.nan(parcels$newv),"newv"] <- NA
      # colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_yoyg_",py,"pya")
      # 
      # 
      # # Lag by one year
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      # parcels <- DataCombine::slide(parcels,
      #                               Var = paste0(voi,"_yoyg_",py,"pya"), 
      #                               TimeVar = "year",
      #                               GroupVar = "lonlat",
      #                               NewVar = paste0(voi,"_yoyg_",py,"pya_lag1"),
      #                               slideBy = -1, 
      #                               keepInvalid = TRUE)  
      # # parcels <- dplyr::arrange(parcels, lonlat, year)
      # 
      # 
      # ## contemporaneous AND pya yoyg mean - OVERALLMEASURE -   
      # 
      # # note that we add voi column (not lagged) in the row mean
      # parcels$newv <- rowMeans(x = parcels[,c(paste0(voi,"_yoyg"), paste0(voi,"_yoyg_lag",c(1:py)))], na.rm = FALSE)
      # parcels[is.nan(parcels$newv),"newv"] <- NA
      # # note that we name it ya (year average) and not past year average (pya). It the average of past years AND
      # # contemporaneous obs..
      # # When e.g. the looping variable py = 2, then pya are computed as averages of t-1 and t-2 values and
      # # ya are computed as averages of t, t-1 and t-2 values. 
      # # Coherently, the names of ya variables have _3ya_ (the py+1) to reflect the average being made over 3 years.    
      # colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_yoyg_",py+1,"ya")
      # 
      # 
      # # Lag by one year
      # parcels <- dplyr::arrange(parcels, lonlat, year)
      # parcels <- DataCombine::slide(parcels,
      #                               Var = paste0(voi,"_yoyg_",py+1,"ya"), 
      #                               TimeVar = "year",
      #                               GroupVar = "lonlat",
      #                               NewVar = paste0(voi,"_yoyg_",py+1,"ya_lag1"),
      #                               slideBy = -1, 
      #                               keepInvalid = TRUE)  
      # # parcels <- dplyr::arrange(parcels, lonlat, year)
      # 
    }
  }# closes the loop on variables  
  
  # remove some variables that were only temporarily necessary - NOT ANYMORE, IN ORDER TO TRY SPECIFICATION WITH ANNUAL PRICE SIGNALS 
  # vars_torm <- names(parcels)[grepl(pattern = "price_", x = names(parcels)) &
  #                               (grepl(pattern = "_lag2", x = names(parcels)) |
  #                                  grepl(pattern = "_lag3", x = names(parcels)) |
  #                                  grepl(pattern = "_lag4", x = names(parcels)) |
  #                                  grepl(pattern = "_lag5", x = names(parcels)))]
  # 
  # parcels <- parcels[,!(names(parcels) %in% vars_torm)]
  
  
  ### PRICE LOGARITHMS
  # Eventhough they are useful, we builod them within the regression function, to make the data sets a bit lighter. 
  # # select prices for which it's relevant/useful to compute the log
  # price_variables <- names(parcels)[grepl(pattern = "price_", x = names(parcels)) &
  #                                     !grepl(pattern = "pko", x = names(parcels)) &
  #                                     !grepl(pattern = "dev", x = names(parcels)) &
  #                                     !grepl(pattern = "yoyg", x = names(parcels))]
  # 
  # 
  # for(var in price_variables){
  #   parcels[,paste0("ln_",var)] <- log(parcels[,var])
  # }
  
  
  ### TIME SERIES & IV
  parcels <- merge(parcels, ts, by = "year")
  
  # Make the SHIFT SHARE INSTRUMENTAL VARIABLES 
  
  for(IMP in c(1,2)){
    # for(SP in c(1:6)){
    # make the instrument based on contemporaneous export shares
    parcels[,paste0("wa_iv_contemp_imp",IMP)] <- parcels[,paste0("wa_prex_cpo_imp",IMP)]*parcels[,paste0("spread")] + parcels[,"dom_blwn_cpo_y"]
    
    # make the instrument based on lagged export shares
    parcels[,paste0("wa_iv_lagged_imp",IMP)] <- parcels[,paste0("wa_lag1_prex_cpo_imp",IMP)]*parcels[,paste0("spread")] + parcels[,"dom_blwn_cpo_y"]
    
    # and on the average export share over time (at mill level)
    parcels[,paste0("wa_iv_avged_imp",IMP)] <- parcels[,paste0("wa_avg_prex_cpo_imp",IMP)]*parcels[,paste0("spread")] + parcels[,"dom_blwn_cpo_y"]
    
    # this is the average export shares of reachable mills over time at the parcel level -yields much less NA
    # not done because less representative of the spread shock at the mill level. 
    # avg_prex_cpo <- ddply(.data = parcels[, c("lonlat", "year", paste0("wa_prex_cpo_imp",IMP))], 
    #                       .variables = "lonlat", 
    #                       .fun = summarise, 
    #                       newvar := mean( !!as.symbol(paste0("wa_prex_cpo_imp",IMP)), na.rm = TRUE))
    # 
    # names(avg_prex_cpo) <- c("lonlat", paste0("wa_prex_cpo_imp",IMP,"_avg"))
    # 
    # parcels <- left_join(parcels, avg_prex_cpo, by = "lonlat")
    # 
    # # and make the instrument based on it
    # parcels[,paste0("iv_imp",IMP,"_avged")] <- parcels[,paste0("wa_prex_cpo_imp",IMP,"_avg")]*parcels[,paste0("spread")] 
    
    #}
  }
  
  # lag the iv variables
  # ivS <- c(paste0("iv",c(1:6),"_imp1"), paste0("iv",c(1:6),"_imp2"))
  ivS <- c(paste0("wa_iv_contemp_imp", c(1,2)), paste0("wa_iv_lagged_imp", c(1,2)), paste0("wa_iv_avged_imp", c(1,2)))
  
  parcels <- dplyr::arrange(parcels, lonlat, year)
  
  for(IV in ivS){
    for(lag in c(1:5)){
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = IV, 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(IV,"_lag", lag),
                                    slideBy = -lag, 
                                    keepInvalid = TRUE)  
      # parcels <- dplyr::arrange(parcels, lonlat, year)
    }
  }
  
  # View(parcels[!is.na(parcels$wa_prex_cpo_imp1_lag1) &
  #                parcels$year>2007 &
  #                parcels$wa_prex_cpo_imp1_lag1!=0 ,c("lonlat" ,"year", paste0("wa_prex_cpo_imp",c(1,2),"_lag1"),
  #                                                      paste0("spread",c(1:4)),
  #                                                      paste0("iv",c(1:4),"_imp1"),
  #                                                      paste0("iv",c(1:4),"_imp2"), 
  #                                                      paste0("iv",c(1:4),"_imp1_lag1"))])
  
  
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))  
  rm(parcels)
  
  # voi <- "wa_ffb_price_imp1"
  # View(parcels[500,c("lonlat", "year",
  #                    voi,#
  #                    paste0(voi,"_lag",c(1:4)),#
  #                    paste0(voi,"_3ya"),
  #                    paste0(voi,"_3ya_lag1"),
  #                    paste0(voi,"_4ya"),
  #                    paste0(voi,"_4ya_lag1"),
  #                    paste0(voi,"_5ya"),
  #                    paste0(voi,"_5ya_lag1"),
  #                    paste0(voi,"_3pya"),#
  #                    paste0(voi,"_3pya_lag1"),#
  #                    paste0(voi,"_dev_3pya"),#
  #                    paste0(voi,"_dev_3pya_lag1"),#
  #                    paste0(voi,"_yoyg"),#
  #                    paste0(voi,"_yoyg_lag", c(1:3)),#
  #                    paste0(voi,"_yoyg_4ya"),
  #                    paste0(voi,"_yoyg_4ya_lag1"),
  #                    paste0(voi,"_yoyg_3pya"),#
  #                    paste0(voi,"_yoyg_3pya_lag1"))])#
  
  
}
rm(ts)




##### ADD RSPO, CONCESSIONS, LEGAL LAND USE, & TIME SERIES - IV ##### 
### RSPO
rspo <- st_read("input_data/RSPO_supply_bases/RSPO-certified_oil_palm_supply_bases_in_Indonesia.shp")
rspo <- st_transform(rspo, crs = indonesian_crs)
# set up a year variable $
# we assume that the year of probable behavior change is the year after the date of the icletdate variable
# BUT THIS MAY CHANGE 
rspo$year <- sub(pattern = ".*, ", 
                 replacement = "", 
                 x = rspo$icletdate)
rspo$year[rspo$year=="no letter"] <- "2015"
rspo$year <- rspo$year %>% as.integer()
rspo$year <- rspo$year +1
# min(rspo$year)
# summary(rspo$year)
# rspo$icdate
# rspo$icletdate

### OIL PALM CONCESSIONS
cns <- st_read(file.path("input_data/oil_palm_concessions"))
cns <- st_transform(cns, crs = indonesian_crs)

### LEGAL LAND USE 
llu <- st_read(file.path("input_data/kawasan_hutan/Greenorb_Blog/final/KH-INDON-Final.shp"))
llu <- st_transform(llu, crs = indonesian_crs)
unique(llu$Fungsi)
names(llu)[names(llu) == "Fungsi"] <- "llu"

# restrict llu to provinces of interest
llu <- llu[llu$Province == "Sumatra Utara" |
             llu$Province == "Riau" |
             llu$Province == "Sumatra Selantan" |
             llu$Province == "Papua Barat" |
             llu$Province == "Kalimantan Timur" |
             llu$Province == "Kalimantan Selatan" |
             llu$Province == "Kalimantan Tengah" |
             llu$Province == "Kalimantan Barat" |
             llu$Province == "Bengkulu" |
             llu$Province == "Lampung" |
             llu$Province == "Jambi" |
             llu$Province == "Bangka Belitung" |
             llu$Province == "Kepuluan Riau" |
             llu$Province == "Sumatra Barat" |
             llu$Province == "Aceh", ]



catchment_radiuseS <- c(1e4, 3e4, 5e4)# 
for(catchment_radius in catchment_radiuseS){ 
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"km_IBS_CR.rds")))
  
  
  parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  
  
  ### RSPO
  parcels$rspo_cert <- rep(FALSE, nrow(parcels))
  
  for(y in min(rspo$year):max(parcels$year)){
    # select parcels from the given year
    parcels_cs <- parcels[parcels$year == y, c("lonlat", "year", "rspo_cert")]
    # select supply bases already certified this year
    rspo_cs <- rspo[rspo$year <= y,] %>% st_geometry()
    
    sgbp <- st_within(x = parcels_cs, y = rspo_cs)
    
    parcels$rspo_cert[parcels$year == y][lengths(sgbp) == 1] <- TRUE
  }
  
  
  ### OIL PALM CONCESSIONS
  # We do not observe whether a grid cell is within a concession annually. 
  # Therefore we only proceed with a cross section
  parcels_cs <- parcels[!duplicated(parcels$lonlat),]
  sgbp <- st_within(parcels_cs, cns)
  parcels_cs$concession <- rep(FALSE, nrow(parcels_cs))
  parcels_cs$concession[lengths(sgbp) > 0] <- TRUE
  
  parcels <- st_drop_geometry(parcels)
  parcels_cs <- st_drop_geometry(parcels_cs)
  
  parcels <- left_join(parcels, parcels_cs[,c("lonlat", "concession")], by = "lonlat")
  
  # note that some parcels fall within more than one concession record. There may be several reasons for concession overlaps 
  # like renewal of concession, with our withour aggrandisement. For our purpose, it only matters that there is at least one 
  # concession record. 
  
  
  ### LEGAL LAND USE 
  parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  # this is quite long (~5min)
  parcels_cs <- parcels[!duplicated(parcels$lonlat),]
  parcels_cs <- st_join(x = parcels_cs, 
                        y = st_make_valid(llu[,"llu"]), # st_make_valid bc thrown error otherwise 
                        join = st_within, 
                        left = TRUE)
  
  # some grid cells seem to fall within overlapping llu shapes though. 
  # It's really marginal (12 instances). Just remove the duplicates it produces.
  parcels_cs <- parcels_cs[!duplicated(parcels_cs$lonlat),]
  
  # merge back with panel 
  parcels <- st_drop_geometry(parcels)
  parcels_cs <- st_drop_geometry(parcels_cs)
  
  parcels <- left_join(parcels, parcels_cs[,c("lonlat", "llu")], by = "lonlat")
  
  unique(parcels$llu)
  ### ILLEGAL LUCFP 
  # one possible link to shed light on accronyms http://documents1.worldbank.org/curated/pt/561471468197386518/pdf/103486-WP-PUBLIC-DOC-107.pdf
  
  parcels <- dplyr::mutate(parcels,
                           illegal1 = (!concession & (llu != "HPK" )), # it's not in concession and not in a convertible forest zone. Don't add the following code, because these NAs are for all places outside the forest estate, and it changes exactly nothing to add this condition | llu == "<NA>"
                           illegal2 = (!concession & (llu == "HL" | # it's not in concession and it's in a permanent forest zone designation
                                                      
                                                      llu == "HP" | # production forest : " these areas may be selectively logged in a normal manner".
                                                      llu == "HPT" | # limited production forest : "These areas be logged less intensively than is permitted in the Permanent Production Forest" 
                                                        
                                                      llu=="HK" | # below are all categories of HK
                                                        llu=="KSA/KPA" |
                                                        llu=="KSA" | 
                                                        llu=="CA" | 
                                                        llu=="SM" | 
                                                        llu=="KPA" | 
                                                        llu=="TN" | 
                                                        llu=="TWA" | 
                                                        llu=="Tahura" | 
                                                        llu=="SML" | 
                                                        llu=="CAL" | 
                                                        llu=="TNL" | 
                                                        llu=="TWAL" | 
                                                        llu=="KSAL" | 
                                                        llu=="TB" | 
                                                        llu=="Hutan Cadangan")))
  
  # yields many missing in illegal because many grid cells are within a mising land use legal classification
  # parcels[!duplicated(parcels$lonlat) & !is.na(parcels$llu), c("lonlat", "concession", "llu", "illegal1", "illegal2")]
  
  
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_land_des_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
}

rm(cns, llu, parcels, parcels_cs, rspo, rspo_cs, sgbp)




#### MERGE ALL ADDITIONAL VARIABLES #### 

for(catchment_radius in c(1e4, 3e4, 5e4)){
  reachable <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_reachable_uml_",
                                        parcel_size/1000,"km_",
                                        catchment_radius/1000,"CR.rds")))
  
  geovars <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  
  extmar <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_extmar_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"CR.rds")))
  
  time_dyna <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                        parcel_size/1000,"km_",
                                        catchment_radius/1000,"CR.rds")))
  
  land_des <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_land_des_",
                                       parcel_size/1000,"km_",
                                       catchment_radius/1000,"CR.rds")))
  
  # the inner joins restrict the panel to only years in [2001, 2015]
  parcels <- inner_join(reachable[,c("lonlat","year",
                                     "n_reachable_uml", "n_reachable_uml_lag1","n_reachable_ibsuml","sample_coverage","reachable")], 
                        geovars[,c("lonlat","year",
                                   "island","province","district", "subdistrict", "village",
                                   "island_year","province_year","district_year", "subdistrict_year", "village_year")],
                        by = c("lonlat","year"))
  
  parcels <- inner_join(parcels, extmar[,c("lonlat","year", "extensive", "extensive_restr")], by = c("lonlat", "year"))
  
  parcels <- inner_join(parcels, time_dyna,
                        by = c("lonlat","year"))
  
  parcels <- inner_join(parcels, 
                        land_des[,c("lonlat","year",
                                    "rspo_cert", "concession", "llu", "illegal1", "illegal2")],
                        by = c("lonlat","year"))
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
  
  
}

rm(reachable, geovars, extmar, time_dyna, land_des, parcels)



#### CONSTRUCT THE OUTCOME VARIABLES - LHS - ####
# including the spatial lag - ~15 hours 
make_spatial_ov_lags <- function(catchment_radius){
  
  # merge lucpfip and lucfip outcome variables data sets together
  lucpfip <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                      parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  lucfip <- readRDS(file.path(paste0("temp_data/processed_parcels/lucfip_panel_",
                                     parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  lucpfsmp <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfsmp_panel_",
                                       parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  lucfsmp <- readRDS(file.path(paste0("temp_data/processed_parcels/lucfsmp_panel_",
                                      parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  lucpfip_dyn <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_dynamics_",
                                          parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  
  # keep only year before 2015 (after they mean nothing since we plantation data are from 2015)
  lucpfip <- lucpfip[lucpfip$year<=2015,] # now runs from 2001-1998
  lucfip <- lucfip[lucfip$year<=2015,] # now runs from 2001-1998
  lucpfsmp <- lucpfsmp[lucpfsmp$year<=2015,] # now runs from 2001-1998
  lucfsmp <- lucfsmp[lucfsmp$year<=2015,] # now runs from 2001-1998
  lucpfip_dyn <- lucpfip_dyn[lucpfip_dyn$year<=2015,]
  
  # # remove coordinates, except in one data frame
  lucfip <- dplyr::select(lucfip, -lat, -lon, -idncrs_lat, -idncrs_lon)
  lucpfsmp <- dplyr::select(lucpfsmp, -lat, -lon, -idncrs_lat, -idncrs_lon)
  lucfsmp <- dplyr::select(lucfsmp, -lat, -lon, -idncrs_lat, -idncrs_lon)
  lucpfip_dyn <- dplyr::select(lucpfip_dyn, -lat, -lon, -idncrs_lat, -idncrs_lon)
  
  # this is necessary because in lucpfip_dyn, coordinates got extracted slightly differently, after the 10th decimal or so, and we want them all equal 
  # to merge based on them.
  # lucpfip$lon <- round(lucpfip$lon, 6)
  # lucpfip$lat <- round(lucpfip$lat, 6)
  # lucfip$lon <- round(lucfip$lon, 6)
  # lucfip$lat <- round(lucfip$lat, 6)
  # lucpfsmp$lon <- round(lucpfsmp$lon, 6)
  # lucpfsmp$lat <- round(lucpfsmp$lat, 6)
  # lucpfip_dyn$lon <- round(lucpfip_dyn$lon, 6)
  # lucpfip_dyn$lat <- round(lucpfip_dyn$lat, 6)
  
  # if(nrow(lucpfip)!=nrow(lucfip) |
  # # nrow(lucpfsmp)!=nrow(lucfsmp) |
  # # nrow(lucpfip)!=nrow(lucfsmp) |
  # nrow(lucpfip)!=nrow(lucpfip_dyn)){print("LHS datasets don't all have the same number of rows, which is a known fact, see notes in script")} 
  # # lucfip has 3 more parcels (45 obs.)
  
  LHS <- inner_join(lucpfip, lucfip, by = c("lonlat", "year")) 
  LHS <- inner_join(LHS, lucpfsmp, by = c("lonlat", "year")) 
  LHS <- inner_join(LHS, lucfsmp, by = c("lonlat", "year"))   
  LHS <- inner_join(LHS, lucpfip_dyn, by = c("lonlat", "year"))# 
  if(nrow(LHS) != nrow(lucpfip)){stop("LHS datasets don't all have the same sets of parcels")}
  
  #   (nrow(lucpfip) - nrow(LHS))/15
  
  rm(lucpfip, lucpfip_dyn, lucpfsmp, lucfip, lucfsmp)
  
  ## make a variable that counts rapid and slow lucfp events (this is only computed for primary forest as of now)
  # THIS IS GOING TO BE THE MAIN OUTCOME VARIABLE INSTEAD OF lucpfip_pixelcount_total
  LHS$lucpfip_pixelcount <- LHS$lucpfip_rapid_pixelcount + LHS$lucpfip_slow_pixelcount
  LHS$lucfip_pixelcount <- LHS$lucfip_pixelcount_30th# pour l'instant on met lucfip_pixelcount_total dans "all producers" et pas rapid + slow, car on n'a 
  # pas calculé rapid et slow pour ce type de forêt encore
  
  # make variable that counts lucfp events on both small and medium sized plantations 
  LHS$lucpfsmp_pixelcount <- LHS$lucpfsp_pixelcount_total + LHS$lucpfmp_pixelcount_total
  LHS$lucfsmp_pixelcount <- LHS$lucfsp_pixelcount_30th + LHS$lucfmp_pixelcount_30th
  
  # The paragraph below is commented out because now the overlap of industrial and smallholders is handled in prepare_lucpfsmp.R
  #-------------------------
  # MAKE A UNIQUE DEPENDENT VARIABLE, WITH DUMMIES FOR ONLY INDUSTRIAL, ONLY SMALLHOLDERS, AND OVERLAPS
  
  # # This is true for obs. with positive lucpfip and lucpfsmp.
  # # It could be that both are positive while not overlapping (i.e. at the same place in the parcel), but we do not try to disentangle that here (would need to be done at the pixel level)
  # LHS$lucp_i_and_sm_bar <- (LHS$lucpfip_pixelcount == 0 | LHS$lucpfsmp_pixelcount == 0)
  # LHS$luc_i_and_sm_bar <- (LHS$lucfip_pixelcount == 0 | LHS$lucfsmp_pixelcount == 0)
  # 
  # ## make a variable that counts lucfp events on all types of plantations
  # LHS <- mutate(LHS, lucpfap_pixelcount = (lucpfip_pixelcount + lucpfsmp_pixelcount)*lucp_i_and_sm_bar) 
  # LHS <- mutate(LHS, lucfap_pixelcount = (lucfip_pixelcount + lucfsmp_pixelcount)*luc_i_and_sm_bar) 
  
  # overlap <- LHS[LHS$lucpfip_pixelcount > 0 & LHS$lucpfsmp_pixelcount > 0,]
  # overlap[overlap$lucpfip_pixelcount == overlap$lucpfsmp_pixelcount ,"lucpfip_pixelcount"] %>% length()
  # ---------------------------  
  
  ## make a variable that counts lucfp events on all types of plantations
  LHS <- mutate(LHS, lucpfap_pixelcount = (lucpfip_pixelcount + lucpfsmp_pixelcount))
  LHS <- mutate(LHS, lucfap_pixelcount = (lucfip_pixelcount + lucfsmp_pixelcount))
  
  
  #   # check that rapid + slow = total ? 
  #   LHS$lucpfip_rapidslowreplace_pixelcount <- LHS$lucpfip_rapid_pixelcount + LHS$lucpfip_slow_pixelcount + LHS$lucpfip_replace_pixelcount
  # 
  #   LHS$diff_0expct <- LHS$lucpfip_pixelcount_total - LHS$lucpfip_rapidslow_pixelcount
  #   
  #   LHS$diff_negexpct <- LHS$lucpfip_pixelcount_total - LHS$lucpfip_rapidslowreplace_pixelcount
  #   
  #   # investigate within cells that have at least one thing happening
  #   lhs <- LHS[LHS$lucpfip_pixelcount_total!=0 | 
  #                               LHS$lucpfip_rapid_pixelcount!=0 |
  #                               LHS$lucpfip_slow_pixelcount!=0 |
  #                               LHS$lucpfip_replace_pixelcount!=0 ,]
  # 
  #   Hmisc::describe(lhs$diff_negexpct)
  #   summary(lhs$diff_negexpct)
  # # this is never positif. Meaning that total is never greater than all three dynamic measures. 
  #   Hmisc::describe(lhs$diff_0expct)
  #   summary(lhs$diff_0expct)
  #   
  # # But some replace seems to be counted in total. This is when total is positive. This is 10% of the cases. 
  #   sum(lhs$diff_0expct>0)/nrow(lhs)
  #   sum(lhs$diff_0expct==0)/nrow(lhs)
  #   sum(lhs$diff_0expct<0)/nrow(lhs)
  #   
  #   View(lhs[lhs$diff_0expct<0,c("lonlat","year","diff_0expct", "lucpfip_pixelcount_total", "lucpfip_rapid_pixelcount",  "lucpfip_slow_pixelcount",  "lucpfip_replace_pixelcount" )])
  #   
  #   unique(no_0expct$year) # no special pattern in time
  #   
  #   # check pattern in space
  #   no_0expct <- st_as_sf(no_0expct, coords = c("lon", "lat"), crs = 4326)
  #   plot(st_geometry(no_0expct)) # no special pattern in space neither
  
  #   length(unique(LHS$diff_0expct))
  #   # 23442 obs. have a different lucpfip_pixelcount_total than rapid + slow, 
  #   # identify them   
  #   no_0expct <- LHS[LHS$diff_0expct != 0,]
  #   
  #   Hmisc::describe(LHS$diff_negexpct)
  #   
  #   
  #   LHS[LHS$diff ==0 &LHS$lucpfip_pixelcount_total!=0,c("lucpfip_rapid_pixelcount", "lucpfip_slow_pixelcount","lucpfip_pixelcount_total" )] %>% nrow()
  # 
  #   rm(lucpfip, lucfip, lucpfsmp, 
  #      #lucfsmp, 
  #      lucpfip_dyn)
  
  
  ## Compute the 4-year lagged total deforestation in neighboring cells.
  # (it does make sense only for the total deforestation to be spatially lagged, even when it's another type that is studied).
  
  # lag the outcome variable
  LHS <- dplyr::arrange(LHS, lonlat, year)
  LHS <- DataCombine::slide(LHS,
                            Var = "lucpfap_pixelcount",
                            TimeVar = "year",
                            GroupVar = "lonlat",
                            NewVar = "lucpfap_pixelcount_lag4",
                            slideBy = -4,
                            keepInvalid = TRUE)
  LHS <- dplyr::arrange(LHS, lonlat, year)
  
  # keep the most general cross section
  LHS_cs <- LHS[!duplicated(LHS$lonlat),c("lonlat", "year", "idncrs_lat", "idncrs_lon", "lucpfap_pixelcount")]
  
  # spatial
  LHS_cs <- st_as_sf(LHS_cs, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  
  # identify neighbors
  nn <- st_nn(st_geometry(LHS_cs), st_geometry(LHS_cs), sparse = TRUE, k = 9, maxdist = 10000)
  # remove self index
  nn <- lapply(nn, FUN = function(x){x[-1]})
  
  LHS_cs <- st_drop_geometry(LHS_cs)
  
  # get the parcel_id corresponding to the sgbp index
  lonlat_list <-  lapply(nn, function(ngbid){LHS_cs[ngbid, "lonlat"]})
  
  LHS[,"ngb_ov_lag4"] <- rep(NA, nrow(LHS))
  
  for(y in unique(LHS[LHS$year>2004,"year"])){
    LHS[LHS$year == y, "ngb_ov_lag4"] <- sapply(lonlat_list,
                                                FUN = function(lonlat_id){mean(LHS[LHS$lonlat %in% lonlat_id[[3]] & LHS$year == y, "lucpfap_pixelcount_lag4"],
                                                                               na.rm = TRUE)})
  }
  
  # empty neighbor sets (those that have no neighbors but themselves) and  end up with NaN as mean(integer(0))
  # turn them into NA
  LHS[is.nan(LHS$ngb_ov_lag4),"ngb_ov_lag4"] <- NA
  
  
  # save LHS, because the spatial operation is long
  saveRDS(LHS, paste0("temp_data/processed_parcels/parcels_lhs_panel_final_", 
                      parcel_size/1000,"km_",catchment_radius/1000,"CR.rds"))
  
}

### Execute
for(CR in c(3e4, 5e4)){
  make_spatial_ov_lags(catchment_radius = CR)
}



# 
# #### ADD N REACHABLE UML TO PARCELS IN ***UML*** CATCHMENT RADIUS #### 
# # this is a different thing, it does not add variables to our sample for analysis, but to 
# # another sample, that of parcels within CR of a **UML** mill, as outputed from prepare_lucpfip.R
# # this is necessary to later compute the aggregation factor in demand for deforestation. 
# 
# # prepare geographic data once before looping
# # island
# island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
# names(island_sf)[names(island_sf)=="island"] <- "shape_des"
# 
# island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)
# #province
# province_sf <- st_read(file.path("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp"))
# province_sf <- dplyr::select(province_sf, NAME_1)
# province_sf_prj <- st_transform(province_sf, crs = indonesian_crs)
# # district
# district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
# district_names <- read.dta13(file.path("temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta"))
# district_names <- district_names[!duplicated(district_names$bps_),]
# district_sf$d__2000 <- district_sf$d__2000 %>% as.character()
# district_names$bps_ <- district_names$bps_ %>% as.character()
# district_sf <- left_join(x = district_sf, y = district_names[,c("name_", "bps_")], 
#                          by = c("d__2000" = "bps_"),
#                          all = FALSE, all.x = FALSE, all.y = FALSE)
# 
# 
# district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)
# 
# catchment_radius <- 10000
# while(catchment_radius < 60000){
#   # read the panel of parcels within CR of a **UML** mill, as outputed from prepare_lucpfip.R
#   parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
#                                       parcel_size/1000,"km_",
#                                       catchment_radius/1000,"km_UML_CR.rds")))
#   
#   # make a spatial cross section of it (parcels' coordinates are constant over time)
#   parcels_centro <- parcels[parcels$year == 2001, c("lonlat", "idncrs_lat", "idncrs_lon")]
#   # (lon lat are already expressed in indonesian crs)
#   parcels_centro <- st_as_sf(parcels_centro, coords = c("idncrs_lon", "idncrs_lat"), remove = T, crs = indonesian_crs)
#   
#   parcels$newv_uml <- rep(0, nrow(parcels))
#   
#   for(t in 1:length(years)){
#     
#     # UML
#     # This is not a panel, so the information on presence or not a given year is whether 
#     # the establishment year is anterior. We impute NA establishment year to be older than 1998. 
#     present_uml <- uml[uml$est_year_imp <= years[t] | is.na(uml$est_year_imp),]
#     
#     annual_reachable_uml <- st_is_within_distance(parcels_centro, present_uml, dist = catchment_radius)
#     parcels[parcels$year == years[t], "newv_uml"] <- lengths(annual_reachable_uml)
#   } 
#   colnames(parcels)[colnames(parcels) == "newv_uml"] <- paste0("n_reachable_uml")
#   
#   ### ADD GEOGRAPHIC VARIABLES
#   parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
#   
#   # ISLAND variable
#   parcels$island <- rep("", nrow(parcels))
#   
#   # make the operation faster by using island bbox (other wise the island polygons make 
#   # the computation very long)
#   # (and this also includes parcel centroids in the sea)
#   island_sf_prj_bbox <- sapply(island_sf_prj$geometry, function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
#   
#   sgbp <- st_within(parcels$geometry, island_sf_prj_bbox)
#   # the bboxes of Sumatra and Kalimantan intersect a bit, so we check that no parcel falls 
#   # in the intersection, this is the case for catchment_radius = 50km, for 2538 parcels, 
#   # these parcel centroids belong to Kalimantan (after visual check). 
#   # intersect <- st_intersection(island_sf_prj_bbox[1], island_sf_prj_bbox[3])
#   # plot(island_sf_prj_bbox[[1]])
#   # plot(island_sf_prj, add = TRUE)
#   # plot(parcels$geometry[parcels$island==4], col = "red", add = TRUE)
#   
#   sgbp[lengths(sgbp)==2] <- 3
#   
#   # island_sf_prj features are in this order : 1 Sumatra; 2 Papua; 3 Kalimantan
#   unique(unlist(sgbp))
#   parcels$island <- unlist(sgbp)
#   
#   parcels$island <- replace(parcels$island, parcels$island == 1, "Sumatra")
#   unique(parcels$island)
#   parcels$island <- replace(parcels$island, parcels$island == 2, "Papua")
#   parcels$island <- replace(parcels$island, parcels$island == 3, "Kalimantan")
#   
#   
#   # PROVINCE variable
#   
#   # Work with a cross section for province and district attribution
#   parcels_cs <- parcels[!duplicated(parcels$lonlat),]
#   
#   
#   # the nearest feature function enables to also grab those parcels which centroids are in the sea.
#   nearest_prov_idx <- st_nearest_feature(parcels_cs, province_sf_prj)
#   
#   parcels_cs$province <- province_sf_prj$NAME_1[nearest_prov_idx]
#   
#   # DISTRICT variable
#   # the nearest feature function enables to also grab those parcels which centroids are in the sea.
#   nearest_dstr_idx <- st_nearest_feature(parcels_cs, district_sf_prj)
#   
#   # 4 parcels are closest to district with no name (NA) 
#   parcels_cs$district <- district_sf_prj$name_[nearest_dstr_idx]
#   
#   parcels <- merge(st_drop_geometry(parcels),
#                    st_drop_geometry(parcels_cs[,c("lonlat", "province", "district")]),
#                    by = "lonlat")
#   
#   # REGIONAL TRENDS VARIABLES
#   parcels$island_year <- paste0(parcels$island,"_",parcels$year)
#   parcels$province_year <- paste0(parcels$province,"_",parcels$year)
#   parcels$district_year <- paste0(parcels$district,"_",parcels$year)
#   
#   saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/lucpfip_panel_reachable_geovars_",
#                                     parcel_size/1000,"km_",
#                                     catchment_radius/1000,"km_UML_CR.rds")))
#   
#   
#   catchment_radius <- catchment_radius + 20000
# }
# 
#   















# # bon il y a ce truc bizarre que 33 mills ne sont les plus proches d'aucune parcelle, ce qui ne devrait pas arriver puisque toutes les mills ont au moins
# # une parcelle qui les recouvre. Cela dit c'est possible qu'une parcelle couvre deux mills et que le centroid soit plus proche de l'une que de l'autre... mais 33 fois ??? 
# # en fait c'est aussi toutes celles à Sulawesi, qui n'ont pas de parcelle qui les recouvre puisque on a sélectionné les parcelles par île.
# # ok ce qu'il se passe c'est que plusieurs mills ont les mêmes coordonnées (de village centroid probablement)
# nearest_mill_idx <- st_nearest_feature(parcels_cs, ibs_cs)
# ibs2 <- ibs_cs[!(ibs_cs$firm_id %in% ibs_cs$firm_id[unique(nearest_mill_idx)]),]
# ibs3 <- ibs2[ibs2$island_name == island,]
# st_crs(ibs3)
# st_within(ibs3, total_ca)
# 
# parcels2 <- parcels_cs[parcels_cs$lon > 97 &
#                          parcels_cs$lon < 98 &
#                          parcels_cs$lat > 1.4 &
#                          parcels_cs$lat < 1.5, ]
# nearest_test <- st_nearest_feature(parcels2, ibs_cs)
# plot(ibs_cs[262,"geometry"], add = TRUE, col = "green")
# ibs[ibs$firm_id==2072,]
# ibs[ibs$firm_id %in% ibs3$firm_id,c("lon","lat")] %>% unique()










#   dlm <- st_contains(district_sf_prj, parcels$geometry, sparse = FALSE)
#   
#   nrow(dlm)
#   class(dlm[,1])
#   row.names(dlm) <- district_sf_prj$name_
#   
#   # to each column of dlm, i.e. to each grid cell, attribute the name of the district 
#   # it is contained by. 
#   districts <- lapply(c(1:ncol(dlm)), FUN = function(col){row.names(dlm)[dlm[,col]==T]})
#   
#   districts[lengths(districts)==0] %>% length()
#   
#   # 6264, 101250 and (for 10 and 30 CR resp.) grid cells have their centroids in the sea. 
#   # We did not want to discard them because they may have informative patterns in the cell's 
#   # part that is on ground. 
#   
#   
#   
#   districts[lengths(districts)==0] <- "sea"
#   
#   parcels$district <- unlist(districts)
# 
#   nearest_dstr_idx <- st_nearest_feature(parcels[parcels$district=="sea", "geometry"], district_sf_prj)
#   # with this operation, some parcels get NA for an index of nearest district (from district_sf_prj)
#   # I don't know where it comes from. It is 72 records (4 parcels). We give them no name (NA). 
#   parcels1[parcels1$district=="sea" & is.na(nearest_dstr_idx1), "district"] <- ""
# 
#   parcels1$district[parcels1$district=="sea"][!is.na(nearest_dstr_idx1)]  <- district_sf_prj$name_[nearest_dstr_idx1][!is.na(nearest_dstr_idx1)]  
# all.equal(district_sf_prj$name_[nearest_dstr_idx1][!is.na(nearest_dstr_idx1)] , district_sf_prj$name_[!is.na(nearest_dstr_idx1)] )
#   
#   nearest_dstr_idx1[is.na(district_sf_prj$name_[nearest_dstr_idx1])]
#   
#   parcels1[parcels1$seadis == "Kab. Aceh Barat", "geometry"] %>% plot()
#   parcels1[parcels1$seadis == "sea" & parcels1$district == "Kab. Aceh Barat", "geometry"] %>% plot(col = "red", add = T)
#   district_sf_prj[district_sf_prj$name_=="Kab. Aceh Barat",]$geometry %>% plot(add =T)