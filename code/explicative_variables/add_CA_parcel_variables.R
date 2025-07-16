### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# This script adds variables to the parcel panel data frame. 
# In particular, it compute the number of UML mills each parcel can reach within 10, 30 and 50km. 
# And it adds geographic variables (provinces and districts). 
#
#
#   Inputs: UML most complete version
#           --> UML_valentin_imputed_est_year.dta
# 
#           parcel panel from previous step (i.e. making weighted averages)
#           --> pattern: wa_panel_parcels_ ; for each parcel_size and catchment_radius combination
#   
#           Province polygons
#           --> input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
#
#           district polygons and names (prepare in prepare_crosswalks.do)
#           --> input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp
#           --> temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
# 
#           baseline forest extents variables in cross-section (prepare in prepare_2000_forest_extents.R)
#           --> pattern: temp_data/processed_parcels/baseline_fc_cs_", for each parcel_size and IBS catchment_radius combination 
#
#
#   Outputs: parcel panel with new columns: the parcel and time varying numbers of UML mills reachable within 10, 30 and 50km. 
#                                           the province, the district, and the pixelcounts and areas (ha) of 2000 forest extents
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
neededPackages = c("plyr", "tidyr","dplyr", "readstata13", "foreign", "sjmisc",
                   "rgdal", "sf", 
                   "DataCombine")
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

# WE DO NOT COMPUTE REACHABLE UML HERE BECAUSE WE ARE NOT USING IT IN THE ANALYSIS SO FAR, NOR SAMPLE COVERAGE, 
# AND IT TAKES QUITE A WHILE TO COMPUTE BECAUSE THE ANNUAL OSRM DURATONS TO UML MILLS WOULD NEED TO BE REQUESTED. 


### PREPARE SPATIAL DATA 
# island
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"
island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

# province
province_sf <- st_read(file.path("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp"))
province_sf <- dplyr::select(province_sf, NAME_1)
province_sf_prj <- st_transform(province_sf, crs = indonesian_crs)

# district
district_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/district2000_sf"))
district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)

### PREPARE IBS AND UML DATA
# # read the sample panel of IBS geolocalized mills
# ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# # keep only those that are matched with UML. 
# ibsuml <- ibs[ibs$uml_matched_sample==1,]
# # make it a cross section
# ibsuml <- ibsuml[!duplicated(ibsuml$firm_id),]
# ibsuml <- ibsuml[!is.na(ibsuml$lat),]
# ibsuml <- st_as_sf(ibsuml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
# ibsuml <- st_transform(ibsuml, crs = indonesian_crs)

# read the most complete version of UML we have. 
# uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
# uml <- uml[!is.na(uml$lat),]
# uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
# uml <- st_transform(uml, crs = indonesian_crs)
# # there is no island column in this dataset, hence we select mills on the specific island geographically
# select_within <- st_within(x = uml, y = island_sf_prj[island_sf_prj$shape_des == "Sumatra",])
# uml_sumatra <- uml[lengths(select_within)>0,]
# select_within <- st_within(x = uml, y = island_sf_prj[island_sf_prj$shape_des == "Kalimantan",])
# uml_kalimantan <- uml[lengths(select_within)>0,]
# rm(select_within)


# island <- "Sumatra"
# travel_time <- 2
# t <- 8


#### ADD N_REACHABLE_UML AND GEOGRAPHIC VARIABLES AND THEIR TRENDS and QUEEN NEIGHBORS #### 

# this needs to be done at the island level, because it involves OSRM durations, that are relevant to compute at this level. 
for(travel_time in c(2)){# ,4,6
  
  # as the computations on n_reachable_uml requires using parcels at the level of island, we store the results in a list to 
  # be able to append parcels over islands after. 
  isl_parcel_list <- list()
  
  for(island in c("Sumatra", "Kalimantan")){
    
    # read the parcel panel as outputted from wa_at_parcels_durations
    parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                        island,"_",
                                        parcel_size/1000,"km_",
                                        travel_time,"h_CA.rds")))
    
    parcels$n_reachable_uml <- rep(NA, nrow(parcels))
    parcels$set_reachable_ibs <- rep(NA, nrow(parcels))
    parcels <- mutate(parcels, lonlat = paste0(lon, lat))
    
    # loop over years to get the annual number of reachable uml mills
    years <- seq(from = 1998, to = 2015, by = 1)
    for(t in 1:length(years)){
      
      # take the annual cross section of it (parcels' coordinates are constant over time)
      parcels_centro <- parcels[parcels$year == years[t], c("year", "lonlat", "lat", "lon")]
      # and of the uml data set (only useful for the checks)
      # if(island == "Sumatra"){
      #   uml_cs <- uml_sumatra[uml_sumatra$est_year_imp <= years[t] | is.na(uml_sumatra$est_year_imp),]
      # }
      # if(island == "Kalimantan"){
      #   uml_cs <- uml_kalimantan[uml_kalimantan$est_year_imp <= years[t] | is.na(uml_kalimantan$est_year_imp),]
      # }
      
      # read the duration matrix. 
      # This is the driving travel time (duration) matrix between each pair of parcel and mill. 
      # no matter the t, this matrix has the same amount of rows (parcels_centro) as parcels_centro
      dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_",travel_time,"h_UML_",years[t])))
      dur_mat <- dur_mat$durations
      
      # below are several checks that we are allocating the right duration to the right pair of parcel-mill. 
      # if(nrow(uml_cs) != ncol(dur_mat)){stop(paste0("duration matrix and mill cross section don't match in ",island,"_",parcel_size/1000,"km_",travel_time,"h_UML_",years[t]))}
      
      # for more safety, merge the dur_mat with the parcels_centro based on coordinates identifiers
      dur_mat <- as.data.frame(dur_mat)
      dur_mat$lonlat <- row.names(dur_mat)
      
      # sort = FALSE is MEGA IMPORTANT because otherwise dur_mat get sorted in a different way than parcels_centro
      dur_mat <- merge(dur_mat, parcels_centro[,c("lonlat", "lon", "lat")], by = "lonlat", all = FALSE, sort = FALSE)
      row.names(dur_mat) <- dur_mat$lonlat
      # so there should be the same amount of parcels_centro which coordinates matched, as the total amount of parcels.
      if(nrow(dur_mat) != nrow(parcels_centro)){stop(paste0("duration matrix and parcel cross section did not merge in ",island,"_",parcel_size/1000,"km_",travel_time,"h_UML_",years[t]))}
      
      # # additional check:
      row.names(parcels_centro) <- parcels_centro$lonlat
      
      parcels_centro <- dplyr::arrange(parcels_centro, lonlat)
      dur_mat <- dplyr::arrange(dur_mat, lonlat)
      
      if(!(all.equal(row.names(parcels_centro), row.names(dur_mat)))){stop()}

      # convert durations (in minutes) into logical whether each parcel-mill pair's travel time is inferior to a threshold travel_time
      dur_mat <- dplyr::select(dur_mat, -lonlat, -lon, -lat)
      dur_mat <- as.matrix(dur_mat)
      
      dur_mat_log <- dur_mat/(60) < travel_time
      
      # replace the few NAs with FALSE
      #anyNA(dur_mat_log)
      # dur_df_log = dur_mat_log %>% as.data.frame() 
      # dur_df_log <- tidyr::replace_na(dur_df_log, replace = FALSE)
      # # dur_mat_log <- tidyr::replace_na(dur_mat_log, replace = FALSE)
      # dur_mat_log <- mutate(dur_mat_log, if_else( = FALSE)
      
      #anyNA(dur_mat_log)
      
      # count reachable mills from each parcel  
      parcels_centro$n_reachable_uml <- sapply(1:nrow(parcels_centro), FUN = function(i){sum(dur_mat_log[i,], na.rm = TRUE)})
      
      # rearrange the panel, so that within cross sections, it's ordred by lonlat. 
      #parcels <- dplyr::arrange(parcels, year, lonlat)
      
      ### PREPARE THE SETS OF REACHABLE IBS MILLS in duration terms
      dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_",years[t])))
      dur_mat <- dur_mat$durations
      dur_mat <- as.data.frame(dur_mat)
      dur_mat$lonlat <- row.names(dur_mat)      
      
      dur_mat <- merge(dur_mat, parcels_centro[,c("lonlat", "lon", "lat")], by = "lonlat", all = FALSE, sort = FALSE)
      row.names(dur_mat) <- dur_mat$lonlat
      # so there should be the same amount of parcels_centro which coordinates matched, as the total amount of parcels.
      if(nrow(dur_mat) != nrow(parcels_centro)){stop(paste0("duration matrix and parcel cross section did not merge in ",island,"_",parcel_size/1000,"km_",travel_time,"h_UML_",years[t]))}
      
      row.names(parcels_centro) <- parcels_centro$lonlat
      parcels_centro <- dplyr::arrange(parcels_centro, lonlat)
      dur_mat <- dplyr::arrange(dur_mat, lonlat)
      
      # convert durations (in minutes) into logical whether each parcel-mill pair's travel time is inferior to a threshold travel_time
      dur_mat <- dplyr::select(dur_mat, -lonlat, -lon, -lat)
      dur_mat <- as.matrix(dur_mat)
      dur_mat_log <- dur_mat/(60) < travel_time
      
      # Collect column names of columns with TRUE
      dur_mat_log %>% colnames() %>% unique() %>% length() == ncol(dur_mat_log)
      
      # parcels_centro$set_reachable = list()
      
      parcels_centro$set_reachable_ibs <- lapply(1:nrow(parcels_centro), 
                                             FUN = function(i){
                                               dur_mat_log[i,dur_mat_log[i,]] %>% names() %>% na.omit() %>% sort()
      })
      
      # MERGE BACK
      # this makes sure that the n_reachable_uml cross section is integrated into the panel in face of the right lonlat. 
      parcels[parcels$year==years[t], c("lonlat", "year", "n_reachable_uml", "set_reachable_ibs")] <- 
        inner_join(parcels[, c("lonlat", "year")], 
                   parcels_centro[,c("lonlat", "year", "n_reachable_uml", "set_reachable_ibs")], 
                   by = c("lonlat", "year"))
      
      rm(parcels_centro, uml_cs, dur_mat, dur_mat_log)
    }# closes loop over years
    
    # set to NA where no mill is reachable  
    parcels$set_reachable_ibs[lengths(parcels$set_reachable_ibs)==0] <- NA
    # MAKE THE SET OF REACHABLE MILLS ID 
    parcels = 
      parcels %>% 
      mutate(reachable = match(set_reachable_ibs, unique(parcels$set_reachable_ibs)), 
             reachable = if_else(is.na(set_reachable_ibs), NA, reachable))
    
    # make the island variable now
    parcels$island <- island

    # store these parcels of a specific island in the dedicated list
    isl_parcel_list[[match(island, c("Sumatra", "Kalimantan"))]] <- parcels
    
    rm(parcels)
  }# closes loop over islands.
  
  # make the object of all parcels across the two islands
  parcels <- bind_rows(isl_parcel_list)
  
  
  
  
  ### ADD GEOGRAPHIC VARIABLES AND THEIR TRENDS 
  
  parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  
  
  ### ISLAND variable
  
  # Already set at the end of wa_at_parcels_durations.R
  
  ### PROVINCE variable
  
  # Work with a cross section for province and district attribution
  parcels_cs <- parcels[!duplicated(parcels$lonlat),]
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_prov_idx <- st_nearest_feature(parcels_cs, province_sf_prj)
  
  parcels_cs$province <- province_sf_prj$NAME_1[nearest_prov_idx]
  
  ### DISTRICT variable
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_dstr_idx <- st_nearest_feature(parcels_cs, district_sf_prj)
  
  # 4 parcels are closest to district with no name (NA) 
  parcels_cs$district <- district_sf_prj$name_[nearest_dstr_idx]
  
  parcels <- merge(st_drop_geometry(parcels),
                   st_drop_geometry(parcels_cs[,c("lonlat", "province", "district")]),
                   by = "lonlat")
  
  
  # SPATIAL TRENDS VARIABLES
  parcels$island_year <- paste0(parcels$island,"_",parcels$year)
  parcels$province_year <- paste0(parcels$province,"_",parcels$year)
  parcels$district_year <- paste0(parcels$district,"_",parcels$year)
  
  
  
  ### NEIGHBORS VARIABLE
  # create a grouping variable at the cross section (9 is to recall the the group id includes the 8 neighbors + the central grid cell.)
  parcels_cs <- parcels[!duplicated(parcels$lonlat),c("lonlat", "year", "idncrs_lat", "idncrs_lon")]
  
  # spatial
  parcels_cs <- st_as_sf(parcels_cs, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  
  # identify neighbors
  # this definition of neighbors includes the 8 closest, surrounding, grid cells.
  parcels_buf <- st_buffer(parcels_cs, dist = parcel_size - 10)
  row.names(parcels_buf) <- parcels_buf$lonlat
  sgbp <- st_intersects(parcels_buf)
  
  neighbors <- list()
  length(neighbors) <- nrow(parcels_cs)
  parcels_cs$neighbors <- neighbors
  parcels_cs$neighbors <- lapply(1:nrow(parcels_cs), FUN = function(i){parcels_cs$lonlat[sgbp[[i]]]}) 
  
  parcels <- left_join(parcels, parcels_cs[,c("lonlat", "neighbors")], by = "lonlat")
  
  
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_CA.rds")))
  
  rm(parcels, parcels_cs)
}





#### TIME DYNAMICS VARIABLES ####
for(travel_time in c(2)){ #2,4,4,6
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                      parcel_size/1000,"km_",
                                      travel_time,"h_CA.rds")))
  
  
  ### EXPORT SHARES FROM PCT TO FRACTION   
  parcels$wa_prex_cpo_imp1 <- parcels$wa_prex_cpo_imp1/100 
  parcels$wa_prex_cpo_imp2 <- parcels$wa_prex_cpo_imp2/100 
  
  
  ### Short lags of other variables than prices 
  
  variables <- c("wa_prex_cpo_imp1","wa_prex_cpo_imp2",       
                 "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp",  "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",     
                 #"wa_concentration_10",     "wa_concentration_30", "wa_concentration_50",     
                 "n_reachable_uml")#"n_reachable_ibs", "n_reachable_ibsuml",   "sample_coverage"
  
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
    parcels <- dplyr::arrange(parcels, lonlat, year)
  }
  
  #parcels1 <- parcels
  
  
  ### Operations relating contemporaneous to past information - on prices only
  variables <- c("wa_ffb_price_imp1", "wa_ffb_price_imp2", 
                 "wa_cpo_price_imp1", "wa_cpo_price_imp2") #,"wa_pko_price_imp1", "wa_pko_price_imp2"
  
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
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
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
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
      
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
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
      
    }
    
    
    
    
    ### Contemporaneous and past-year averaged YEAR-ON-YEAR GROWTH
    
    ## contemporaneous yoyg - SHORT RUN MEASURE - (invalid for at least the first record of each lonlat)
    parcels <- mutate(parcels,
                      !!as.symbol(paste0(voi,"_yoyg")) := 100*(!!as.symbol(paste0(voi)) - 
                                                                 !!as.symbol(paste0(voi,"_lag1"))) /
                        !!as.symbol(paste0(voi,"_lag1")))
    
    ## Lagged yoyg 
    # (the first lag is invalid for at least two first records of each lonlat;    
    # the fourth lag is invalid for at least 5 first records)
    for(lag in c(1:4)){
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_yoyg"), 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_yoyg_lag",lag),
                                    slideBy = -lag, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, lonlat, year)
    }
    
    
    ## past-years averaged yoyg - LONG RUN MEASURE - 
    for(py in c(2,3,4)){
      parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_yoyg_lag",c(1:py))], na.rm = FALSE)
      # treat NaNs that arise from means over only NAs when na.rm = T 
      parcels[is.nan(parcels$newv),"newv"] <- NA
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_yoyg_",py,"pya")
      
      
      # Lag by one year
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_yoyg_",py,"pya"), 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_yoyg_",py,"pya_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
      
      ## contemporaneous AND pya yoyg mean - OVERALLMEASURE -   
      
      # note that we add voi column (not lagged) in the row mean
      parcels$newv <- rowMeans(x = parcels[,c(paste0(voi,"_yoyg"), paste0(voi,"_yoyg_lag",c(1:py)))], na.rm = FALSE)
      parcels[is.nan(parcels$newv),"newv"] <- NA
      # note that we name it ya (year average) and not past year average (pya). It the average of past years AND
      # contemporaneous obs..
      # When e.g. the looping variable py = 2, then pya are computed as averages of t-1 and t-2 values and
      # ya are computed as averages of t, t-1 and t-2 values. 
      # Coherently, the names of ya variables have _3ya_ (the py+1) to reflect the average being made over 3 years.    
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_yoyg_",py+1,"ya")
      
      
      # Lag by one year
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_yoyg_",py+1,"ya"), 
                                    TimeVar = "year",
                                    GroupVar = "lonlat",
                                    NewVar = paste0(voi,"_yoyg_",py+1,"ya_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
    }
  }# closes the loop on variables  
  
  # remove some variables that were only temporarily necessary
  vars_torm <- names(parcels)[grepl(pattern = "price_", x = names(parcels)) &
                                (grepl(pattern = "_lag2", x = names(parcels)) |
                                   grepl(pattern = "_lag3", x = names(parcels)) |
                                   grepl(pattern = "_lag4", x = names(parcels)) |
                                   grepl(pattern = "_lag5", x = names(parcels)))]
  
  parcels <- parcels[,!(names(parcels) %in% vars_torm)]
  
  
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
  
  
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_CA.rds")))  
  rm(parcels)
} # closes the loop on travel_time


##### ADD RSPO, CONCESSIONS, LEGAL LAND USE, & TIME SERIES - IV ##### 


for(travel_time in c(2)){ #, 4, 6 is too heavy for now
  # load data in the loop, not time efficient but here the constraint is rather the memory so we want to be able to remove objects asap
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
               #llu$Province == "Papua Barat" |
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
  
  
  ### TIME SERIES 
  ts <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
  ts <- dplyr::select(ts, year, 
                      taxeffectiverate,
                      ref_int_cpo_price,
                      cif_rtdm_cpo,
                      dom_blwn_cpo,
                      fob_blwn_cpo,
                      spread_int_dom_paspi,
                      rho,
                      dom_blwn_pko,
                      cif_rtdm_pko,
                      spread1, spread2, spread3, spread4, spread5, spread6)
  # we only need the time series
  ts <- ts[!duplicated(ts$year),]
  
  
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                      parcel_size/1000,"km_",
                                      travel_time,"h_CA.rds")))
  
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
    
    rm(parcels_cs, rspo_cs, sgbp)
  }
  rm(rspo)
  
  ### OIL PALM CONCESSIONS
  
  # We do not observe whether a grid cell is within a concession annually. 
  # Therefore we only proceed with a cross section
  parcels_cs <- parcels[!duplicated(parcels$lonlat),]
  sgbp <- st_within(parcels_cs, cns)
  parcels_cs$concession <- rep(FALSE, nrow(parcels_cs))
  parcels_cs$concession[lengths(sgbp) > 0] <- TRUE
  rm(sgbp)
  parcels <- st_drop_geometry(parcels)
  parcels_cs <- st_drop_geometry(parcels_cs)
  
  parcels <- left_join(parcels, parcels_cs[,c("lonlat", "concession")], by = "lonlat")
  rm(cns, parcels_cs)
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
  rm(parcels_cs)
  
  #unique(parcels$llu)
  ### ILLEGAL LUCFP 
  # one possible link to shed light on accronyms http://documents1.worldbank.org/curated/pt/561471468197386518/pdf/103486-WP-PUBLIC-DOC-107.pdf
  
  parcels <- dplyr::mutate(parcels,
                           illegal1 = (!concession & (llu != "HPK" | llu == "<NA>")), # it's not in concession and not in a convertible forest zone
                           illegal2 = (!concession & (llu == "KSA/KPA" | # it's not in concession and it's in a permanent forest zone designation
                                                        llu == "KSA" |
                                                        llu == "KPA" |
                                                        llu == "KSAL" |
                                                        llu == "HP" |
                                                        llu == "HPT" |
                                                        llu == "HL")))
  rm(llu)
  # yields many missing in illegal because many grid cells are within a mising land use legal classification
  # parcels[!duplicated(parcels$lonlat) & !is.na(parcels$llu), c("lonlat", "concession", "llu", "illegal1", "illegal2")]
  
  
  
  ### TIME SERIES & IV
  # parcels <- merge(parcels, ts, by = "year")
  # 
  # # Make the SHIFT SHARE INSTRUMENTAL VARIABLES 
  # for(IMP in c(1,2)){
  #   for(SP in c(1:6)){
  #     parcels[,paste0("iv",SP,"_imp",IMP)] <- parcels[,paste0("wa_prex_cpo_imp",IMP,"_lag1")]*parcels[,paste0("spread",SP)] 
  #   }
  # }
  # rm(IMP, SP)
  # # lag the iv variables
  # ivS <- c(paste0("iv",c(1:6),"_imp1"), paste0("iv",c(1:6),"_imp2"))
  # 
  # for(IV in ivS){
  #   parcels <- dplyr::arrange(parcels, lonlat, year)
  #   parcels <- DataCombine::slide(parcels,
  #                                 Var = IV, 
  #                                 TimeVar = "year",
  #                                 GroupVar = "lonlat",
  #                                 NewVar = paste0(IV,"_lag1"),
  #                                 slideBy = -1, 
  #                                 keepInvalid = TRUE)  
  #   parcels <- dplyr::arrange(parcels, lonlat, year)
  # }
  # rm(IV, ivS)
  # # View(parcels[!is.na(parcels$wa_prex_cpo_imp1_lag1) &
  # #                parcels$year>2007 &
  # #                parcels$wa_prex_cpo_imp1_lag1!=0 ,c("lonlat" ,"year", paste0("wa_prex_cpo_imp",c(1,2),"_lag1"),
  # #                                                      paste0("spread",c(1:4)),
  # #                                                      paste0("iv",c(1:4),"_imp1"),
  # #                                                      paste0("iv",c(1:4),"_imp2"), 
  # #                                                      paste0("iv",c(1:4),"_imp1_lag1"))])
  # 
  # 
  # rm(ts)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                    parcel_size/1000,"km_",
                                    travel_time,"h_CA.rds")))
  
  rm(parcels)
}

