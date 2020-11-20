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
#           --> input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp
#           --> temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta
# 
#           baseline forest extents variables in cross-section (prepare in prepare_2000_forest_extents.R)
#           --> pattern: temp_data/processed_parcels/baseline_fc_cs_", for each parcel_size and IBS catchment_radius combination 
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

# This takes ~1h 
#### ADD N REACHABLE UML AND SAMPLE COVERAGE ####

# The sample coverage has to be computed as the number of reachable IBS-UML matched sample
# related to the number of reachable UML mills, with the all the former being included in the latter. 
# We cannot compute a ratio of all mills in the analysis sample - i.e. also including IBS not matched with 
# UML but with a desa centroid - because we would not know whether a mill from the latter group 
# is an additional mill that is not in UML or if it is in UML but it was not matched. 

# read the sample panel of IBS geolocalized mills
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# keep only those that are matched with UML. 
ibsuml <- ibs[ibs$uml_matched_sample==1,]
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

make_n_reachable_uml <- function(parcel_size, catchment_radius){
  
  # read the parcel panel
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  

  
  # make a spatial cross section of it (parcels' coordinates are constant over time)
  parcels_centro <- parcels[parcels$year == 1998, c("lonlat", "idncrs_lat", "idncrs_lon")]
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
  
  return(parcels)
}

catchment_radius <- 10000
while(catchment_radius < 60000){
  
  make_n_reachable_uml(parcel_size, catchment_radius) %>% 
    saveRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_reachable_uml_",
                             parcel_size/1000,"km_",
                             catchment_radius/1000,"CR.rds")))
  
  catchment_radius <- catchment_radius + 20000
}
  

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


catchment_radiuseS <- c(1e4, 3e4, 5e4)#
for(catchment_radius in catchment_radiuseS){
  # this is prepared in this script's previous part (add n reachable uml... )
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_reachable_uml_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  
  
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
 
  
  ### NEIGHBOR VARIABLE
  # create a grouping variable at the cross section (9 is to recall the the group id includes the 8 neighbors + the central grid cell.)
  parcels_cs <- parcels[!duplicated(parcels$lonlat),c("lonlat", "year", "idncrs_lat", "idncrs_lon")]
  
  # spatial
  parcels_cs <- st_as_sf(parcels_cs, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  
  # identify neighbors
  # this definition of neighbors includes the 8 closest, surounding, grid cells. 
  parcels_buf <- st_buffer(parcels_cs, dist = parcel_size - 10)
  row.names(parcels_buf) <- parcels_buf$lonlat
  sgbp <- st_intersects(parcels_buf)
  
  
  
  
  
  
  # SPATIAL TRENDS VARIABLES
  parcels$island_year <- paste0(parcels$island,"_",parcels$year)
  parcels$province_year <- paste0(parcels$province,"_",parcels$year)
  parcels$district_year <- paste0(parcels$district,"_",parcels$year)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                     parcel_size/1000,"km_",
                                     catchment_radius/1000,"CR.rds")))
}

  # parcels_list[[match(catchment_radius, catchment_radiuseS)]] <- parcels




#### TIME DYNAMICS VARIABLES ####
catchment_radiuseS <- c(1e4, 3e4, 5e4)# 
for(catchment_radius in catchment_radiuseS){ 
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_geovars_",
                                                  parcel_size/1000,"km_",
                                                  catchment_radius/1000,"CR.rds")))
  
  
  ### EXPORT SHARES FROM PCT TO FRACTION   
  parcels$wa_prex_cpo_imp1 <- parcels$wa_prex_cpo_imp1/100 
  parcels$wa_prex_cpo_imp2 <- parcels$wa_prex_cpo_imp2/100 

  
  ### LAGS AND LEADS OF A LARGE SET OF VARIABLES 
  
  variables <- c("wa_ffb_price_imp1", "wa_ffb_price_imp2", 
                "wa_cpo_price_imp1", "wa_cpo_price_imp2", "wa_prex_cpo_imp1","wa_prex_cpo_imp2",       
                #"wa_pko_price_imp1",       "wa_pko_price_imp2",       "wa_prex_pko_imp1",        "wa_prex_pko_imp2",       
                "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp",  "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",     
                #"wa_concentration_10",     "wa_concentration_30", "wa_concentration_50",     
                "n_reachable_ibs", "n_reachable_uml", "n_reachable_ibsuml",   "sample_coverage")
  
  for(voi in variables){
    ## lags
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
  }
  
  #parcels1 <- parcels
  
  
  ### Operations relating contemporaneous to past information - on prices only
  variables <- c("wa_ffb_price_imp1", "wa_ffb_price_imp2", 
                 "wa_cpo_price_imp1", "wa_cpo_price_imp2") #,"wa_pko_price_imp1", "wa_pko_price_imp2"
  
  for(voi in variables){
    
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
      
      
      ## and absolute deviation - SHORT RUN MEASURE -
      parcels <- mutate(parcels,
                        !!as.symbol(paste0(voi,"_dev_",py,"pya")) := !!as.symbol(paste0(voi)) - 
                                                                     !!as.symbol(paste0(voi,"_",py,"pya")))
      # # and relative deviation
      # parcels <- mutate(parcels,
      #                   !!as.symbol(paste0(voi,"_rdev_",py,"pya")) := (!!as.symbol(paste0(voi)) - 
      #                                                               !!as.symbol(paste0(voi,"_",py,"pya"))) /
      #  
      
      # Lag these deviations by one year
      parcels <- dplyr::arrange(parcels, lonlat, year)
      parcels <- DataCombine::slide(parcels,
                                  Var = paste0(voi,"_dev_",py,"pya"), 
                                  TimeVar = "year",
                                  GroupVar = "lonlat",
                                  NewVar = paste0(voi,"_dev_",py,"pya_lag1"),
                                  slideBy = -1, 
                                  keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, lonlat, year)
      
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
  }  
  
  # remove some variables that were only temporarily necessary
  vars_torm <- names(parcels)[grepl(pattern = "price_", x = names(parcels)) &
                              (grepl(pattern = "_lag2", x = names(parcels)) |
                              grepl(pattern = "_lag3", x = names(parcels)) |
                              grepl(pattern = "_lag4", x = names(parcels)) |
                              grepl(pattern = "_lag5", x = names(parcels)))]
  
  parcels <- parcels[,!(names(parcels) %in% vars_torm)]
  
  
  ### PRICE LOGARITHMS
  # select prices for which it's relevant/useful to compute the log
  price_variables <- names(parcels)[grepl(pattern = "price_", x = names(parcels)) &
                                      !grepl(pattern = "pko", x = names(parcels)) &
                                      !grepl(pattern = "dev", x = names(parcels)) &
                                      !grepl(pattern = "yoyg", x = names(parcels))]
  
  
  for(var in price_variables){
    parcels[,paste0("ln_",var)] <- log(parcels[,var])
  }
  
  
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))  
}





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


catchment_radiuseS <- c(1e4, 3e4, 5e4)# 
for(catchment_radius in catchment_radiuseS){ 
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))

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
  
  parcels <- merge(parcels, parcels_cs[,c("lonlat", "concession")], by = "lonlat")
  
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
  
  parcels <- merge(parcels, parcels_cs[,c("lonlat", "llu")], by = "lonlat")
  
  unique(parcels$llu)
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

  # yields many missing in illegal because many grid cells are within a mising land use legal classification
  # parcels[!duplicated(parcels$lonlat) & !is.na(parcels$llu), c("lonlat", "concession", "llu", "illegal1", "illegal2")]
  
  
  
  ### TIME SERIES & IV
  parcels <- merge(parcels, ts, by = "year")
  
  # Make the SHIFT SHARE INSTRUMENTAL VARIABLES 
  for(IMP in c(1,2)){
    for(SP in c(1:6)){
      parcels[,paste0("iv",SP,"_imp",IMP)] <- parcels[,paste0("wa_prex_cpo_imp",IMP,"_lag1")]*parcels[,paste0("spread",SP)] 
    }
  }
  
  # lag the iv variables
  ivS <- c(paste0("iv",c(1:6),"_imp1"), paste0("iv",c(1:6),"_imp2"))
  
  for(IV in ivS){
    parcels <- dplyr::arrange(parcels, lonlat, year)
    parcels <- DataCombine::slide(parcels,
                                  Var = IV, 
                                  TimeVar = "year",
                                  GroupVar = "lonlat",
                                  NewVar = paste0(IV,"_lag1"),
                                  slideBy = -1, 
                                  keepInvalid = TRUE)  
    parcels <- dplyr::arrange(parcels, lonlat, year)
  }
  
  # View(parcels[!is.na(parcels$wa_prex_cpo_imp1_lag1) &
  #                parcels$year>2007 &
  #                parcels$wa_prex_cpo_imp1_lag1!=0 ,c("lonlat" ,"year", paste0("wa_prex_cpo_imp",c(1,2),"_lag1"),
  #                                                      paste0("spread",c(1:4)),
  #                                                      paste0("iv",c(1:4),"_imp1"),
  #                                                      paste0("iv",c(1:4),"_imp2"), 
  #                                                      paste0("iv",c(1:4),"_imp1_lag1"))])
  
  
  
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"CR.rds")))
}

rm(cns, llu, parcels, parcels_cs, rspo, rspo_cs, sgbp, ts)


#### ADD N REACHABLE UML TO PARCELS IN ***UML*** CATCHMENT RADIUS #### 
# this is a different thing, it does not add variables to our sample for analysis, but to 
# another sample, that of parcels within CR of a **UML** mill, as outputed from prepare_lucpfip.R
# this is necessary to later compute the aggregation factor in demand for deforestation. 

# prepare geographic data once before looping
# island
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)
#province
province_sf <- st_read(file.path("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp"))
province_sf <- dplyr::select(province_sf, NAME_1)
province_sf_prj <- st_transform(province_sf, crs = indonesian_crs)
# district
district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
district_names <- read.dta13(file.path("temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta"))
district_names <- district_names[!duplicated(district_names$bps_),]
district_sf$d__2000 <- district_sf$d__2000 %>% as.character()
district_names$bps_ <- district_names$bps_ %>% as.character()
district_sf <- left_join(x = district_sf, y = district_names[,c("name_", "bps_")], 
                         by = c("d__2000" = "bps_"),
                         all = FALSE, all.x = FALSE, all.y = FALSE)


district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)

catchment_radius <- 10000
while(catchment_radius < 60000){
  # read the panel of parcels within CR of a **UML** mill, as outputed from prepare_lucpfip.R
  parcels <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"km_UML_CR.rds")))
  
  # make a spatial cross section of it (parcels' coordinates are constant over time)
  parcels_centro <- parcels[parcels$year == 2001, c("lonlat", "idncrs_lat", "idncrs_lon")]
  # (lon lat are already expressed in indonesian crs)
  parcels_centro <- st_as_sf(parcels_centro, coords = c("idncrs_lon", "idncrs_lat"), remove = T, crs = indonesian_crs)
  
  parcels$newv_uml <- rep(0, nrow(parcels))
  
  for(t in 1:length(years)){
    
    # UML
    # This is not a panel, so the information on presence or not a given year is whether 
    # the establishment year is anterior. We impute NA establishment year to be older than 1998. 
    present_uml <- uml[uml$est_year_imp <= years[t] | is.na(uml$est_year_imp),]
    
    annual_reachable_uml <- st_is_within_distance(parcels_centro, present_uml, dist = catchment_radius)
    parcels[parcels$year == years[t], "newv_uml"] <- lengths(annual_reachable_uml)
  } 
  colnames(parcels)[colnames(parcels) == "newv_uml"] <- paste0("n_reachable_uml")
  
  ### ADD GEOGRAPHIC VARIABLES
  parcels <- st_as_sf(parcels, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  
  # ISLAND variable
  parcels$island <- rep("", nrow(parcels))
  
  # make the operation faster by using island bbox (other wise the island polygons make 
  # the computation very long)
  # (and this also includes parcel centroids in the sea)
  island_sf_prj_bbox <- sapply(island_sf_prj$geometry, function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
  
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
  
  
  # PROVINCE variable
  
  # Work with a cross section for province and district attribution
  parcels_cs <- parcels[!duplicated(parcels$lonlat),]
  
  
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_prov_idx <- st_nearest_feature(parcels_cs, province_sf_prj)
  
  parcels_cs$province <- province_sf_prj$NAME_1[nearest_prov_idx]
  
  # DISTRICT variable
  # the nearest feature function enables to also grab those parcels which centroids are in the sea.
  nearest_dstr_idx <- st_nearest_feature(parcels_cs, district_sf_prj)
  
  # 4 parcels are closest to district with no name (NA) 
  parcels_cs$district <- district_sf_prj$name_[nearest_dstr_idx]
  
  parcels <- merge(st_drop_geometry(parcels),
                   st_drop_geometry(parcels_cs[,c("lonlat", "province", "district")]),
                   by = "lonlat")
  
  # REGIONAL TRENDS VARIABLES
  parcels$island_year <- paste0(parcels$island,"_",parcels$year)
  parcels$province_year <- paste0(parcels$province,"_",parcels$year)
  parcels$district_year <- paste0(parcels$district,"_",parcels$year)
  
  saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/lucpfip_panel_reachable_geovars_",
                                    parcel_size/1000,"km_",
                                    catchment_radius/1000,"km_UML_CR.rds")))
  
  
  catchment_radius <- catchment_radius + 20000
}

  
  
voi <- "wa_ffb_price_imp1"
View(parcels[500,c("lonlat", "year",
                                        voi,#
                                        paste0(voi,"_lag",c(1:4)),#
                                        paste0(voi,"_3ya"),
                                        paste0(voi,"_3ya_lag1"),
                                        paste0(voi,"_4ya"),
                                        paste0(voi,"_4ya_lag1"),
                                        paste0(voi,"_5ya"),
                                        paste0(voi,"_5ya_lag1"),
                                        paste0(voi,"_3pya"),#
                                        paste0(voi,"_3pya_lag1"),#
                                        paste0(voi,"_dev_3pya"),#
                                        paste0(voi,"_dev_3pya_lag1"),#
                                        paste0(voi,"_yoyg"),#
                                        paste0(voi,"_yoyg_lag", c(1:3)),#
                                        paste0(voi,"_yoyg_4ya"),
                                        paste0(voi,"_yoyg_4ya_lag1"),
                                        paste0(voi,"_yoyg_3pya"),#
                                        paste0(voi,"_yoyg_3pya_lag1"))])#


  





















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