# This script uses dataframes of parcels built in prepare_lucpfip.R 
# These dataframes have parcels of a given size (3x3km only for now), 
# within respectively 10km, 30km and 50km catchment radius of IBS mills
# For each combination of parcel size and catchment radius, there is a different parcel dataframe. 

# Now we would like to add explanatory variables to these dataframes.  

# So we repeat the following main steps for each parcel_size X catchment_radius combination. 

# 1. Take one set of parcels - their geographic characteristics are constant. 
# Have it as an sf obect with point data (coordinates of centroids)

# 2. import the panel of IBS geolocalized mills. 

# Add to the sf dataframe parcel object 18 list columns (for years 1998 to 2015). Each row has one list element in each of these list columns.
# Each of these elements is a dataframe with one record for each mill that is less than 40 km from the element's point. 

# Each of the records within each of the 18 annual elements of every parcels is an IBS mill, with its attributes from IBS_UML_panel_final
# + its distance to this particular parcel's centroid. 
# + attribute specific weights based on this distance and the distance of other reachable mills that have a non missing value for this attribute.   
# + its weighted value for each variable of interest. 

# Add to the sf dataframe 18*K columns (with K explanatory variables of interest). Each row has the sum of the column with weighted values for this variable in the annual element

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("data.table", "dplyr", "readstata13", "readxl",
                   "rgdal", "sf",
                   "doParallel", "foreach", "parallel")
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

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 



### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/processed_parcels/temp_cs_wa_explanatory")


### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"


### IBS YEARS
years <- seq(from = 1998, to = 2015, by = 1)

##### PREPARE IBS DATA ##### 

ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))  

# keep only geolocalized mills
ibs <- ibs[ibs$analysis_sample == 1,]
length(unique(ibs$firm_id))

# keep only the variables that identify mills and those which we want to distribute over parcels. 
ibs <- ibs[, c("firm_id", "year", "trase_code", "uml_id", "mill_name", "parent_co", "lat", "lon",
               "min_year","est_year", "est_year_imp", "max_year", 
               "ffb_price_imp1", #"ffb_price_imp2", "in_ton_ffb_imp1", "in_ton_ffb_imp2", "in_val_ffb_imp1", "in_val_ffb_imp2",
               "cpo_price_imp1"  #,"cpo_price_imp2", "out_ton_cpo_imp1", "out_ton_cpo_imp2", "out_val_cpo_imp1", "out_val_cpo_imp2",
               # "prex_cpo_imp1", "prex_cpo_imp2", "lag1_prex_cpo_imp1", "lag1_prex_cpo_imp2", "avg_prex_cpo_imp1", "avg_prex_cpo_imp2", 
               # "pko_price_imp1","pko_price_imp2", "out_ton_pko_imp1", "out_ton_pko_imp2", "out_val_pko_imp1", "out_val_pko_imp2",
               # "prex_pko_imp1", "prex_pko_imp2",
               # "export_pct_imp", "revenue_total", "workers_total_imp3",
               # "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp"
)]#,
# "iv2_imp1", "iv2_imp2", "iv3_imp1", "iv3_imp2", "iv4_imp1", "iv4_imp2", 
# "concentration_10", "concentration_30", "concentration_50"

# we don't keep the logs because we don't want to compute means of logs, but logs of means. 
# "ffb_price_imp1_ln", "ffb_price_imp2_ln", "cpo_price_imp1_ln", "cpo_price_imp2_ln",        
#           "pko_price_imp1_ln", "pko_price_imp2_ln", "out_val_cpo_imp1_ln", "out_val_cpo_imp2_ln", "out_val_pko_imp1_ln",      
#          "out_val_pko_imp2_ln", "revenue_total_ln" 


# split the panel into sf cross sections 
class(ibs$year)
ibs_cs <- lapply(years, FUN = function(x) ibs[ibs$year == x,]) 
ibs_cs <- lapply(ibs_cs, FUN = st_as_sf, coords =  c("lon", "lat"), remove = TRUE, crs = 4326)
ibs_cs <- lapply(ibs_cs, FUN = st_transform, crs = indonesian_crs)
# later we need these coordinates to turn nested ibs_cs to sf objects. 
ibs_cs <- lapply(ibs_cs, FUN = function(cs){mutate(cs, indo_crs_lon = st_coordinates(cs)[,"X"])})
ibs_cs <- lapply(ibs_cs, FUN = function(cs){mutate(cs, indo_crs_lat = st_coordinates(cs)[,"Y"])})


##### COMPUTE WEIGHTED AVERAGES OF MILL VARIABLES FOR PARCEL OF A GIVEN SET OF PARCELS ##### 
# 
parcel_size <- 5000
catchment_radius <- 30000
t <- 1

# note on addressing geometries: the two first of the following calls are equivalent; however, the third is different 
# parcels$geometry[parcels$lonlat == i] 
# parcels[["geometry"]][parcels$lonlat == i]
# parcels[parcels$lonlat == i, "geometry"]

parcel_set_w_average <- function(parcel_size, catchment_radius){
  
  #### Prepare parcel panel ####
  # Import the parcel panel (for IBS)
  parcels_centro <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_dynamics_",
                                             parcel_size/1000,"km_",
                                             catchment_radius/1000,"km_IBS_CR.rds")))
  # keep only one cross-section, no matter which. 
  parcels_centro <- dplyr::filter(parcels_centro, year == 2001)
  # turn it into a sf object (idncrs_lon idncrs_lat are already expressed in indonesian crs)
  parcels_centro <- st_as_sf(parcels_centro, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  parcels_centro <- dplyr::select(parcels_centro, lonlat, geometry)
  
  
  
  
  #### Compute panel of weighted averages ####
  
  ### Function description 
  # 
  # 
  # 
  parallel_w_averages <- function(ncores){
    
    ## Sequence over which to execute the task 
    years <- seq(from = 1998, to = 2015, by = 1)
    
    ## read the input to the task
    # has been done beforehand because it is not specific to function variables
    
    ## Define the task 
    # (Computes cross-sectional weighted averages of mill variables at parcels)
    annual_w_averages <- function(t){
      #let's not be year specific in this function, and we will rename and append everything after. 
      
      ## Attribute to each parcel centroid the sf data frame of reachable mills 
      
      # this is a data frame of pairs of parcel and year t mill points that are within catchment_radius
      # ***the geometry kept is from x ***
      d <- st_join(x = parcels_centro, 
                   y = ibs_cs[[t]], 
                   join = st_is_within_distance, 
                   dist = catchment_radius, 
                   left = FALSE) # to return only records that spatially match (i.e. an inner join and not a left join)
      
      # nest the sets of reachable mills within each parcel row.
      # they need to be no sf object for that. 
      d <- st_set_geometry(d, NULL)
      parcels <- nest_join(parcels_centro, d, 
                           by = "lonlat", 
                           keep = T, # keep = T keeps lonlat in the nested df.
                           name = "reachable") %>% st_as_sf() # (bc the nest_join removes the sf class)
      rm(d)
      
      # select non empty reachable nested data frames (data frames of reachable mills) - programing purpose
      s <- sapply(parcels$reachable, FUN = nrow)>0
      
      # compute the number of reachable mills at each parcel - informative purpose
      parcels[,"n_reachable_ibs"] <- sapply(parcels$reachable, 
                                            FUN = nrow)
      # turn (non-empty) nested tibbles into sf data frames. 
      parcels$reachable[s] <- lapply(parcels$reachable[s], 
                                     FUN = st_as_sf, 
                                     coords =  c("indo_crs_lon", "indo_crs_lat"), 
                                     remove = FALSE, 
                                     crs = indonesian_crs)
      # # at this stage let's store our working object somewhere 
      # parcels_saved <- parcels
      
      # make a distance column in the reachable data frames 
      parcels$reachable[s] <- lapply(parcels$lonlat[s], 
                                     FUN = function(i){
                                       mutate(parcels$reachable[parcels$lonlat == i][[1]], 
                                              distance = st_distance(x = parcels$geometry[parcels$lonlat == i], # parcel geom
                                                                     y = parcels$reachable[parcels$lonlat == i][[1]]$geometry, # reachable mill geom
                                                                     by_element = FALSE) %>% as.numeric()
                                       )
                                     })
      # #  parcels$lonlat[s][52:length(parcels$lonlat[s])]
      # for(i in parcels$lonlat[s]){
      #   mutate(parcels$reachable[parcels$lonlat == i][[1]], 
      #          distance = st_distance(x = parcels$geometry[parcels$lonlat == i], 
      #                                 y = parcels$reachable[parcels$lonlat == i][[1]]$geometry,
      #                                 by_element = FALSE) %>% as.numeric()
      #   )
      # }
      # match(i, parcels$lonlat[s])
      
      # make the inverse of distance column
      parcels$reachable[s] <- lapply(parcels$lonlat[s], 
                                     FUN = function(i){
                                       mutate(parcels$reachable[parcels$lonlat == i][[1]], 
                                              w = 1/distance
                                       )
                                     })
      
      # Define the variables of interest we want to compute the weighted averages of. 
      # variables <- c("ffb_price_imp1", "ffb_price_imp2")      
      variables <- c("est_year_imp", # "min_year", "est_year", "max_year",
                     "ffb_price_imp1", # "ffb_price_imp2", # "in_ton_ffb_imp1", "in_ton_ffb_imp2", "in_val_ffb_imp1", "in_val_ffb_imp2",
                     "cpo_price_imp1"  #,"cpo_price_imp2", # "out_ton_cpo_imp1", "out_ton_cpo_imp2", "out_val_cpo_imp1", "out_val_cpo_imp2",
                     # "prex_cpo_imp1", "prex_cpo_imp2", "lag1_prex_cpo_imp1", "lag1_prex_cpo_imp2", "avg_prex_cpo_imp1", "avg_prex_cpo_imp2", 
                     #"pko_price_imp1","pko_price_imp2", # "out_ton_pko_imp1", "out_ton_pko_imp2", "out_val_pko_imp1", "out_val_pko_imp2",
                     #"prex_pko_imp1", "prex_pko_imp2",
                     # "export_pct_imp", "revenue_total", "workers_total_imp3",
                     #"pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp"
      )
      
      
      # make the variable specific sum of the inverse of distance over all the reachable mills that have no missing on this variable.
      for(voi in variables){
        parcels$reachable[s] <- lapply(parcels$lonlat[s], 
                                       FUN =function(i){
                                         voi_missing <- parcels$reachable[parcels$lonlat == i][[1]][,voi] %>% st_drop_geometry() %>% is.na() %>% as.vector() 
                                         
                                         mutate(parcels$reachable[parcels$lonlat == i][[1]],
                                                !!as.symbol(paste0("sum_w_",voi)) := ifelse(voi_missing, 
                                                                                            yes = NA, 
                                                                                            no = sum(w[!voi_missing])))
                                       })
        
        # for(i in parcels$lonlat[s]){
        #   voi_missing <- parcels$reachable[parcels$lonlat == i][[1]][,voi] %>% st_drop_geometry() %>% is.na() %>% as.vector() 
        #   parcels$reachable[parcels$lonlat == i][[1]][voi_missing, paste0("sum_w_",voi)] <- NA
        # } # this solution takes 4.85s. rather than 7.3 with "full looping"
        # # the whole time is taken by the for loop, but I cannot think of a way to do this 
        # # replace with NA within the lapply  
        
        
        parcels$reachable[s] <- lapply(parcels$lonlat[s], 
                                       FUN =function(i){
                                         mutate(parcels$reachable[parcels$lonlat == i][[1]],
                                                # make the standardized weights. They are variable specific too.
                                                # this ratio indeed is NA if sum_w_voi is NA
                                                !!as.symbol(paste0("std_w_", voi)) := w/!!as.symbol(paste0("sum_w_", voi)),
                                                # compute an intermediate column in reachable in which each mill gets the product of its weight and its attribute
                                                !!as.symbol(paste0("w_var_", voi)) := !!as.symbol(paste0("std_w_", voi))*(!!as.symbol(voi))
                                         )
                                       })
        
        # build the column in parcels in which every cell is the weighted average of voi
        # sum these weighted terms and hence compute the weighted averages
        # it makes sure that parcels whose reachable mills are all NA on a voi don't get a 0 but a NA for weighted mean)
        parcels$var_template <- rep(NA, nrow(parcels)) 
        colnames(parcels)[colnames(parcels) == "var_template"] <- paste0("wa_", voi)
        parcels[s,paste0("wa_", voi)] <- sapply(parcels$reachable[s], 
                                                FUN = function(x) x[,paste0("w_var_",voi)] 
                                                %>% st_drop_geometry() 
                                                %>% is.na()
                                                %>% all() 
                                                %>% ifelse(yes = NA, no = sum(st_drop_geometry(x[,paste0("w_var_",voi)]),na.rm = T)))
        
      }# closes the loop on variables
      
      parcels$idncrs_lon <- st_coordinates(parcels)[,"X"]
      parcels$idncrs_lat <- st_coordinates(parcels)[,"Y"]
      parcels <- st_drop_geometry(parcels)
      
      parcels <- mutate(parcels, 
                        reachable = NULL)
      
      # give year specific variable names to the variables built in this function 
      names(parcels) <- names(parcels) %>% paste0(".", years[t])
      
      saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/temp_cs_wa_explanatory/cs_wa_explanatory_",
                                        parcel_size/1000,"km_",catchment_radius/1000,"CR_",years[t],".rds")))
    }# closes annual_w_averages
    
    
    
    ## register cluster
    registerDoParallel(cores = ncores) 
    
    ## define foreach object 
    foreach(t = 1:length(years), 
            #.combine = cbind,
            .multicombine = TRUE, # not necessary if .combine = cbind because then the default multicombine is TRUE anyways
            .inorder = FALSE, # we don't care that the results be combined in the same order they were submitted
            .export = c("parcels_centro", "ibs_cs", "parcel_size", "catchment_radius", "indonesian_crs"), 
            .packages = c("sf", "dplyr")) %dopar% 
      annual_w_averages(t) # the function that is parallelly applied to different years. 
  }# closes parallel_w_averages
  
  ### Execute it
  parallel_w_averages(min(detectCores() - 1, 4))
  
  
  ### Read annual cross-sections and bind them 
  annual_parcel_paths <- list.files(path = file.path("temp_data/processed_parcels/temp_cs_wa_explanatory/"), 
                                    pattern = paste0("cs_wa_explanatory_",parcel_size/1000,"km_",catchment_radius/1000,"CR_"), 
                                    full.names = TRUE) 
  
  wide_parcels <- lapply(annual_parcel_paths, 
                         FUN = function(x){readRDS(x)}) %>% bind_cols()
  
  # manage repetitions of lonlat variables over years
  wide_parcels$lonlat <- wide_parcels$lonlat.1998
  wide_parcels <- dplyr::select(wide_parcels, lonlat, everything()) # just changes columns order
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("lonlat."))  
  
  # manage repetitions of lat and lon variables over years - NO BECAUSE THEY ARE NOT THERE ANYMORE
  # wide_parcels$lat <- wide_parcels$lat.1998
  # wide_parcels <- dplyr::select(wide_parcels, -starts_with("lat."))
  # wide_parcels$lon <- wide_parcels$lon.1998
  # wide_parcels <- dplyr::select(wide_parcels, -starts_with("lon."))
  
  # manage repetitions of idncrs_lat and idncrs_lon variables over years
  wide_parcels$idncrs_lat <- wide_parcels$idncrs_lat.1998
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("idncrs_lat."))
  wide_parcels$idncrs_lon <- wide_parcels$idncrs_lon.1998
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("idncrs_lon."))
  
  # saveRDS(wide_parcels, file = here(paste0("build/output/wa_wide_panel_parcels_",parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))
  
  # reshape to long 
  varying_vars <- wide_parcels %>% dplyr::select(-lonlat, -idncrs_lat, -idncrs_lon) %>% colnames()   # -lat, -lon, 
  
  long_parcels <- reshape(wide_parcels, 
                          varying = varying_vars, 
                          #v.names =  
                          timevar = "year",
                          idvar = "lonlat",
                          ids = "lonlat",
                          direction = "long",
                          sep = ".")
  
  
  rm(varying_vars)
  
  long_parcels <- dplyr::arrange(long_parcels, lonlat, year)
  
  saveRDS(long_parcels, file = file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))
  
}

##### EXECUTE FUNCTION #####
# parcel_size = 1000 is all the point of running this script. 
CR <- 30
while(CR < 60){
  parcel_set_w_average(parcel_size = 5000, catchment_radius = CR*1000)
  CR <- CR + 20                   
}




