# This script uses dataframes of parcels built in make_osrm_CA.R 
# These dataframes have parcels of a given size (3x3km only for now), 
# within respectively 2h, 4h and 6h driving catchment areas of IBS mills
# For each combination of parcel size and catchment area, there is a different parcel data frame. 

# IN THIS SCRIPT, CONTRARY TO wa_at_parcels_distances.R EVERYTHING IS DONE WITHIN ISLANDS
# this is because the driving travel times need to be computed first for every possible parcel-mill pair. Therefore,  
# Now we would like to add explanatory variables to these dataframes.  

# So we repeat the following main steps for each parcel_size X catchment area combination. 

# 1. Take one set of parcels - their geographic characteristics are constant. 
# Have it as an sf obect with point data (coordinates of centroids)

# 2. import the panel of IBS geolocalized mills. 

# Add to the sf dataframe parcel object 18 list columns (for years 1998 to 2015). Each row has one list element in each of these list columns.
# Each of these elements is a dataframe with one record for each mill that is less than X hours driving from the element's point. 

# Each of the records within each of the 18 annual elements of every parcels is an IBS mill, with its attributes from IBS_UML_panel_final
# + its travel time to this particular parcel's centroid. 
# + attribute specific weights based on this time and the time to other reachable mills that have a non missing value for this attribute.   
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
neededPackages = c("data.table", "tidyr", "dplyr", "readstata13", "readxl",
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
dir.create("temp_data/processed_parcels/temp_cs_wa_explanatory2")


### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE 
parcel_size <- 3000


### IBS YEARS
years <- seq(from = 1998, to = 2015, by = 1)




##### COMPUTE WEIGHTED AVERAGES OF MILL VARIABLES FOR PARCEL OF A GIVEN SET OF PARCELS ##### 
# 
# 
# catchment_radius <- 30000
# t <- 1

# note on addressing geometries: the two first of the following calls are equivalent; however, the third is different 
# parcels$geometry[parcels$parcel_id == i] 
# parcels[["geometry"]][parcels$parcel_id == i]
# parcels[parcels$parcel_id == i, "geometry"]
# 
island <- "Sumatra"
parcel_size <- 3000
travel_time <- 2
t <- 7


parcel_set_w_average <- function(island, parcel_size, travel_time){
  ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))  
  
  # keep only geolocalized mills on the island of interest
  ibs <- ibs[ibs$analysis_sample == 1,]
  ibs <- ibs[ibs$island_name == island,]
  
  # keep only the variables that identify mills and those which we want to distribute over parcels. 
  ibs <- ibs[, c("firm_id", "year", "lat", "lon",#"trase_code", "uml_id", "mill_name", "parent_co", 
                "est_year_imp", "min_year","est_year",  "max_year",#  
                 "ffb_price_imp1", "ffb_price_imp2", #"in_ton_ffb_imp1", "in_ton_ffb_imp2", "in_val_ffb_imp1", "in_val_ffb_imp2",
                 "cpo_price_imp1","cpo_price_imp2", #"out_ton_cpo_imp1", "out_ton_cpo_imp2", "out_val_cpo_imp1", "out_val_cpo_imp2",
                 "prex_cpo_imp1", "prex_cpo_imp2",
                 "pko_price_imp1","pko_price_imp2", #"out_ton_pko_imp1", "out_ton_pko_imp2", "out_val_pko_imp1", "out_val_pko_imp2",
                 "prex_pko_imp1", "prex_pko_imp2",
                 "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp")]
  
  # we don't keep the logs because we don't want to compute means of logs, but logs of means. 
  # "ffb_price_imp1_ln", "ffb_price_imp2_ln", "cpo_price_imp1_ln", "cpo_price_imp2_ln",        
  #           "pko_price_imp1_ln", "pko_price_imp2_ln", "out_val_cpo_imp1_ln", "out_val_cpo_imp2_ln", "out_val_pko_imp1_ln",      
  #          "out_val_pko_imp2_ln", "revenue_total_ln" 
  
  
  # split the panel into sf cross sections 
  ibs_cs <- lapply(years, FUN = function(x) ibs[ibs$year == x,]) 

  # ibs_cs <- lapply(ibs_cs, FUN = st_as_sf, coords =  c("lon", "lat"), remove = TRUE, crs = 4326)
  # ibs_cs <- lapply(ibs_cs, FUN = st_transform, crs = indonesian_crs)
  # later we need these coordinates to turn nested ibs_cs to sf objects. 
  # ibs_cs <- lapply(ibs_cs, FUN = function(cs){mutate(cs, indo_crs_lon = st_coordinates(cs)[,"X"])})
  # ibs_cs <- lapply(ibs_cs, FUN = function(cs){mutate(cs, indo_crs_lat = st_coordinates(cs)[,"Y"])})


  #### Prepare parcel panel ####
  # Import the parcel panel (for IBS)
  parcels_centro <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_CA_total.rds")))
  # keep only one cross-section, no matter which as this is a balanced panel
  parcels_centro <- parcels_centro[!duplicated(parcels_centro$parcel_id),]
  parcels_centro <- mutate(parcels_centro, lonlat = paste0(lon, lat))
  parcels_centro <- dplyr::select(parcels_centro, -lucpfip_pixelcount_total, -lucpfip_ha_total)
  
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

      dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_",years[t])))
      dur_mat <- dur_mat$durations

      if(nrow(ibs_cs[[t]]) != ncol(dur_mat)){stop(paste0("duration matrix and mill cross section don't match in ",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_",years[t]))}


      # This is the driving travel time (duration) matrix between each pair of parcel and mill.
      # no matter the t, this matrix has the same amount of rows (parcels) as parcels_centro

      # for more safety, merge the dur_mat with the parcels centro based on coordinates identifiers
      dur_mat <- as.data.frame(dur_mat)
      dur_mat$lonlat <- row.names(dur_mat)
      # sort = FALSE is MEGA IMPORTANT because otherwise dur_mat get sorted in a different way than
      dur_mat <- merge(dur_mat, parcels_centro[,c("lonlat", "parcel_id")], by = "lonlat", all = FALSE, sort = FALSE)
      row.names(dur_mat) <- dur_mat$lonlat
      # so there should be the same amount of parcels which coordinates matched, as the total amount of parcels.
      if(nrow(dur_mat) != nrow(parcels_centro)){stop(paste0("duration matrix and parcel cross section did not merge in ",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_",years[t]))}

      # additional check:
      parcels <- parcels_centro
      row.names(parcels) <- parcels$lonlat

      parcels <- dplyr::arrange(parcels,parcel_id)
      dur_mat <- dplyr::arrange(dur_mat,parcel_id)

      if(!(all.equal(row.names(parcels), row.names(dur_mat)))){stop()}

      # parcels <- parcels[1:2500,]
      # dur_mat <- dur_mat[1:2500,]

      # nest the sets of reachable mills within each parcel row.
      dur_mat <- dplyr::select(dur_mat, -lonlat, -parcel_id)
      dur_mat <- as.matrix(dur_mat)

      # convert durations (in minutes) into logical whether each parcel-mill pair's travel time is inferior to a threshold travel_time
      dur_mat_log <- dur_mat/(60) < travel_time

      anyNA(dur_mat_log)
      dur_mat_log <- replace_na(dur_mat_log, replace = FALSE)
      anyNA(dur_mat_log)

      list_col <- list()
      length(list_col) <- nrow(parcels)
      parcels$reachable <- list_col
      rm(list_col)


      # nest reachable mills within each parcel record.
      parcels$reachable <- lapply(1:nrow(parcels), FUN = function(i){ibs_cs[[t]][dur_mat_log[i,],]})

      # select non empty reachable nested data frames (data frames of reachable mills) - programming purpose
      s <- sapply(parcels$reachable, FUN = nrow)>0

      # compute the number of reachable mills at each parcel - informative purpose
      parcels[,"n_reachable_ibs"] <- sapply(parcels$reachable, FUN = nrow)


      # for(i in 1:nrow(parcels[s,])){
      #   # add to this reachable mill data set a variable for their respective durations to parcel i.
      #   parcels[s,][i,]$reachable[[1]]$durations <- dur_mat[s,][i,dur_mat_log[s,][i,]] #%>% as.vector()
      # }

      parcels$reachable[s] <- lapply(row.names(parcels)[s],
                                  FUN = function(ids){
                                    mutate(parcels[ids,]$reachable[[1]],
                                           durations = dur_mat[ids,dur_mat_log[ids,]]) %>% as.vector()
                                  })

      # compute the inverse of this duration
      parcels$reachable[s] <- lapply(parcels$reachable[s], mutate, w = 1/durations)

      
      parcels["98.00524087692564.42566807199589",]$reachable[[1]] %>% ncol()
      parcels[ids,]$reachable[[1]]
      ids <- "103.048821259831-1.59157153981874"
      dur_mat[ids,dur_mat_log[ids,]]

      
       for(ids in row.names(parcels)[s]){
        if(!all.equal(as.character(parcels[ids,]$reachable[[1]]$firm_id), names(dur_mat[ids,])[dur_mat_log[ids,]])){stop("PB")}
       }     
      

      
      
      # Define the variables of interest we want to compute the weighted averages of.
      variables <- "ffb_price_imp1"
      # variables <- c("est_year_imp", # "min_year", "est_year", "max_year",
      #                "ffb_price_imp1", "ffb_price_imp2", # "in_ton_ffb_imp1", "in_ton_ffb_imp2", "in_val_ffb_imp1", "in_val_ffb_imp2",
      #                "cpo_price_imp1", "cpo_price_imp2", # "out_ton_cpo_imp1", "out_ton_cpo_imp2", "out_val_cpo_imp1", "out_val_cpo_imp2",
      #                "prex_cpo_imp1", "prex_cpo_imp2",
      #                "pko_price_imp1","pko_price_imp2", # "out_ton_pko_imp1", "out_ton_pko_imp2", "out_val_pko_imp1", "out_val_pko_imp2",
      #                "prex_pko_imp1", "prex_pko_imp2",
      #                # "export_pct_imp", "revenue_total", "workers_total_imp3",
      #                "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp")
      #
      #
      # make the variable specific sum of the inverse of durations over all the reachable mills that have no missing on this variable.
      for(voi in variables){
        parcels$reachable[s] <- lapply(parcels$parcel_id[s],
                                       FUN =function(i){
                                         voi_missing <- parcels$reachable[parcels$parcel_id == i][[1]][,voi] %>% is.na() %>% as.vector()

                                         mutate(parcels$reachable[parcels$parcel_id == i][[1]],
                                                !!as.symbol(paste0("sum_w_",voi)) := ifelse(voi_missing,
                                                                                            yes = NA,
                                                                                            no = sum(w[!voi_missing])))
                                       })

        # for(i in parcels$parcel_id[s]){
        #   voi_missing <- parcels$reachable[parcels$parcel_id == i][[1]][,voi] %>% st_drop_geometry() %>% is.na() %>% as.vector()
        #   parcels$reachable[parcels$parcel_id == i][[1]][voi_missing, paste0("sum_w_",voi)] <- NA
        # } # this solution takes 4.85s. rather than 7.3 with "full looping"
        # # the whole time is taken by the for loop, but I cannot think of a way to do this
        # # replace with NA within the lapply


        parcels$reachable[s] <- lapply(parcels$parcel_id[s],
                                       FUN =function(i){
                                         mutate(parcels$reachable[parcels$parcel_id == i][[1]],
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
                                                %>% is.na()
                                                %>% all()
                                                %>% ifelse(yes = NA, no = sum(x[,paste0("w_var_",voi)],na.rm = T)))

      }# closes the loop on variables

      parcels <- mutate(parcels,
                        reachable = NULL)

      # remove the lonlat ids
      parcels <- dplyr::select(parcels, -lonlat)
      row.names(parcels) <- NULL

      # give year specific variable names to the variables built in this function
      names(parcels) <- names(parcels) %>% paste0(".", years[t])



      saveRDS(parcels, file.path(paste0("temp_data/processed_parcels/temp_cs_wa_explanatory2/cs_wa_explanatory_",
                                        island,"_",parcel_size/1000,"km_",travel_time,"h_CA_",years[t],".rds")))
      
    }# closes annual_w_averages
    
    
    
    ## register cluster
    registerDoParallel(cores = ncores) 
    
    ## define foreach object 
    foreach(t = 1:4, #length(years)
            #.combine = cbind,
            .multicombine = TRUE, # not necessary if .combine = cbind because then the default multicombine is TRUE anyways
            .inorder = FALSE, # we don't care that the results be combined in the same order they were submitted
            .export = c("parcels_centro", "ibs_cs", "island", "parcel_size", "travel_time", "indonesian_crs"), 
            .packages = c("tidyr","dplyr")) %dopar% 
      annual_w_averages(t) # the function that is parallely applied to different years. 
  }# closes parallel_w_averages
  
  ### Execute it
  system.time(
  parallel_w_averages(detectCores() - 1)
  )
  
  
  ### Read annual cross-sections and bind them 
  annual_parcel_paths <- list.files(path = file.path("temp_data/processed_parcels/temp_cs_wa_explanatory2/"), 
                                    pattern = paste0("cs_wa_explanatory_",island,"_",parcel_size/1000,"km_",travel_time,"h_CA_"), 
                                    full.names = TRUE) 
  
  wide_parcels <- lapply(annual_parcel_paths, 
                         FUN = function(x){readRDS(x)}) %>% bind_cols()
  
  # manage repetitions of parcel_id variables over years
  wide_parcels$parcel_id <- wide_parcels$parcel_id.1998
  wide_parcels <- dplyr::select(wide_parcels, parcel_id, everything())
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("parcel_id."))  
  # manage repetitions of lat and lon variables over years
  wide_parcels$lat <- wide_parcels$lat.1998
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("lat."))
  wide_parcels$lon <- wide_parcels$lon.1998
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("lon."))
  # manage repetitions of idncrs_lat and idncrs_lon variables over years
  wide_parcels$idncrs_lat <- wide_parcels$idncrs_lat.1998
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("idncrs_lat."))
  wide_parcels$idncrs_lon <- wide_parcels$idncrs_lon.1998
  wide_parcels <- dplyr::select(wide_parcels, -starts_with("idncrs_lon."))
  
  # saveRDS(wide_parcels, file = here(paste0("build/output/wa_wide_panel_parcels_",parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))
  
  # reshape to long 
  varying_vars <- wide_parcels %>% dplyr::select(-parcel_id, -lat, -lon, - idncrs_lat, -idncrs_lon) %>% colnames()
  
  long_parcels <- reshape(wide_parcels, 
                          varying = varying_vars, 
                          #v.names =  
                          timevar = "year",
                          idvar = "parcel_id",
                          ids = "parcel_id",
                          direction = "long",
                          sep = ".")
  
  
  rm(varying_vars)
  
  long_parcels <- dplyr::arrange(long_parcels, parcel_id, year)
  
  saveRDS(long_parcels, file = file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                                island,"_",parcel_size/1000,"km_",travel_time,"h_CA.rds")))
  
}

##### EXECUTE FUNCTION #####

for(ISL in c("Sumatra", "Kalimantan")){
  for(TT in c(2,4,6)){
    parcel_set_w_average(island = ISL, parcel_size = 3000, travel_time = TT)
  }
}


### Gather the lucfip variables for each parcel_size and catchment area combinations. ####
PS <- 3000
# for each catchment area (2, 4, and 6 hours) stack the RHS parcel data from the two islands.  
for(travel_time in c(2,4,6)){
  # For each Island, join columns of 
  pf_df_list <- list()
  islandS <- c("Sumatra", "Kalimantan")
  for(island in islandS){

        df_total    <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                                island,"_",parcel_size/1000,"km_",travel_time,"h_CA.rds")))
        
        # leave the island information in a variable 
        df_total$island <- island
        
        pf_df_list[[match(island, islandS)]] <- df_total

  }

  # stack the islands together
  indo_df <- bind_rows(pf_df_list)

  indo_df <- dplyr::select(indo_df, parcel_id, year,
                            everything())

  saveRDS(indo_df, file = file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                                           parcel_size/1000,"km_",travel_time,"h_CA.rds")))

  rm(indo_df, pf_df_list)

}