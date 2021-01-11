### Merge panel data frames of parcels with LHS and RHS variables 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("plyr", "dplyr", "foreign", 
                   "DataCombine")

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
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# There are several merges to do because several data frames
# Each data frame is a different set of parcels, depending on parcel size and catchment radius. 

# For each set of parcels, there are also different sets of outcome variables 
# (lucfip, lucpfip, emissions, lucpfsmp...)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# parcel_size <- 3000
# catchment_radius <- 3e4

merge_lhs_rhs <- function(parcel_size, catchment_radius){
  
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

  ## make a variable that counts rapid and slow lucfp events (this is only computed for primary forest as of now)
  # THIS IS GOING TO BE THE MAIN OUTCOME VARIABLE INSTEAD OF lucpfip_pixelcount_total
  LHS$lucpfip_pixelcount <- LHS$lucpfip_rapid_pixelcount + LHS$lucpfip_slow_pixelcount
  LHS$lucfip_pixelcount <- LHS$lucfip_pixelcount_30th# pour l'instant on met lucfip_pixelcount_total dans "all producers" et pas rapid + slow, car on n'a 
  # pas calculé rapid et slow pour ce type de forêt encore
  
  # make variable that counts lucfp events on both small and medium sized plantations 
  LHS$lucpfsmp_pixelcount <- LHS$lucpfsp_pixelcount_total + LHS$lucpfmp_pixelcount_total
  LHS$lucfsmp_pixelcount <- LHS$lucfsp_pixelcount_30th + LHS$lucfmp_pixelcount_30th
  
  ## make a variable that counts lucfp events on all types of plantations
  LHS$lucpfap_pixelcount <- LHS$lucpfip_pixelcount + LHS$lucpfsmp_pixelcount
  LHS$lucfap_pixelcount <- LHS$lucfip_pixelcount + LHS$lucfsmp_pixelcount 


  
  
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
  
  
  ### EXPLICATIVE VARIABLES (runs from 1998-2015)
  RHS <-  readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                    parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))
  
  # # this is temporary necessary, but in the new version of wa_at_parcels_distances, lonlat is already available in RHS and merging can be directly executed with  #  
  # library(sf)
  # indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
  # RHS <- st_as_sf(RHS, coords = c("lon", "lat"), crs = indonesian_crs, remove = FALSE)
  # RHS <- st_transform(RHS, crs = 4326)
  # RHS$lon <- st_coordinates(RHS)[,"X"]%>% round(6) # the rounding is bc otherwise there are very little differences in the decimals of the coordinates... 
  # RHS$lat <- st_coordinates(RHS)[,"Y"]%>% round(6) 
  # RHS <- mutate(RHS, lonlat = paste0(lon, lat))
  # RHS <- st_drop_geometry(RHS)
  # RHS <- dplyr::select(RHS, -lat, -lon, -parcel_id)
  
  # Final code if wa_at_CR_parcels.R and add_CR_parcels.R are run again
  RHS <- dplyr::select(RHS, -lat, -lon, -idncrs_lat, -idncrs_lon)

  # MERGE
  # years 1998 - 2000 from RHS will not match, we don't need to keep them because the information 
  # from these years is captured in add_parcel_variables.R within lag variables. Hence the inner_join
  # one can check that without Papua nor 1998-2000 obs., RHS has the same number of rows 
  # RHS <- RHS[RHS$island != "Papua",]
  # RHS <- RHS[RHS$year > 2000,]  
  parcels <- inner_join(LHS, RHS, by = c("lonlat", "year")) #  
  

  if(nrow(parcels) != nrow(LHS)){stop("LHS and RHS don't have the exact same set of parcels")}
  rm(LHS, RHS)  

  
  
  
  ### OUTCOME VARIABLE TIME DYNAMICS
  
  # outcome_variables <- c("lucpfip_ha_intact", "lucpfip_ha_degraded", "lucpfip_ha_total",
  #                        "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total",
  #                        "lucfip_ha_30th", "lucfip_ha_60th", "lucfip_ha_90th", 
  #                        "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th")
  
  # # retirer cette ligne a la fin, mais pour l'instant fait gagner du temps, on a pas besoin de toutes les autres
  # outcome_variables <- c("lucpfip_rapidslow_pixelcount_total", "lucpfsmp_pixelcount_total")
  # 
  # for(voi in outcome_variables){
  #   ## different lags
  #   for(lag in c(1:4)){
  #     parcels <- dplyr::arrange(parcels, lonlat, year)
  #     parcels <- DataCombine::slide(parcels,
  #                                   Var = voi, 
  #                                   TimeVar = "year",
  #                                   GroupVar = "lonlat",
  #                                   NewVar = paste0(voi,"_lag",lag),
  #                                   slideBy = -lag, 
  #                                   keepInvalid = TRUE)
  #     parcels <- dplyr::arrange(parcels, lonlat, year)
  #   }
  #   
  #   
  #   for(py in c(2,3,4)){
  #     ## Past-year average (2, 3 and 4 years) 
  #     parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
  #     parcels[is.nan(parcels$newv),"newv"] <- NA
  #     colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_",py,"pya")
  #     
  #     # Lag it
  #     # note that 3pya_lag1 is different from 4pya. 
  #     parcels <- dplyr::arrange(parcels, lonlat, year)
  #     parcels <- DataCombine::slide(parcels,
  #                                   Var = paste0(voi,"_",py,"pya"), 
  #                                   TimeVar = "year",
  #                                   GroupVar = "lonlat",
  #                                   NewVar = paste0(voi,"_",py,"pya_lag1"),
  #                                   slideBy = -1, 
  #                                   keepInvalid = TRUE)  
  #     parcels <- dplyr::arrange(parcels, lonlat, year)
  #     
  #     ## YOY growth rate
  #     parcels <- mutate(parcels,
  #                       !!as.symbol(paste0(voi,"_yoyg")) := 100*(!!as.symbol(paste0(voi)) - 
  #                                                                !!as.symbol(paste0(voi,"_lag1"))) /
  #                                                                !!as.symbol(paste0(voi,"_lag1")))
  #     # lag it
  #     parcels <- dplyr::arrange(parcels, lonlat, year)
  #     parcels <- DataCombine::slide(parcels,
  #                                   Var = paste0(voi,"_yoyg"), 
  #                                   TimeVar = "year",
  #                                   GroupVar = "lonlat",
  #                                   NewVar = paste0(voi,"_yoyg_lag1"),
  #                                   slideBy = -1, 
  #                                   keepInvalid = TRUE)  
  #     parcels <- dplyr::arrange(parcels, lonlat, year)
  #   }
  # }
  
  ### ADD BASELINE FOREST EXTENT 
  # this is a cross section, computed in prepare_2000_forest_extents.R
  bfe <- readRDS(file.path(paste0("temp_data/processed_parcels/baseline_fc_cs_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"km_IBS_CR.rds")))
  
  bfe <- dplyr::select(bfe, -lat, -lon, -idncrs_lat, -idncrs_lon)
  
  if(!all.equal(sort(unique(bfe$lonlat)), sort(unique(parcels$lonlat)))){stop("baseline forest extent and main data frames do not have the same set of grid cells")}  
  # merge it with our parcel panel
  parcels <- left_join(parcels, bfe, by = "lonlat")
  
  rm(bfe)
  # test that the lonlat have been attributed to parcels equally in the two processes 
  # (prepare_lucpfip.R and prepare_2000_forest_extents.R)
  # 
  # parcels2 <- base::merge(parcels, bfe, by = c("lonlat", "lat","lon"))
  # setorder(parcels1, lonlat, year)
  # setorder(parcels2, lonlat, year)
  # row.names(parcels1) <- NULL
  # row.names(parcels2) <- NULL
  # all.equal(st_drop_geometry(parcels1[,c("lonlat", "lon","lat")]), 
  #           st_drop_geometry(parcels2[,c("lonlat", "lon","lat")]))
  # 
  # all(names(parcels1)==names(parcels2))
  # all.equal(st_drop_geometry(parcels1), st_drop_geometry(parcels2))
  # rm(parcels2)
  
  
  ### COMPUTE ESTIMATED ANNUAL FOREST REMAINING
  
  if(file.exists(file.path(paste0("temp_data/processed_parcels/remaining_forest_panel_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"CR.rds")))){
    
    remaining <- readRDS(file.path(paste0("temp_data/processed_parcels/remaining_forest_panel_",
                                          parcel_size/1000,"km_",
                                          catchment_radius/1000,"CR.rds")))
    
    if(nrow(parcels)!=nrow(remaining)){stop("remaining does not have the same number of parcels as main panel")}
    parcels <- inner_join(parcels, remaining, by = c("lonlat", "year"))
    
  } else { 
  
    # anyNA(parcels$lucpfap_pixelcount) returns FALSE
    
    # then annual lucfp accumulated over past years
    
    year_list <- list()
    
    # in the first year (2001), the past year accumulated lucfp is null. 
    year_list[["2001"]] <- parcels[parcels$year == 2001, c("lonlat", "year")] 
    # names(year_list[["2001"]]) <- "lonlat"%>% as.data.frame() 
    year_list[["2001"]][,"accu_lucfp_since2k"] <- 0
    year_list[["2001"]][,"accu_lucpfp_since2k"] <- 0
    
    # then, each year's lucfp accumulated in the past is the sum of *past years'* lucpfap_pixelcount
    years <- 2002:max(parcels$year)
    for(y in years){
      sub_ <- parcels[parcels$year < y,]
      year_list[[as.character(y)]] <- ddply(sub_, "lonlat", summarise,
                                              accu_lucfp_since2k = sum(lucfap_pixelcount, na.rm = TRUE),
                                              accu_lucpfp_since2k = sum(lucpfap_pixelcount, na.rm = TRUE))
      year_list[[as.character(y)]][,"year"] <- y
    }
    
    # data <- parcels[parcels$year < y,]
    # t <- ddply(data, "lonlat", summarise,
    #            #past_accu_lucfp = sum(total_lucfp_30th, na.rm = TRUE),
    #            accu_lucpfp_since2k = sum(lucpfap_pixelcount, na.rm = TRUE))
    
    accu_lucfp_df <- bind_rows(year_list)
    
    parcels <- inner_join(parcels, accu_lucfp_df, by = c("lonlat", "year"))
    
  
    # summary(parcels$accu_lucpfp_since2k)
    parcels <- dplyr::mutate(parcels, 
                             remain_f30th_pixelcount = fc2000_30th_pixelcount - accu_lucfp_since2k,
                             remain_pf_pixelcount = pfc2000_total_pixelcount - accu_lucpfp_since2k)
  
    # parcels[parcels$lonlat == 1267,c("lonlat", "year", "lucpfap_pixelcount", "accu_lucpfp_since2k", "remain_pf_pixelcount")] 
    
    
    remaining <- parcels[,c("lonlat", "year", "remain_f30th_pixelcount", "remain_pf_pixelcount", "accu_lucfp_since2k", "accu_lucpfp_since2k" )]
    saveRDS(remaining, file.path(paste0("temp_data/processed_parcels/remaining_forest_panel_",
                                        parcel_size/1000,"km_",
                                        catchment_radius/1000,"CR.rds")))
    
    rm(accu_lucfp_df, year_list)
  }
  
  
  
  ### MAKE THE PARCEL_ID 
  uni_lonlat <- unique(parcels$lonlat)
  parcels <- mutate(parcels, 
                    parcel_id = match(lonlat, uni_lonlat)) 
  rm(uni_lonlat)
  
  ## some arrangements
  parcels <- dplyr::arrange(parcels, parcel_id, lonlat, year)
  row.names(parcels) <- seq(1,nrow(parcels))
  
  
  saveRDS(parcels, file.path(paste0("temp_data/panel_parcels_ip_final_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"CR.rds")))  
  
  # this does not write the largest data frames (catchment_radius of 50km)
  # (Erreur : Error in libxlsxwriter: 'Worksheet row or column index out of range.')
  # run it to get the two other ones in xlsx format still. 
  # write_xlsx(parcels, file.path(paste0("temp_data/panel_parcels_ip_final_",
  #                                   parcel_size/1000,"km_",
  #                                   catchment_radius/1000,"CR.xlsx")))  
  
  rm(parcels)
}

PS <- 3000 
catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
for(CR in catchment_radiuseS){
  merge_lhs_rhs(parcel_size = PS, 
                catchment_radius = CR)
}
rm(merge_lhs_rhs)


### Write to Stata if analyses/tests are to be done there...
# parcel_size <- 3000
# catchment_radius <- 3e4
# catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
# for(CR in catchment_radiuseS){
# toc <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
#                                   parcel_size/1000,"km_",
#                                   catchment_radius/1000,"CR.rds")))  
# write.dta(toc, file.path(paste0("temp_data/panel_parcels_ip_final_",
#                                 parcel_size/1000,"km_",
#                                 catchment_radius/1000,"CR.dta")))
# }



#### TESTING ZONE #### 
# names(LHS)  
# names(RHS)  
# 
# length(unique(LHS$lonlat))
# length(unique(RHS_2001$lonlat))
# 
# # LHS was not ordered as expected, hence 
# RHS_2001 <- RHS[RHS$year >= 2001,]
# all.equal(RHS_2001$lonlat, LHS$lonlat) # returns FALSE
# 
# LHS_ordered <- dplyr::arrange(LHS, lonlat, year)
# all.equal(LHS_ordered, LHS)
# # RHS was ordered as expected
# RHS_ordered <- dplyr::arrange(RHS, lonlat, year)
# all.equal(RHS_ordered, RHS)
# 
# # once reordered, we indeed have the same set of parcels in both data frames. 
# all.equal(RHS_2001$lonlat, LHS_ordered$lonlat)
# all.equal(RHS_2001[,c("lat", "lon")], LHS_ordered[,c("lat", "lon")]) 
# coordinates match too (only the row names are different). 



