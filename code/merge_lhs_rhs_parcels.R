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


merge_lhs_rhs <- function(parcel_size, catchment_radius){
  
  # megre lucpfip and lucfip outcome variables data sets together
  lucpfip <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                  parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  lucfip <- readRDS(file.path(paste0("temp_data/processed_parcels/lucfip_panel_",
                                      parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  lucpfsmp <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfsmp_panel_",
                                       parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR.rds")))
  
  # keep only year before 2015 (after they mean nothing since we plantation data are from 2015)
  lucpfip <- lucpfip[lucpfip$year<=2015,] # now runs from 2001-1998
  lucfip <- lucfip[lucfip$year<=2015,] # now runs from 2001-1998
  lucpfsmp <- lucpfsmp[lucpfsmp$year<=2015,] # now runs from 2001-1998
  
  # remove coordinates, they are already in RHS
  lucpfip <- dplyr::select(lucpfip, -lat, -lon)
  lucfip <- dplyr::select(lucfip, -lat, -lon)
  lucpfsmp <- dplyr::select(lucpfsmp, -lat, -lon)
  
  
  nrow(lucpfip)==nrow(lucfip)
  LHS <- base::merge(lucpfip, lucfip, by = c("parcel_id", "year"))
  LHS <- base::merge(LHS, lucpfsmp, by = c("parcel_id", "year"))
  
  # make variable that counts lucfp events on both small and medium sized plantations 
  LHS$lucpfsmp_pixelcount_total <- LHS$lucpfsp_pixelcount_total + LHS$lucpfmp_pixelcount_total
  LHS$lucpfsmp_ha_total <- LHS$lucpfsp_ha_total + LHS$lucpfmp_ha_total
  
  # LHS$lucfsmp_pixelcount_30th <- LHS$lucfsp_pixelcount_30th + LHS$lucfmp_pixelcount_30th
  # LHS$lucfsmp_ha_30th <- LHS$lucfsp_ha_30th + LHS$lucfmp_ha_30th
  
  
  # explicative variables (runs from 1998-2015)
  RHS <-  readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                    parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))
  
  # MERGE
  # years 1998 - 2000 from RHS will not match, we don't need to keep them because the information 
  # from these years is captured in add_parcel_variables.R within lag variables. Hence all = FALSE
  parcels <- base::merge(LHS, RHS, by = c("parcel_id", "year"), all = FALSE)  
  
  rm(LHS, RHS, lucpfip, lucfip, lucpfsmp)
  
  
  ### OUTCOME VARIABLE TIME DYNAMICS
  
  # outcome_variables <- c("lucpfip_ha_intact", "lucpfip_ha_degraded", "lucpfip_ha_total",
  #                        "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total",
  #                        "lucfip_ha_30th", "lucfip_ha_60th", "lucfip_ha_90th", 
  #                        "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th")
  
  # retirer cette ligne a la fin, mais pour l'instant fait gagner du temps, on a pas besoin de toutes les autres
  outcome_variables <- c("lucpfip_pixelcount_total", "lucpfsmp_pixelcount_total", "lucpfmp_pixelcount_total", "lucpfsp_pixelcount_total")
  
  for(voi in outcome_variables){
    ## different lags
    for(lag in c(1:4)){
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = voi, 
                                    TimeVar = "year",
                                    GroupVar = "parcel_id",
                                    NewVar = paste0(voi,"_lag",lag),
                                    slideBy = -lag, 
                                    keepInvalid = TRUE)
      parcels <- dplyr::arrange(parcels, parcel_id, year)
    }
    
    
    for(py in c(2,3,4)){
      ## Past-year average (2, 3 and 4 years) 
      parcels$newv <- rowMeans(x = parcels[,paste0(voi,"_lag",c(1:py))], na.rm = FALSE)
      parcels[is.nan(parcels$newv),"newv"] <- NA
      colnames(parcels)[colnames(parcels)=="newv"] <- paste0(voi,"_",py,"pya")
      
      # Lag it
      # note that 3pya_lag1 is different from 4pya. 
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_",py,"pya"), 
                                    TimeVar = "year",
                                    GroupVar = "parcel_id",
                                    NewVar = paste0(voi,"_",py,"pya_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      
      ## YOY growth rate
      parcels <- mutate(parcels,
                        !!as.symbol(paste0(voi,"_yoyg")) := 100*(!!as.symbol(paste0(voi)) - 
                                                                 !!as.symbol(paste0(voi,"_lag1"))) /
                                                                 !!as.symbol(paste0(voi,"_lag1")))
      # lag it
      parcels <- dplyr::arrange(parcels, parcel_id, year)
      parcels <- DataCombine::slide(parcels,
                                    Var = paste0(voi,"_yoyg"), 
                                    TimeVar = "year",
                                    GroupVar = "parcel_id",
                                    NewVar = paste0(voi,"_yoyg_lag1"),
                                    slideBy = -1, 
                                    keepInvalid = TRUE)  
      parcels <- dplyr::arrange(parcels, parcel_id, year)
    }
  }
  
  ### ADD BASELINE FOREST EXTENT 
  # this is a cross section, computed in prepare_2000_forest_extents.R
  bfe <- readRDS(file.path(paste0("temp_data/processed_parcels/baseline_fc_cs_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"km_IBS_CR.rds")))
  bfe <- dplyr::select(bfe, parcel_id, lat, lon, everything())
  
  # merge it with our parcel panel
  parcels <- base::merge(parcels, 
                         dplyr::select(bfe, -lat, -lon, -fc2000_30th_ha, -pfc2000_total_ha), 
                         by = "parcel_id")
  
  rm(bfe)
  # test that the parcel_id have been attributed to parcels equally in the two processes 
  # (prepare_lucpfip.R and prepare_2000_forest_extents.R)
  # 
  # parcels2 <- base::merge(parcels, bfe, by = c("parcel_id", "lat","lon"))
  # setorder(parcels1, parcel_id, year)
  # setorder(parcels2, parcel_id, year)
  # row.names(parcels1) <- NULL
  # row.names(parcels2) <- NULL
  # all.equal(st_drop_geometry(parcels1[,c("parcel_id", "lon","lat")]), 
  #           st_drop_geometry(parcels2[,c("parcel_id", "lon","lat")]))
  # 
  # all(names(parcels1)==names(parcels2))
  # all.equal(st_drop_geometry(parcels1), st_drop_geometry(parcels2))
  # rm(parcels2)
  
  
  ### COMPUTE ESTIMATED ANNUAL FOREST REMAINING - ONLY FOR PRIMARY AS OF NOW
  
  # first need to compute annual lucfp
  # parcels <- dplyr::mutate(parcels, 
  #                          total_lucfp_30th = lucfip_pixelcount_30th + lucfsmp_pixelcount_30th)
  
  parcels <- dplyr::mutate(parcels, 
                           total_lucpfp_total = lucpfip_pixelcount_total + lucpfsmp_pixelcount_total)
  
  # anyNA(parcels$total_lucpfp_total) returns FALSE
  
  # then annual lucfp accumulated over past years
  
  year_list <- list()
  
  # in the first year (2001), the past year accumulated lucfp is null. 
  year_list[["2001"]] <- parcels[parcels$year == 2001, c("parcel_id", "year")] 
  # names(year_list[["2001"]]) <- "parcel_id"%>% as.data.frame() 
  year_list[["2001"]][,"accu_lucpfp_since2k"] <- 0
  
  # then, each year's lucfp accumulated in the past is the sum of *past years'* total_lucpfp_total
  years <- 2002:max(parcels$year)
  for(y in years){
    sub_ <- parcels[parcels$year < y,]
    year_list[[as.character(y)]] <- ddply(sub_, "parcel_id", summarise,
                                            #accu_lucfp_since2k = sum(total_lucfp_30th, na.rm = TRUE),
                                            accu_lucpfp_since2k = sum(total_lucpfp_total, na.rm = TRUE))
    year_list[[as.character(y)]][,"year"] <- y
  }
  
  # data <- parcels[parcels$year < y,]
  # t <- ddply(data, "parcel_id", summarise,
  #            #past_accu_lucfp = sum(total_lucfp_30th, na.rm = TRUE),
  #            accu_lucpfp_since2k = sum(total_lucpfp_total, na.rm = TRUE))
  
  accu_lucfp_df <- bind_rows(year_list)
  
  parcels <- merge(parcels, accu_lucfp_df, by = c("parcel_id", "year"))
  

  # summary(parcels$accu_lucpfp_since2k)
  parcels <- dplyr::mutate(parcels, 
                           #remain_f30th_pixelcount = fc2000_30th_pixelcount - accu_lucfp_since2k,
                           remain_pf_pixelcount = pfc2000_total_pixelcount - accu_lucpfp_since2k)
  
  # parcels[parcels$parcel_id == 1267,c("parcel_id", "year", "total_lucpfp_total", "accu_lucpfp_since2k", "remain_pf_pixelcount")] 
  
  
  
  
  ## some arrangements
  parcels <- dplyr::arrange(parcels, parcel_id, year)
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
  
  rm(accu_lucfp_df, year_list, parcels)
}

PS <- 3000 
catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
for(CR in catchment_radiuseS){
  merge_lhs_rhs(parcel_size = PS, 
                catchment_radius = CR)
}



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
# length(unique(LHS$parcel_id))
# length(unique(RHS_2001$parcel_id))
# 
# # LHS was not ordered as expected, hence 
# RHS_2001 <- RHS[RHS$year >= 2001,]
# all.equal(RHS_2001$parcel_id, LHS$parcel_id) # returns FALSE
# 
# LHS_ordered <- dplyr::arrange(LHS, parcel_id, year)
# all.equal(LHS_ordered, LHS)
# # RHS was ordered as expected
# RHS_ordered <- dplyr::arrange(RHS, parcel_id, year)
# all.equal(RHS_ordered, RHS)
# 
# # once reordered, we indeed have the same set of parcels in both data frames. 
# all.equal(RHS_2001$parcel_id, LHS_ordered$parcel_id)
# all.equal(RHS_2001[,c("lat", "lon")], LHS_ordered[,c("lat", "lon")]) 
# coordinates match too (only the row names are different). 



