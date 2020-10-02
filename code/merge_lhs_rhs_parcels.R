### Merge panel data frames of parcels with LHS and RHS variables 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##### 0. PACKAGES, WD, OBJECTS #####

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "foreign", 
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
# keep only year before 2015 (after they mean nothing since we plantation data are from 2015)
lucpfip <- lucpfip[lucpfip$year<=2015,] # now runs from 2001-1998
lucfip <- lucfip[lucfip$year<=2015,] # now runs from 2001-1998

# remove coordinates, they are already in RHS
lucpfip <- dplyr::select(lucpfip, -lat, -lon)
lucfip <- dplyr::select(lucfip, -lat, -lon)


nrow(lucpfip)==nrow(lucfip)
LHS <- base::merge(lucpfip, lucfip, by = c("parcel_id", "year"))

# explicative variables (runs from 1998-2015)
RHS <-  readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                                  parcel_size/1000,"km_",catchment_radius/1000,"CR.rds")))

# MERGE
# years 1998 - 2000 from RHS will not match, we don't need to keep them because the information 
# from these years is captured in add_parcel_variables.R within lag variables. Hence all = FALSE
parcels <- base::merge(LHS, RHS, by = c("parcel_id", "year"), all = FALSE)  




### OUTCOME VARIABLE TIME DYNAMICS

# outcome_variables <- c("lucpfip_ha_intact", "lucpfip_ha_degraded", "lucpfip_ha_total",
#                        "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount_total",
#                        "lucfip_ha_30th", "lucfip_ha_60th", "lucfip_ha_90th", 
#                        "lucfip_pixelcount_30th", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th")

# retirer cette ligne à la fin, mais pour l'instant ça fait gagner du temps, on a pas besoin de toutes les autres
outcome_variables <- "lucpfip_pixelcount_total"

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

}

PS <- 3000 
catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
for(CR in catchment_radiuseS){
  merge_lhs_rhs(parcel_size = PS, 
                catchment_radius = CR)
}


# Write to Stata if analyses/tests are to be done there...
parcel_size <- 3000
catchment_radius <- 3e4
catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
for(CR in catchment_radiuseS){
toc <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"CR.rds")))  
write.dta(toc, file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                catchment_radius/1000,"CR.dta")))
}
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



