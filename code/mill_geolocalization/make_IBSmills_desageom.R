######################################################################
#                                                                    #
#   Give IBS mills a desa geometry                                   #
#                                                                    #
#   Input:  - IBS mills, each of the year it has the most recent     #
#             valid desa id                                          #
#             --> IBSmills_valid_desa.dta                            #
#                                                                    #
#           - crosswalked annual village shapefiles and 2010 map:    #
#             --> desa_1998_crosswalked.shp to 2009                  #
#             --> indo_by_desa_2010.shp                              #
#                                                                    #
#                                                                    #
#   Output: - IBS mill with the smallest desa geometry possible      #
#             (for further UML matching works)                       #
#           --> IBSmills_desageom.Rdata                              #
#
#           - IBS mills with smallest desa geometry's centroid
#             (for those that could not be matched)
#           --> IBSmills_desacentro.dta                             
#   
######################################################################
######################################################################

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "readxl", "foreign", "readstata13",
                   "sf", "rgdal")
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


#######################################################################


# read IBS mills
ibs_desa <- read.dta13(file.path("temp_data/processed_mill_geolocalization/IBSmills_valid_desa.dta"))


################################# ASSOCIATE DESA GEOMETRIES TO IBS DESA_ID ###########################################
#templates:
ibs_desa[,"geom"] <- c(1:nrow(ibs_desa))
#annual_n_badcodes <- data.frame(year = c(1998:2010))
#annual_n_badcodes$n_obs <- c(1998:2010)


t <- 1998
while(t < 2010){
  #read year t desa shapefile
  annualpolypath <- paste0("input_data/indonesia_spatial/village_shapefiles/desa_",t,"_crosswalked.shp")
  annual_desa_poly <- st_read(file.path(annualpolypath))

  # give the column in "using" the same name as in master
  names(annual_desa_poly)[names(annual_desa_poly) == paste0("d__",t)] <- "desa_id"

  mills_t <- filter(ibs_desa, year == t)

  #bc <- 0
  for(i in mills_t$firm_id){
    # desa code of the mill from year t.
    desa_of_i <- mills_t[mills_t$firm_id == i, "desa_id"]

    # some desa codes are not found in village maps (holes in crosswalk probably, or misreported desa codes in ibs*)
    if(nrow(annual_desa_poly[annual_desa_poly$desa_id == desa_of_i, "geometry"])!= 0){
      # give the corresponding desa geometry.
        ibs_desa[ibs_desa$firm_id == i, "geom"] <- annual_desa_poly[annual_desa_poly$desa_id == desa_of_i, "geometry"]

    } else {
      ibs_desa[ibs_desa$firm_id == i, "geom"] <- NA

      #bc <- bc + 1
    }
  }
  #accounts for mismatches peryear
  #annual_n_badcodes[annual_n_badcodes$year == t, "n_obs"] <- bc
  rm(annual_desa_poly)
  t <- t+1
}


### 2010
desa_map_2010 <- st_read(file.path("input_data/indonesia_spatial/village_shapefiles/desa_map_2010/indo_by_desa_2010.shp"))
names(desa_map_2010)[names(desa_map_2010) == "IDSP2010"] <- "desa_id"
mills_t <- filter(ibs_desa, year == 2010)
bc <- 0
ndu_desa <- 0
#ibs_desa$from_du_desa <- c(1:nrow(ibs_desa))*0
for(i in mills_t$firm_id){
  # desa code of the mill from year t.
  desa_of_i <- mills_t[mills_t$firm_id == i, "desa_id"]

  # All desa codes are found in village maps (no holes issue since this is the actual map of 2010 - not crosswalked). But pb is that
  # there are duplicates in IDSP2010 so that several polygons may match with one desa_id.
  if(nrow(desa_map_2010[desa_map_2010$desa_id == desa_of_i, "geometry"]) == 1){
    # give the corresponding desa geometry.
    ibs_desa[ibs_desa$firm_id == i, "geom"] <- desa_map_2010[desa_map_2010$desa_id == desa_of_i, "geometry"]

  } else if(nrow(desa_map_2010[desa_map_2010$desa_id == desa_of_i, "geometry"]) > 1){
    ibs_desa[ibs_desa$firm_id == i, "geom"] <- NA
    ndu_desa <- ndu_desa + 1
    #ibs_desa[ibs_desa$firm_id == i, "from_du_desa"] <- 1
  } else if(nrow(desa_map_2010[desa_map_2010$desa_id == desa_of_i, "geometry"]) == 0){
    ibs_desa[ibs_desa$firm_id == i, "geom"] <- NA
    bc <- bc + 1
  }
}
#ibs_desa[ibs_desa$from_du_desa == 1,]
#accounts for mismatches peryear
#annual_n_badcodes[annual_n_badcodes$year == 2010, "n_obs"] <- bc

# * we don't try to match these desa codes with n-1 or n-2 years of village maps because there is a risk that these codes correspond to
# different areas in other years.

saveRDS(ibs_desa, file = "temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata")

### all those with a desa geometry, we attribute them desa centroid coordinates
ibs_desa_sf <- st_as_sf(ibs_desa[!is.na(ibs_desa$geom),])
# Accounting for the likelihood that a mill may also be not in the largest polygon, we set the option to FALSE
ibs_desa_centro <- st_centroid(ibs_desa_sf, of_largest_polygon = FALSE)

# of_largest_polygon option does not change much, as multi parts are not very far from each others. 
# ibs_desa_centro2 <- st_centroid(ibs_desa_sf, of_largest_polygon = TRUE)
# diff_centro_poly <- ibs_desa_sf[ibs_desa_centro$geom != ibs_desa_centro2$geom,]
# diff_centro <- st_centroid(diff_centro_poly, of_largest_polygon = FALSE)
# diff_centro2 <- st_centroid(diff_centro_poly, of_largest_polygon = TRUE)
# plot(diff_centro_poly$geom[4])
# plot(diff_centro1$geom, add = TRUE, col = "green")
# plot(diff_centro2$geom, add = TRUE, col = "red")

ibs_desa_centro1 <- cbind(ibs_desa_centro$firm_id, st_coordinates(ibs_desa_centro)) 
colnames(ibs_desa_centro1) <- c("firm_id", "lon", "lat")
ibs_desa_centro1 <- as.data.frame(ibs_desa_centro1)

write.dta(ibs_desa_centro1, "temp_data/processed_mill_geolocalization/IBSmills_desacentro.dta")
