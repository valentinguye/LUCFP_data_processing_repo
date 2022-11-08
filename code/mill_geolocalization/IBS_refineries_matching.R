######################################################################
#                                                                    #
#   Match IBS mill dataset with Trase refineries                     #
#                                                                    #
#   Input:  - IBS mills, each of the year it has the most recent     #
#             valid desa id, with the corresponding geometries       #
#             --> IBSmills_desageom.Rdata                            #
#                                                                    #
#           - Trase georeferenced refineries                         #
#             --> trase_refineries.xlsx                              #
#                                                                    #
#   Output: several subsets of IBS, depending on how they            #
#           matched with refineries.                                 #
#   
######################################################################
######################################################################

#### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### TEMP SOLUTION TO RELOCATE WD OF RFSFOOD PROJECT FROM GLOBILUC R PROJECT (AND THUS REPO)
setwd(dir = "C:/Users/GUYE/Desktop/LUCFP/data_processing")
# (necessaire de changer le wd avant de loader here)

# Create a directory to store outputs from this script on REFINERIES, in isolation from the outputs from the mill geolocalization work. 
dir.create("temp_data/processed_refinery_geolocalization") 

#### OR CALLED FROM LUCFP PROJECT master.do FILE.
#### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 

### PACKAGES 
# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "readxl", "foreign", "readstata13",
                   "sf", "rgdal", 
                   "here") 

renv::restore(packages = neededPackages)

# Load them
lapply(neededPackages, library, character.only = TRUE)

here()
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### ### ### 
# READ AND PREPARE DATA 

# read UML mills
# uml <- read_excel(file.path("input_data/uml/traseMills_capEstyear.xlsx"))

# read refineries
refine <- read.csv(here("input_data", "trase_refineries.csv"))
# read IBS mills - note that this is a cross section 
ibs <- readRDS(file.path("temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata"))
nrow(ibs[!duplicated(ibs$firm_id),]) == nrow(ibs)

# Remove IBS manufactories that are matched with Universal Mill List (UML)
ibs_uml <- read.dta13(here("temp_data/IBS_UML_panel_final.dta"))

ibs_uml_cs <- ibs_uml[!duplicated(ibs_uml$firm_id), c("firm_id", "mill_name", "uml_matched_sample", "is_mill") ]

sum(ibs_uml_cs$uml_matched_sample)

ibs_NO_uml <- left_join(ibs, ibs_uml_cs, by = "firm_id")

# among the 17 firm_id that are manually changed in merging_geolocalization_works.do, 
# 13 have a desa geom and thus, are present in IBSmills_desageom.Rdata BUT NOT in IBS_UML_panel_final.dta  
# but we have identified them as being UML mills already (so we remove them here)
nrow(ibs_NO_uml[is.na(ibs_NO_uml$uml_matched_sample),]) 
# --> KEEP only those firm_id that were NOT matched with UML, be they in the IBS_UML_panel_final (==0) or not (NA)
ibs_NO_uml <- dplyr::filter(ibs_NO_uml, ibs_NO_uml$uml_matched_sample == 0 & !is.na(ibs_NO_uml$uml_matched_sample)) 

unique(ibs_NO_uml$mill_name)
# mill_name is actually useless (it has been wiped out when a firm_id was not UML)
ibs_NO_uml <- dplyr::select(ibs_NO_uml, -mill_name)

rm(ibs)

### Make spatial 

## Indonesian CRS
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

#prepare the sfc of ibs mills and their most recently known villages. 
ibs_na <- filter(ibs_NO_uml, is.na(geom))
ibs_NO_uml <- filter(ibs_NO_uml, is.na(geom)== F)   

ibs_NO_uml <- st_as_sf(ibs_NO_uml, crs = 4326)
ibs_NO_uml <- st_transform(ibs_NO_uml, crs = indonesian_crs)

ibs_NO_uml

# prepare UML dataset
refine$latitude <- as.numeric(refine$latitude)
refine$longitude <- as.numeric(refine$longitude)
refine$lat <- refine$latitude
refine$lon <- refine$longitude
refine <- st_as_sf(refine,	coords	=	c("longitude",	"latitude"), crs = 4326)
refine <- st_transform(refine, crs = indonesian_crs)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### MATCH IBS AND UML MILLS ####

ibs_j <- st_join(x = ibs_NO_uml, y = refine, join = st_contains, left = T) # is equivalent to ibs_j <- st_join(x = ibs_NO_uml, y = refine) bc default is st_intersect.

#plot(st_geometry(ibs_NO_uml[ibs_NO_uml$firm_id == 1761,]))
#plot(st_geometry(uml[uml$mill_name == "PT. Socfin - Seumanyam",]), add = TRUE, col = "red")


########################################### IDENTIFY DIFFERENTLY MATCHED SUBSETS ######################################################

# variable lon (or lat) comes from trase_refineries data 
no_nas <- ibs_j[is.na(ibs_j$lon)==F,] # 76 non-missings
nrow(no_nas[unique(no_nas$firm_id),])
# 45 IBS "manufactories" are in a village where at least one Trase refinery stands. 

pot_otm <- no_nas[!(duplicated(no_nas$firm_id) | duplicated(no_nas$firm_id, fromLast = TRUE)), ] 
# these are the 29 potential one-to-many matches (the IBS manufactory is the sole one from IBS (we know the desa of) to be in this village and 
# there may be several Trase refineries in this same village.). 
# Among them, there are:
oto <- pot_otm[!(duplicated(pot_otm$company) | duplicated(pot_otm$company, fromLast = TRUE)), ] 
# 10 IBS manufactories we know the desa of, are the only IBS manufactories in their villages while there is only one Trase refinery in each village.  

otm <- pot_otm[duplicated(pot_otm$company) | duplicated(pot_otm$company, fromLast = TRUE), ]
# 19 IBS mills with identified desa in which they are the only IBS manufactories, but where there is more than one Trase refinery. 
nrow(distinct(otm, lat, lon)) # 12 Trase refineries are involved.

du <- no_nas[(duplicated(no_nas$firm_id) | duplicated(no_nas$firm_id, fromLast = TRUE)), ] 
# these are the 47 IBS manufactories with identified desa that are in the same desa with at least another IBS one.
# These villages with many IBS may encompass either one or many Trase refineries. 

noto <- rbind(otm, du)
# these are the 66 IBS manufactories (19 + 47) that have an identified desa and are either m:m, o:m or m:o with Trase refineries. 

total_potential <- no_nas[!duplicated(no_nas$lon),]
# these are the 45 manufactories of which we know the coordinates, and that are within the villages of IBS manufactories. 


## So, the maximum number of pre-2011 IBS manufactories we can hope to geolocalize with this desa matching technique is 45
## The minimum number is 10 (unless we double-check them with workers and some don't match). 


#### Those that have a desa polygon but match with no Trase refinery (414)
# It is useless to try to find them manually with the directory number of workers and uml's list. 
# we can still find their names with directories and google-search them, and/or directly spot them manually within their villages. 
ibs_unref <- left_join(x = ibs_NO_uml, y = st_set_geometry(oto[,c("firm_id","lat")], NULL), by = "firm_id")
ibs_unref <- filter(ibs_unref, is.na(ibs_unref$lat))

ibs_unref <- left_join(x = ibs_unref, y = st_set_geometry(noto[,c("firm_id","lon")], NULL), by = "firm_id")
ibs_unref <- filter(ibs_unref, is.na(ibs_unref$lon))

ibs_unref <-dplyr::select(ibs_unref, -lon, -lat)
# (the use of "lat" and "lon" was just an arbitrary choice of non-empty variables to flag oto and noto resp.)

# They are in 455 different polygons, of which 359 do not intersect with another one
# (intersections but not equal likely when there is a split and a mill is associated with the polygon of the village before, and one or more 
#other mills are associated with the villages after the split.)
#ibs_unref$grp = sapply(st_equals(ibs_unref), max)
#ibs_unref <- arrange(ibs_unref, grp, firm_id)
#nrow(ibs_unref[unique(ibs_unref$grp),])

#unique_unref <-ibs_unref[unique(ibs_unref$grp),]

#ibs_unref$grp_int = sapply(st_intersects(ibs_unref), max)
#nrow(ibs_unref[unique(ibs_unref$grp_int),])

#plot(st_geometry(ibs_unref[ibs_unref$grp == 1 | ibs_unref$grp == 355,]), add = F)
#plot(st_geometry(ibs_unref[ibs_unref$grp == 355,]), add = F)

# also, this is 473 different year_desa_id, meaning that there are cases where two mills have equal polygons but no equal year_desa_id: these are the 
# cases of two mills being associated to the same polygon but with desa_ids from different years. 
#ibs_unref$year_desa_id <- paste(ibs_unref$year, ibs_unref$desa_id, sep = "")
#nrow(ibs_unref[unique(ibs_unref$year_desa_id),])

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### EXPORT THEM ACCORDINGLY ####

### Those that have no desa polygon though they had an un-flagged desa_id (46)
# they return to stata to get merged with the panel and geolocalized manually. 
ibs_na <-dplyr::select(ibs_na, firm_id)
write.dta(ibs_na, "temp_data/processed_refinery_geolocalization/pre2011_bad_desa_id.dta")


### Those that have a desa polygon but match with no Trase refinery (414)
ibs_unref$desa_id <- as.character(ibs_unref$desa_id)
ibs_unref <- st_transform(ibs_unref, crs = 4326)

## For R 
saveRDS(ibs_unref, file = "temp_data/processed_refinery_geolocalization/ibs_unref.Rdata") 

## For GEE
ibs_unref <-dplyr::select(ibs_unref, firm_id, year, min_year, max_year, industry_code, out_ton_cpo_imp2,
                          workers_total_imp3, avg_in_tot_ton_cpo_imp2, avg_out_ton_cpo_imp2, avg_cpo_price_imp2, 
                          out_ton_rpo_imp1, out_ton_rpo_imp2, out_ton_rpko_imp1, out_ton_rpko_imp2,
                          desa_id, geom)

st_write(ibs_unref, "temp_data/processed_refinery_geolocalization/ibs_unref", driver = "ESRI Shapefile", delete_dsn = TRUE)
#st_write(ibs_unref,"unreferenced_mill_desa.gpkg", driver = "GPKG", delete_dsn = TRUE)

## For Stata 
ibs_unref2 <-dplyr::select(ibs_unref, firm_id)
write.dta(st_set_geometry(ibs_unref2, NULL), "temp_data/processed_refinery_geolocalization/ibs_unref.dta")


#### Those who matched (oto & noto)
## oto
oto <- st_set_geometry(oto, NULL)
oto <-dplyr::select(oto, firm_id, company, company_id, group, ref_id, ref_name, cap_mt, lat, lon)
write.dta(oto, "temp_data/processed_refinery_geolocalization/oto.dta")

## noto 
# sort by desa geometry to ease manual matching
noto$grp = sapply(st_intersects(noto), max)
noto <- arrange(noto, grp, firm_id)
noto <- st_set_geometry(noto, NULL)
noto <-dplyr::select(noto, firm_id, company, company_id, group, ref_id, ref_name, cap_mt, lat, lon, grp)
write.dta(noto, "temp_data/processed_refinery_geolocalization/noto.dta")

## Trase refineries (for GEE)
st_crs(refine) <- 4326
refine <- st_transform(refine, crs = 4326)
refine <-dplyr::select(refine, -lat, -lon)
st_write(refine, "temp_data/processed_refinery_geolocalization/trase_refineries", driver = "ESRI Shapefile", delete_dsn = TRUE)

