##### prepare shapefile of Indonesian Islands. 


### PACKAGES ###

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "sf", "foreign", "readstata13")

# Install them in their project-specific versions
renv::restore(packages = neededPackages)

# Load them
lapply(neededPackages, library, character.only = TRUE)

### /!\ IF renv::restore() FAILS TO INSTALL SOME PACKAGES FROM neededPackages /!\ ### 

# 1. Remove these package names from neededPackages above, and rerun renv::restore(packages = neededPackages)
# 2. Write them in troublePackages below, uncomment, and run the following code chunk: 
  
  # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
  # troublePackages <- c("")
  # # Attempt to load packages from user's default libraries. 
  # lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ‘’") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### Aggregate Indonesian province shapes in islands. ####
provinces <- st_read("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp")
provinces <- dplyr::select(provinces, NAME_1)

provinces$island <- NA

provinces$island[provinces$NAME_1 == "Aceh" |
                        provinces$NAME_1 == "Bangka-Belitung" |
                        provinces$NAME_1 == "Bengkulu" |
                        provinces$NAME_1 == "Jambi" |
                        provinces$NAME_1 == "Kepulauan Riau" |
                        provinces$NAME_1 == "Lampung" |
                        provinces$NAME_1 == "Riau" |
                        provinces$NAME_1 == "Sumatera Barat" |
                        provinces$NAME_1 == "Sumatera Selatan" |
                        provinces$NAME_1 == "Sumatera Utara" ] <- "Sumatra"

provinces$island[provinces$NAME_1 == "Kalimantan Barat" |
                        provinces$NAME_1 == "Kalimantan Selatan" |
                        provinces$NAME_1 == "Kalimantan Tengah" |
                        provinces$NAME_1 == "Kalimantan Timur" |
                        provinces$NAME_1 == "Kalimantan Utara" ] <- "Kalimantan"

provinces$island[provinces$NAME_1 == "Papua" |
                        provinces$NAME_1 == "Irian Jaya Barat" ] <- "Papua"


# provinces$island[provinces$NAME_1 == "Gorontalo" |
#                    provinces$NAME_1 == "Sulawesi Selatan" |
#                    provinces$NAME_1 == "Sulawesi Tengah" |
#                    provinces$NAME_1 == "Sulawesi Barat" |
#                    provinces$NAME_1 == "Sulawesi Utara" |
#                    provinces$NAME_1 == "Sulawesi Tenggara"] <- "Sulawesi"

island_sf <- provinces[!is.na(provinces$island),c("island", "geometry")]

IslandS <- c("Sumatra", "Kalimantan", "Papua")#, "Sulawesi"
for(Island in IslandS){
  island_sf$geometry[island_sf$island == Island] <- st_union(island_sf$geometry[island_sf$island == Island])
}

island_sf <- island_sf[!duplicated(island_sf$island),]

## export 
dir.create("temp_data/processed_indonesia_spatial/island_sf")
st_write(island_sf, "temp_data/processed_indonesia_spatial/island_sf", driver = "ESRI Shapefile", delete_dsn = TRUE)


#### Join 2000 district shapes and names ####
## shapes
district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
# there is no dupicate in the district shapes. However, 
# the _new category includes several places (3 in Sumatra and 3 in Java)... which are actually lakes! 
# so remove them
district_sf <- district_sf[district_sf$d__2000!="_new",]

## names
district_names <- read.dta13(file.path("temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta"))
# since district_sf is in 2000 polygons, we should us names from 2000! 
district_names2000 <- district_names[district_names$year == 2000, c("year", "name_", "bps_")]
district_names2000[,c("name_","bps_")] %>% anyDuplicated() # none is a duplicate in both name and code
district_names2000[,c("name_")] %>% anyDuplicated() # neither in name only
district_names2000[,c("bps_")] %>% anyDuplicated() # nor in code only (duplicated codes in years
# have been removed in prepare_crosswalks.do)

## join names to shapes
district_sf$d__2000 <- district_sf$d__2000 %>% as.character()
district_names2000$bps_ <- district_names2000$bps_ %>% as.character()
nrow(district_names2000) == nrow(district_sf)

district_sf <- left_join(x = district_sf, y = district_names2000[,c("name_", "bps_")], 
                         by = c("d__2000" = "bps_"),
                         all = FALSE, all.x = FALSE, all.y = FALSE)

dir.create("temp_data/processed_indonesia_spatial/district2000_sf")
st_write(district_sf, "temp_data/processed_indonesia_spatial/district2000_sf", 
         driver = "ESRI Shapefile", delete_dsn = TRUE)

rm(district_names, district_names2000, district_sf)

### GIVE ISLAND VARIABLE IN DISTRICT SF 
# district_sf <- st_read(file.path("input_data/indonesia_spatial/district_shapefiles/district_2015_base2000.shp"))
# district_names <- read.dta13(file.path("temp_data/processed_indonesia_spatial/province_district_code_names_93_2016.dta"))
# district_names <- district_names[!duplicated(district_names$bps_),]
# district_sf$d__2000 <- district_sf$d__2000 %>% as.character()
# district_names$bps_ <- district_names$bps_ %>% as.character()
# district_sf <- left_join(x = district_sf, y = district_names[,c("name_", "bps_")], 
#                          by = c("d__2000" = "bps_"),
#                          all = FALSE, all.x = FALSE, all.y = FALSE)
# district_sf_prj <- st_transform(district_sf, crs = indonesian_crs)
# 
# island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
# names(island_sf)[names(island_sf)=="island"] <- "shape_des"
# island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)
# 
# districts_sgbp <- st_within(district_sf_prj, island_sf_prj)
# copy <- districts_sgbp
# # all districts are covered by only one island
# districts_sgbp[lengths(districts_sgbp)>1]
# # or another island than the three islands of interest. Attribute them 0  
# districts_sgbp[lengths(districts_sgbp)==0] <- 0
# # 55 districts are in one of our three islands of interest
# unlist(districts_sgbp)[unlist(districts_sgbp)!=0] %>% length()
# 
# unique(unlist(districts_sgbp))
# district_sf$island <- unlist(districts_sgbp) %>% as.character()
# district_sf$island <- island_sf$shape_des[match(district_sf$island, row.names(island_sf))]
# 
# district_sf_within <- district_sf[!is.na(district_sf$island),]
# 
# # with alternative methods
# districts_sgbp_overlaps <- st_overlaps(district_sf_prj, island_sf_prj)
# districts_sgbp_overlaps[lengths(districts_sgbp_overlaps)==0] <- 0
# unlist(districts_sgbp_overlaps)[unlist(districts_sgbp_overlaps)!=0] %>% length() # 96 districts 
# 
# district_sf$island <- unlist(districts_sgbp_overlaps) %>% as.character()
# district_sf$island <- island_sf$shape_des[match(district_sf$island, row.names(island_sf))]
# 
# district_sf_overla <- district_sf[!is.na(district_sf$island),]
# 
# plot(st_geometry(district_sf_overla))
# plot(st_geometry(district_sf_within))
