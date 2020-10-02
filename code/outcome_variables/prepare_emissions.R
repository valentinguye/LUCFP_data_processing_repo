
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###

# Installs all the packages required in this project, if not already installed in LUCFP/data_processing/renv/library
renv::restore()

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "readstata13", 
                   "raster", "rgdal", "sp", "sf","gfcanalysis",
                   "doParallel", "foreach", "parallel")

# Load them
lapply(neededPackages, library, character.only = TRUE)



### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/processed_emissions")



### RASTER OPTIONS ### 
# Do not change chunksize, as it is not safe in combinatin with clusterR and machine specific. 
rasterOptions(timer = TRUE, 
              progress = "text",
              tmpdir = "temp_data/raster_tmp")

### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"




### Prepare some sample croping polygons 
provinces <- st_read("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp")
provinces <- dplyr::select(provinces, NAME_1)
aoi <- provinces[provinces$NAME_1 == "Irian Jaya Barat",]
aoi_sp <- as(aoi, "Spatial")


### CIFOR PEAT MAPS 
cifor_global_peat <- raster("input_data/peat/cifor_gumbricht_swamp/TROP-SUBTROP_PeatDepthV2_2016_CIFOR.tif")

crop(cifor_global_peat, aoi_sp, 
     filename = file.path("temp_data/test/cifor_peat_westpapua.tif"),
     dataType = "INT1U", 
     overwrite = TRUE)

wpapua <- raster(file.path("temp_data/test/cifor_peat_westpapua.tif"))

# According to the readme provided with the data (TROP-SUBTROP_PeatDepthV2_2016_CIFOR.docx):
# Pixel data represent wetland depths in meters
# 0 = none
# 1-10 = peat depth in meters

# "A validation of our depth map against ground measured peat depths 
# (i.e. soil profiles) suggests that our deepest values (>10m) overestimate depth. 
# For this reason, all depths >10m have been thresholded to 10m."

# However, 
# wpapua %>% getValues %>% unique() 
# returns the value 15 and does not return the value 0 
# when plotting wpapua, 15 is the main value, (covering the sea in particular). 
# and 
# cifor_global_peat %>% sampleRegular(1e6) %>% unique()
# returns 

# therefore we reclassify 15 to 0 
raster::reclassify(wpapua, 
                   rcl = cbind(15, 0), 
                   filename = file.path("temp_data/test/cifor_peat_westpapua_reclassified.tif"), 
                   dataType = "INT1U", 
                   overwrite = TRUE)






##### MOA PEAT MAP 
moa_sea <- st_read("input_data/peat/moa_ritung_peatmap/SEA_Peatland.shp")



























