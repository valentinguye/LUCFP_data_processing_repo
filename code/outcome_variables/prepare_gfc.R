### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#   This script downloads and prepares (extracts and thresholds) GFC data, using package gfcanalysis.
#    
#   Input: - Indonesia shape file
#         ---> input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp
#   
#          - INTERNET CONNEXION FOR GFC DATA TILES TO BE DOWNLOADED (Hansen et al. 2013)
#    
#   Output: gfc_data_Indonesia_30th.tif
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### PACKAGES ###

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "sp", "sf", "gfcanalysis")

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



### NEW FOLDERS USED IN THIS SCRIPT 

# this is where GFC tiles are stored 
dir.create("temp_data/GFC_tiles")

# this is where will be stored temporary raster function outputs. 
dir.create("temp_data/processed_lu")


### Define Area of Interest (aoi): Indonesia
provinces <- st_read("input_data/indonesia_spatial/province_shapefiles/IDN_adm1.shp")
provinces <- provinces$geometry

indonesian_extent <- st_as_sfc(st_bbox(provinces))

aoi_sp <- as(indonesian_extent, "Spatial")

rm(provinces, indonesian_extent)

### Download GFC data

#define where all tiles are going to be stored
data_folder <- file.path("temp_data/GFC_tiles")

# Calculate tiles needed to cover the AOI
tiles <- calc_gfc_tiles(aoi_sp)
# length(tiles) 

# version of GFC used here.
gfc_version <- "GFC-2018-v1.6"
# //!\\ this script is not written flexibly to adjust for other versions of GFC. One should check every "18" entries in this script for instance. 

#download tiles - with all layers otherwise later extract_gfc does not work
download_tiles(tiles, 
               output_folder = data_folder, 
               images = c("treecover2000", "lossyear", "gain", "datamask"), 
               dataset = gfc_version)

rm(tiles)


### extract GFC data 

# (can only extract all layers with default stack=change)
# to better understand extract_gfc see https://rdrr.io/cran/gfcanalysis/src/R/extract_gfc.R
extract_gfc(aoi_sp, data_folder,
            stack = "change",
            to_UTM = FALSE,
            dataset = gfc_version,
            filename = file.path(paste0("temp_data/processed_lu/gfc_data_Indonesia.tif")),
            overwrite = TRUE )

rm(aoi_sp)

### Threshold GFC data 

# The forest loss layer is based on a forest definition. This requires to specify a 
# pixel-level canopy cover percentage as the threshold between forest and non-forest state in 2000.
# Here, this is 30%, 60% and 90%.  
# function threshold_gfc is already defined in package gfc_analysis

gfc_data <- brick(file.path(paste0("temp_data/processed_lu/gfc_data_Indonesia.tif")))

# 30%
threshold_gfc(gfc_data,
              forest_threshold=30,
              filename= file.path("temp_data/processed_lu/gfc_data_Indonesia_30th.tif"),
              overwrite = TRUE)
# 60%
threshold_gfc(gfc_data,
              forest_threshold=60,
              filename= file.path("temp_data/processed_lu/gfc_data_Indonesia_60th.tif"),
              overwrite = TRUE)
# 90%
threshold_gfc(gfc_data,
              forest_threshold=90,
              filename= file.path("temp_data/processed_lu/gfc_data_Indonesia_90th.tif"),
              overwrite = TRUE)

rm(gfc_data)