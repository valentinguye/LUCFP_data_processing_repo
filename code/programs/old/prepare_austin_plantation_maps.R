######################################################################
#                                                                    #
#   Palm oil remote sensing data on Desa level                        #
#                                                                    #
#   Input: Palm oil remote sensing raw data                          #
#   Output: Palm oil remote sensing processed data                   #   
#        This file has !(a larger extent and) values in every cell   #
######################################################################
######################################################################

#### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
#### OR CALLED FROM LUCFP PROJECT master.do FILE.
#### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES 
# Installs all the packages required in this project, if not already installed in LUCFP/data_processing/renv/library
renv::restore()

# These are the packages needed in this particular script. 
neededPackages = c("dplyr", "raster", "sf", "foreign", "sp", "data.table",
                   "rgdal")

# Load them
lapply(neededPackages, library, character.only = TRUE)

#######################################################################


dir.create("temp_data/processed_plantation_maps")


######################################################################
######################################################################
### LOAD NECESSARY PALM OIL REMOTE SENSING DATA 

# Load palm oil remote sensing data (raster)
tif.all.names     <- list.files(path = "input_data/austin_plantation_maps/IIASA_indo_oilpalm_map", pattern = ".tif")
other.names       <- list.files(path = "input_data/austin_plantation_maps/IIASA_indo_oilpalm_map", pattern = ".tif.")
tif.names         <- tif.all.names[!(tif.all.names %in% other.names)]
data.file         <- list()

for (i in 1: length(tif.names)){
  
  file.path      <- paste0("input_data/austin_plantation_maps/IIASA_indo_oilpalm_map/", tif.names[i]) %>% file.path()
  data.file[[i]] <- raster(file.path)
  
  # Change name of raster to identify year
  names(data.file[[i]]) <- tif.names[i]
}

# Check CRS of rasters
crs(data.file[[1]]) #WGS84

##################################################################
### PREPARE AND EXPORT DATA SET

# Transform raster layers such that all cells have a value (0 or 1)
for (i in 1:length(data.file)){
  
  # # Raster value
  # val.raster <- getValues(data.file[[i]])
  # 
  # # Find empty cells
  # IND.empty <- is.na(val.raster)
  # 
  # # Create binary vector
  # new.value             <- rep(1, length(val.raster))
  # new.value[IND.empty]  <- 0
  # 
  # # Replace original raster values with binary cell vector
  # data.file[[i]] <-  setValues(data.file[[i]], values = new.value)
  
  # Save raster files
  name.vec <- file.path(paste0("temp_data/processed_plantation_maps/reclassified_" ,names(data.file[[i]])))
  #writeRaster(x = data.file[[i]], filename = name.vec, overwrite=TRUE, datatype = "INT1U") 
  raster::reclassify(data.file[[i]], 
                     rcl = cbind(NA,0), 
                     filename = name.vec, 
                     overwrite=TRUE, 
                     datatype = "INT1U")
}

##################################################################