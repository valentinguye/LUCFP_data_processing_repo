
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