### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# THIS SCRIPT installs (with renv::restore()) project-wide necessary packages. 
# it is called in master.do, but it can (and not must) also be executed separately, before individual R runs.

# This project's data processing reproducibility is based on package version tracking enabled 
# by the package {renv}. 

# This is how it works, briefly, in this project. More info at https://rstudio.github.io/renv/articles/renv.html

# The packages used in this project's data processing are loaded from a project local library renv/library. 
# This library hosts packages in their project-specific versions, hence contributing to data processing reproducibility.   

# These packages and their respective versions are recorded in the renv.lock file provided with the project. 

# renv/library is not provided with the project, but created the first time R is launched from this project, either 
# from master.do, or within R_project_for_individual_runs.Rproj 
# because the project .Rprofile calls renv/activate.R which, in particular, 
# downloads, installs, or updates, if necessary, the renv package; and creates renv/library. 

# Then, upon renv::restore() command, packages recorded in renv.lock are installed in renv/library, if necessary 
# (i.e. if not available from the renv cache, see https://rstudio.github.io/renv/articles/renv.html#cache-1 ) 

# Note that each of this project's R scripts installs the packages it requires, whether the present
# script has been already executed or not. 

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#"kableExtra",
# Specify project-wide packages to install
neededPackages <- c("data.table", "plyr", "tidyr", "dplyr",  "Hmisc", "sjmisc", "stringr",
                    "readstata13", "foreign", "readxl", "writexl",
                    "httr", "haven",
                    "raster", "sp", "spdep", "sf","gfcanalysis", "osrm", "osrmr", "nngeo", # "rgdal"
                    "lubridate","exactextractr",
                    "doParallel", "foreach", "snow", 
                    "knitr", "kableExtra",
                    "DataCombine", 
                    "multcomp", "car", "fixest", "sandwich", "broom",  "urca", # "glmhdfe",
                    "ggeffects",
                    "ggplot2", "leaflet")

# Note that from ?renv::restore 
# "Any required recursive dependencies of the requested packages will be restored as well."

# This is just to provide the user with some information on what renv does
renv::status() 

# This installs, if necessary, in renv/library, the project-specific packages as recorded in renv.lockpackages = neededPackages
renv::restore()

# To install rgdal, see https://stackoverflow.com/questions/78746226/how-to-install-rgdal-in-r-version-4-4-1

#### /!\ IF renv::restore() FAILS TO INSTALL SOME PACKAGES FROM neededPackages /!\ #### 

# For instance sf could cause trouble https://github.com/r-spatial/sf/issues/921 
# or magick, as a dependency of raster and rgdal. 

# FOLLOW THESE STEPS:

# 1. Remove these package names from neededPackages above, and rerun renv::restore(packages = neededPackages)
# 2. Write them in troublePackages below, uncomment, and run the following code chunk: 

# # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("")
# # Attempt to load packages from user's default libraries. 
# # object_default libraries was saved in .Rprofile before renv was activated
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


### LOAD ALL PROJECT PACKAGES ### 
lapply(neededPackages, library, character.only = TRUE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### DISCLOSE PROJECT DEPENDENCIES
# The following code chunk's purpose was uniquely to disclose the packages used in the project.
# This was necessary for renv::dependencies to see the packages. 
# As a result of a one-time execution of renv::init(), the renv.lock file provided 
# with this project records all these packages in the versions installed at that time on the 
# local machine from where I executed the final, reproducible, data processing. 

if(FALSE){
  # NOT RUN
  library(data.table)
  library(plyr)
  library(tidyr)
  library(dplyr)
  library(Hmisc)
  library(sjmisc)
  library(stringr)
  library(readstata13)
  library(foreign)
  library(readxl)
  library(writexl)
  library(raster)
  library(sp)
  library(spdep)
  library(sf)
  library(gfcanalysis)
  library(osrm)
  library(osrmr)
  library(nngeo)
  library(doParallel)
  library(foreach)
  library(snow)
  library(kableExtra)
  library(knitr)
  library(leaflet)
  library(htmltools)
  #library(velox)
  #library(ExPanDaR)
  library(DataCombine)
  library(car)
  library(multcomp)
  library(fixest)
  library(urca)
  #library(effects)
  #library(glmhdfe)
  #library(sandwich)
  #library(lmtest)
  #library(pglm)
  #library(multiwayvcov)
  library(clusterSEs)
  #library(clubSandwich)
  #library(boot)
  #library(Countr)
  library(ggplot2)
  #library(doBy)
  library(lubridate)
  library(httr)
  library(haven)
  library(exactextractr)
  library(ggeffects)
  #library(emmeans)
  library(dreamerr)
  library(stargazer)
}

# necessary currently because renv changes the downloaded method at every new session to curl... 
if(renv:::renv_download_method() != getOption("download.file.method")){
  Sys.setenv(RENV_DOWNLOAD_METHOD = getOption("download.file.method"))
}

# If new packages are needed along the project data processing workflow, one should: 
new_pck <- c("stargazer")
# 1. install the packages in the project library (the default if you are within the project)
install.packages(new_pck)
# devtools::install_github("julianhinz/R_glmhdfe")
# 2. add 'library(package)' in the list above, and ',"package"' in the neededPackages vector above, and save the present file.
# 3. perform an implicit (the default) snapshot, that writes to the renv.lock file the packages at
# the intersection between packages found in the project by renv::dependencies() (hence the 'library(package)' line from 2. above.) 
# and packages installed the specified library (the default being the project's) (hence point 1.)
renv::snapshot() 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


#### DEPRECATED 
# Deprecated way to deal with renv::restore() failure for some packages. 

### SECOND-BEST SOLUTION for packages for which installation through renv::restore() fails
  # For instance sf https://github.com/r-spatial/sf/issues/921 
  # or magick, as a dependency of raster and rgdal. 

# /!\ Note that these packages won't necessary be installed in the same version as they were on the project initial data processing
#     and therefore reproducibility is not guaranteed in this case. 

### Follow the steps: 
# 1. In the present script, above, comment out the lines disclosing (by library(package)) the packages 
#   causing renv::restore() to fail; and save the changes made to the script. 
# 2. Run command of line 27 of the present script: renv::snapshot() to remove these packages from the renv.lock file. 
# 3. Run command of line 28 of the present script: renv::restore() to check if all remaining packages are installed properly. 
# 3. Report all the packages you had to remove from the renv.lock for renv::restore() not to fail, 
# #   in this trouble_packages vector 
#     trouble_packages <- c("")
# #   and install them if not available
#     allPackages    = c(trouble_packages %in% installed.packages()[ , "Package"]) 
#     if(!all(allPackages)) {
#       missingIDX = which(allPackages == FALSE)
#       needed     = trouble_packages[missingIDX]
#       lapply(needed, install.packages)
#     }
    


