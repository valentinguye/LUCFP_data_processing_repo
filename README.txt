IF USER HAS Stata/MP: 
The complete data work can be reproduced running master.do with Stata/MP > 14.0
Notes on outside-Stata programs: 
	- GEE: Google Earth Engine scripts are not called from master.do ; they can be provided upon request (they extract descriptive statistics). 
	- R: master.do calls R scripts; it transfers them its current working directory, and each R script invokes paths relative to this. 


IF USER DOES NOT HAVE Stata/MP and wants to run only R scripts. 
*** The R_project_for_individual_runs.Rproj file should be opened first *** 
(in particular, for each R script to run in the correct working directory, and with the correct packages and packages' versions). 
Then R scripts can be sourced within this project environment;
There is no master R script because it would make no sense to run all R scripts without running inserted Stata scripts.
How Stata and R scripts are interrelated is described below, as a paste of master.do. 


DATA: 
The original input data are stored in the read only directory ~/LUCFP/data_processing/input_data
This is a several Gb folder, that is not available from Github, but is made available upon request. 
Some data are also downloaded within scripts (e.g. GFC tiles, from prepare_lucpfip.R) and internet connexion is therefore required for those scripts and for the complete replication run. 

The temporary data prepared within this project are stored in the write and read directory  ~/LUCFP/data_processing/temp_data
This folder is also large, hence is not available from Github either, but available upon request. 


RUN TIME 
Several scripts can take a long time to be executed, this is the case in particular for: 
prepare_lucpfip.R
prepare_lucfip.R
wa_at_parcels.R
Indicative run times should be given in the master.do where relevant. 


PACKAGES
Stata: 


R: 
At start up, the project .Rprofile calls ~renv/activate.R which, in particular, installs, loads, or updates, if necessary, the renv package. 
Then, the first thing any of this project's R scripts does is renv:restore() which installs the project required 
library (packages in specific versions) if it is not up to date in ~/renv/library/ (basically at first run on a new machine). 
These packages and their respective versions are registered in the renv.lock file provided with the project (the project local library is not). 

The packages required in this project are disclosed to renv::dependencies(), used by renv::snapshot() to update the renv.lock. Therefore:
If the user needs to load new packages, she should add a "library(package)" line in the disclose_dependencies.R script, save it, and run the
renv::snapshot() command. See ?renv::dependencies to know more. 
It the user needs to install new packages, she might have error messages of the kind 'package "..." is not available.' 
A work around is to remove all renv architecture (renv/ and renv.lock), install the new package (in default user lib), 
write a line library(new.package) in disclose_dependencies.R and run renv::init(). The new renv/library should contain the new package. 
But this would completely change the renv.lock file: packages would be written there in their versions from user default libraries, and not from project anymore. 


***********************************************************************************************************************************************
MASTER FILE PASTED BELOW TO INFORM USERS WITHOUT STATA/MP ON PROJECT DEPENDENCIES 
***********************************************************************************************************************************************





