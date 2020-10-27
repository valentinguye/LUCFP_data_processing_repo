# Save the default user library paths before renv is activated
default_libraries <- .libPaths()

# # Install renv if not available (activate.R fails to do that if package is not already installed)
# if(!("renv" %in% utils::installed.packages()[,"Package"])){
#   if (!requireNamespace("remotes")){
#     install.packages("remotes")
#   }
# 
#   remotes::install_github("rstudio/renv")
#   print("renv was not already installed. It has been installed in default library")
# } else print("User already has renv installed in default library")
# 
# # We run this before source("renv/activate.R") because otherwise error 35 is thrown
# # if renv is not available in user default library. 
# renv::activate()
# # 
# # # Activate renv
source("renv/activate.R")

# renv::activate()

# Inform about this .Rprofile completion
print("The above outputs come from the execution of LUCFP/data_processing/.Rprofile at R start up")
