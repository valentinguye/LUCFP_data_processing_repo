# The idea of this program is to double check automatic spatial IBS-UML matches. 
# Such matches were labelled "one-to-one", abbreviated oto. 
# 
# The version of IBS_UML_panel_final.dta is an old one, that was typically the final version 
# before this double check occured, and corrections were made in merging_geolocalization_works.do
# Therefore, this script is not reproduceable from the project folder. 
# It is just a program that was used to look at some data, but without changing anything. 
# The data modifications implied are made in merging_geolocalization_works.do. 


# input:  IBS_UML_panel_final.dta (old version)
#         direktori_industri_merged_cleaned.xlsx

# No output. 


#### WORKING DIRECTORY SHOULD BE LUCFP PROJECT DATA WORK TOP LEVEL (~/LUCFP/data_processing) ####
# if this script is open within R_project_for_individual_runs


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. 
neededPackages = c("data.table","dplyr","sjmisc", "stringr","Hmisc",
                   "readxl",  "foreign", "readstata13",
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

#### PREPARE DATA ####

# read in the digitized manufacturing directories.
md <- read_excel(file.path("input_data/manufacturing_directories/direktori_industri_merged_cleaned.xlsx"))
# read in the full IBS panel, in its latest version (i.e. post all manual conflict resolutions)
ibs <- read.dta13(file.path("C:/Users/GUYE/Desktop/opalval/build/output/IBS_UML_panel_final.dta"))

# get only firm ids of oto "mills". 
oto_cross <- read.dta13(file.path("temp_data/processed_mill_geolocalization/oto.dta"))
oto_cross <- st_as_sf(oto_cross, coords = c("lon", "lat"), crs = 4326)
oto_cross <- oto_cross[,"firm_id"]
oto_cross <- st_set_geometry(oto_cross, NULL)

# get every annual records of each oto mill.
oto <- merge(ibs, oto_cross, by = "firm_id")
setorder(oto, firm_id, year)

length(unique(oto[,"firm_id"]))  == nrow(oto_cross)
oto_cross$firm_id[!is.element(oto_cross$firm_id, base::intersect(oto$firm_id, oto_cross$firm_id))]
# These five mills were oto, and were found to be the same mill has another one, with a lower firm_id 
# (i.e. established earlier) therefore they are not in IBS_UML_panel_final.dta anymore and not in oto. 

oto[as.numeric(oto$merge_oto)==1,c("firm_id", "year", "merge_noto", "merge_unref_geo")]
# these guys were not noto, initially, and had a different firm_id, that was changed to the one currently
# displayed, as they were found to be the same mill as those 3 oto mills. 

# flag those mills who never sell CPO
oto_dt <- data.table(oto)
cpo <- oto_dt[, mean(out_ton_cpo, na.rm = TRUE), by = firm_id]
names(cpo)[names(cpo) != "firm_id"] <- "any_cpo_output"
cpo[,"any_cpo_output"] <- !is.na(cpo[,"any_cpo_output"])

oto <- merge(oto, cpo, by = "firm_id", all = TRUE)
length(unique(oto[oto$any_cpo_output == TRUE,"firm_id"]))
# 2 of them have never produced CPO. 
oto[oto$any_cpo_output == FALSE,"firm_id"]
# this is 2036 (5 years) and 55630. There coordinates point to actual mills, so we don't remove them from IBS-UML


#### MATCH WITH MD #### 

## Match to each oto record, all the records of MD that have the same number of workers (i.e. be it the same year or not, same district or not)
#this adds a list column. There is one list element for each oto row. 
match <- nest_join(oto, md, by = c("workers_total_imp3" = "no_workers"), keep = TRUE)
# each list element is a dataframe
class(match$y[[1]])

## Restrict to matches that are in the same district
#some cleaning of the district name variable (no need to bother about case, bc its handled in the function)
match$district_name <- str_replace(string = match$district_name, pattern = "Kab. ", replacement = "")
# not elegant but enables that empty district_names are not matched with adresses in MD
match$district_name[match$district_name == ""] <- "123456789xyz"


match$n_match <- NA
for(i in 1:nrow(match)){
  #specified as such (with switch = T), the function checks whether x is in any of the patterns (wierd phrasing but that's the way to go with this function)
  kab_filter <- str_contains(x = match$district_name[i], pattern = match$y[[i]]$address, ignore.case = TRUE, switch = TRUE) 
  # keep only the matches that are in the same district
  match$y[[i]] <- filter(match$y[[i]], kab_filter)
  # report the number of different mills that matched
  match$n_match[i] <- length(unique(match$y[[i]]$company_name)) 
}


# In oto, and hence in match, the mill_name and parent_co variables come from UML 
# (through automatic spatial matching)
# we check here whether these names correspond to manufacturing directory company names of establishment
# matching at least one year in terms of number of workers and district. 

match$mill_name_match <- NA
for(id in unique(match$firm_id)){
  
  all_matched_md <- bind_rows(match[match$firm_id == id,"y"])
  matched_md_names <- unique(all_matched_md$company_name)
  
  mill_name_parts <- match[match$firm_id == id, "mill_name"] %>% unique() %>%  strsplit(" ") %>% unlist()
  parent_co_name_parts <- match[match$firm_id == id, "parent_co"] %>% unique() %>%  strsplit(" ") %>% unlist()
  uml_name_parts <- c(mill_name_parts, parent_co_name_parts) %>% as.list()
  
  uml_name_parts_matching <- lapply(uml_name_parts, FUN = function(x){str_contains(matched_md_names, x)}) %>% unlist() %>% sum()
  
  match[match$firm_id == id, "mill_name_match"] <- uml_name_parts_matching > 0
}

oto_to_check <- match[match$mill_name_match == FALSE,]
# these guys need to be checked deeper (manually)
oto_to_check <- dplyr::select(oto_to_check, y, firm_id, year, workers_total_imp3, trase_code, mill_name, parent_co, 
                              district_name, kec_name, village_name, lat, lon)

oto_to_check$firm_id %>% unique()
# So these 16 cases are conflicts between spatial matching (SM) and 
# manufacturing directory matching (MD)

# 50453 and 54271 do not have a mill name in this version, which means that other IBS firms
# were found to be the mills they were spatially matched with. 
# may be they can be matched with mills from MD?

# 2042
# https://www.google.com/maps/place/PT.+Varem+Sawit+Cemerlang+(VSC)/@2.6631912,99.6651437,844m/data=!3m2!1e3!4b1!4m10!1m4!3m3!1s0x0:0x0!2zMsKwMzknNTQuMCJOIDk5wrA0MCcwMS4yIkU!3b1!3m4!1s0x3032696d0a5bf5d3:0xd49167b3c06c5425!8m2!3d2.6631912!4d99.6673324
# seems to confirm SM


firm_idS <- unique(oto_to_check$firm_id)
all_matched_md <- list()
for(id in firm_idS){
 all_matched_md[[match(id,firm_idS)]] <- bind_rows(oto_to_check[oto_to_check$firm_id==id,"y"])
 all_matched_md[[match(id,firm_idS)]] <- all_matched_md[[match(id,firm_idS)]] %>% as.data.frame()
 if(nrow(all_matched_md[[match(id,firm_idS)]])>0){
   all_matched_md[[match(id,firm_idS)]]$firm_id <- id
   #all_matched_md[[match(id,firm_idS)]]$uml_mill_name <- oto_to_check[oto_to_check$firm_id == id, "mill_name"] %>% unique()
   #all_matched_md[[match(id,firm_idS)]]$uml_parent_co <- oto_to_check[oto_to_check$firm_id == id, "parent_co"] %>% unique()
   all_matched_md[[match(id,firm_idS)]] <- dplyr::select(all_matched_md[[match(id,firm_idS)]], firm_id, everything())
 }
}


all_matched_md_df <- bind_rows(all_matched_md)
oto_to_check <- dplyr::select(oto_to_check, -y)
oto_to_check <- nest_join(oto_to_check, all_matched_md_df, by="firm_id", keep = TRUE)
oto_to_check <- dplyr::select(oto_to_check, y, everything())


# So these 16 cases are conflicts between spatial matching (SM) and manufacturing directory matching (MD)
# 50453 and 54271 do not have a mill name in this version, which means that other IBS firms
# were found to be the mills they were spatially matched with. 
# may be they can be matched with mills from MD?



