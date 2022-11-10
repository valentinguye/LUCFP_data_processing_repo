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
neededPackages = c(neededPackages = c("tibble", "plyr", "dplyr", "data.table", "sjmisc", "stringr","Hmisc",
                                      "foreign", "readstata13", "readxl",
                                      "raster", "rgdal",  "sp", "spdep", "sf",
                                      "DataCombine",
                                      "knitr", "kableExtra",
                                      "car",  "fixest", "sandwich", "boot",
                                      "ggplot2", 
                                      "here") )

renv::restore(packages = neededPackages)

# Load them
lapply(neededPackages, library, character.only = TRUE)

here()
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


### ### ### 
# READ AND PREPARE DATA 
# read refinerie
refine <- read.csv(here("input_data", "trase_refineries.csv"))

# final IBS-UML panel 
ibs_uml <- read.dta13(here("temp_data/IBS_UML_panel_final.dta"))

# Remove IBS manufactories that are matched with Universal Mill List (UML)
# --> KEEP only those firm_id that were NOT matched with UML, be they in the IBS_UML_panel_final (==0) or not (NA)
ibs_noUML <- dplyr::filter(ibs_uml, ibs_uml$uml_matched_sample == 0) 
# sum(ibs_noUML_cs$uml_matched_sample)  # 470
# sum(ibs_noUML_cs$is_mill) # 930 (they sourced FFB at least once or sold CPO or PKO at least once, and that they are not located in Java or in Bali)

#### CONTEXTUAL VARIABLES ####  
# collapse to cross section, based on some summarized values
ibs_noUML_cs <- ddply(ibs_noUML, "firm_id", summarise, 
                      uml_matched_sample = unique(uml_matched_sample), 
                      is_mill = unique(is_mill),
                      mill_name = unique(mill_name),
                      avg_in_ton_cpo_imp1 = mean(in_ton_cpo_imp1, na.rm = TRUE), 
                      avg_out_ton_rpo_imp1 = mean(out_ton_rpo_imp1, na.rm = TRUE),
                      avg_out_ton_rpko_imp1 = mean(out_ton_rpko_imp1, na.rm = TRUE),
                      max_in_ton_cpo_imp1 = max(in_ton_cpo_imp1, na.rm = TRUE), 
                      max_out_ton_rpo_imp1 = max(out_ton_rpo_imp1, na.rm = TRUE),
                      max_out_ton_rpko_imp1 = max(out_ton_rpko_imp1, na.rm = TRUE))
# THE WARNINGS are for plants that have only NAs on some of these variables

ibs_noUML_cs <- dplyr::mutate(ibs_noUML_cs, 
                              any_rpo = is.finite(avg_out_ton_rpo_imp1) & max_out_ton_rpo_imp1>0, 
                              any_rpko = is.finite(avg_out_ton_rpko_imp1) & max_out_ton_rpko_imp1>0, 
                              any_incpo = is.finite(avg_in_ton_cpo_imp1) & max_in_ton_cpo_imp1>0)
ibs_noUML_cs <- dplyr::mutate(ibs_noUML_cs, any_refined = any_rpo | any_rpko)

ibs_noUML <- merge(ibs_noUML,ibs_noUML_cs, by = "firm_id")


#### MATCH WITH MD #### 
md <- read_excel(file.path("input_data/manufacturing_directories/direktori_industri_merged_cleaned.xlsx"))
unique(md$main_product)

sjmisc::str_contains(unique(md$main_product), "olahan", ignore.case = T)
sjmisc::str_contains(unique(md$main_product), "goreng", ignore.case = T)

md[md$main_product=="MINYAK GORENG",]
# ibs_noUML$workers_total_imp3[ibs_noUML$workers_total_imp3 %in% ]

## Match to each unref record, all the records of MD that have the same number of workers (i.e. be it the same year or not, same district or not)
#this adds a list column. There is one list element for each unref row. 
match <- dplyr::nest_join(ibs_noUML, md, by = c("workers_total_imp3" = "no_workers"), keep = TRUE, name = "y")
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
  # it's null if address is empty. In this case, we don't want to match
  if(length(match$y[[i]]$address)==0){kab_filter <- FALSE}
  # keep only the matches that are in the same district
  match$y[[i]] <- dplyr::filter(match$y[[i]], kab_filter)
  # report the number of different mills that matched
  match$n_match[i] <- length(unique(match$y[[i]]$company_name)) 
}


# 
# match$n_match2 <- NA
# f <- function(i){
#   kab_filter <- str_contains(x = match$district_name[i], pattern = match$y[[i]]$address, ignore.case = TRUE, switch = TRUE) 
#   # keep only the matches that are in the same district
#   match$y[[i]] <- filter(match$y[[i]], kab_filter)
#   # report the number of different mills (based on COMPANY NAME) that matched
#   match$n_match2[i] <- length(unique(match$y[[i]]$company_name)) 
#   return(match$n_match2[i])
# }
# match$n_match2 <- sapply(1:nrow(match), f)
# 
# 
# table(match$n_match2, match$any_cpo_output)


## Make different categories of establishments, depending on how many different matches they have over their records with MD mills. 

# make groups of establishments, based on how many matches they have repeatedly
grp_n_match <- ddply(match, "firm_id", summarise, 
                     # those establishments that never match
                     no_match = length(n_match[n_match == 0])==length(year), 
                     # one_match category allows for some records, but not for all, to have zero match. 
                     # single matches can be different from one year to another. 
                     one_match = length(n_match[(n_match == 1 | n_match == 0) & no_match == FALSE])==length(year), 
                     # svl_match is true as soon as their is as least one year with 2 different matches. 
                     svl_match = no_match == FALSE & one_match == FALSE)

# ddply note: summarise is the .fun argument and then in ... we give the ... argument of summarise, 
# i.e. we give the "name = value pairs" (see ?summarise). 
match <- merge(match, grp_n_match, by = "firm_id") 


# substract from the one_match category those that have different company name matches across years.
for(i in unique(match[match$one_match == TRUE, "firm_id"])){
  # extract the names of all the MD matches of establishment i within the one_match category 
  names <- lapply(match[match$firm_id == i, "y"], function(i.elmt) i.elmt$company_name)
  names <- unlist(names)
  # for those who have matched different company names across years, switch one_match from TRUE to FALSE 
  new_logicals <- rep(FALSE, nrow(match[match$firm_id == i,]))
  match[match$firm_id == i, "one_match"][length(unique(names)) > 1] <- new_logicals
  # for those who have matched different company names across years, switch svl_match from FALSE to TRUE
  new_logicals <- rep(TRUE, nrow(match[match$firm_id == i,]))
  match[match$firm_id == i, "svl_match"][length(unique(names)) > 1] <- new_logicals
}
# checks: 
# when the company name is the same across annual matches, the one_match status remains unchanged
match[match$firm_id == 2028, "y"]
grp_n_match[grp_n_match$firm_id == 2028,]
match[match$firm_id == 2028, c("no_match","one_match", "svl_match")]

# when the company name is not the same across annual matches, the one_match status changes from TRUE to FALSE and the svl_match from FALSE to TRUE
match[match$firm_id == 2076, "y"]
grp_n_match[grp_n_match$firm_id == 2076,]
match[match$firm_id == 2076, c("no_match","one_match", "svl_match")]

# descriptive part
describe(match[match$no_match == TRUE, "firm_id"]) # 821 establishments (2426 records)
describe(match[match$one_match == TRUE, "firm_id"]) # 118 establishments (472 records)
describe(match[match$svl_match == TRUE, "firm_id"]) # 64 establishments (356 records)

match$i_n_match <- "never matches with anything"
match$i_n_match[match$one_match == TRUE] <- "matches always with the same company name or with nothing"
match$i_n_match[match$svl_match == TRUE] <- "matches with several company names, either the same year or across years"

ddply(match, c("i_n_match","any_refined"), summarise, 
      n_mills = length(unique(firm_id)))

# la question est est-ce qu'on décide de valider systématiquement les cas où il n'y a zéro ou qu'un seul match toujours identique entre les années d'une mill ibs. 
# on pourrait dire : oui a condition qu'il y ait au moins deux occurrences de ce match. 
# ou même pas, manuellement on avait validé même quand il n'y avait qu'une obs. qui matchait.
# une partie de ces cas sont écartés ensuite pendant la phase de résolution des conflits. 

wtn_cfl <- match[match$svl_match == TRUE & match$any_refined == TRUE, ]
length(unique(wtn_cfl$firm_id))

#### WITHIN CONFLICTS ####

### Prepare the data to resolve "within" conflicts, i.e. conflicts in matched company names within each firm_id. 

# add variables on the total different matches for one ibs establishment.  
wtn_cfl$diff_names <- rep(NA, nrow(wtn_cfl))
wtn_cfl$n_diff_names <- rep(NA, nrow(wtn_cfl))
for(i in unique(wtn_cfl$firm_id)){
  names <- lapply(wtn_cfl[wtn_cfl$firm_id == i, "y"], function(i.elmt) i.elmt$company_name)
  wtn_cfl[wtn_cfl$firm_id == i, "diff_names"] <- paste(unique(unlist(names)), collapse = "; ")
  wtn_cfl[wtn_cfl$firm_id == i, "n_diff_names"] <- length(unique(unlist(names)))
}

#wtn_cfl[wtn_cfl$firm_id == 1763, "y"]


# extract all the matches from the list column 
l <- list()
for(i in 1:nrow(wtn_cfl)){
  s <- wtn_cfl[i, "y"][[1]]
  s$matched_firm_id <- rep(wtn_cfl[i, "firm_id"], nrow(wtn_cfl[i, "y"][[1]]))
  s$matched_year <- rep(wtn_cfl[i, "year"], nrow(wtn_cfl[i, "y"][[1]]))
  l[[i]]<- s
}
md_matches <- bind_rows(l)
rm(l)


# and merge them with the panel wtn_cfl
# rename first the year variable in md_matches
names(md_matches)[names(md_matches) == "year"] <- "md_year"
wtn_cfl <- merge(wtn_cfl, md_matches, by.x = c("firm_id","year"), by.y = c("matched_firm_id", "matched_year"), all = TRUE)

# export relevant variables to manual work
wtn_cfl <- dplyr::select(wtn_cfl,	firm_id, min_year, year,	workers_total_imp3,	n_diff_names, diff_names,  
                         md_year, company_name, no_workers, 
                         district_name, kec_name,	village_name, address,
                         main_product, 
                         in_ton_ffb_imp1,	in_ton_ffb_imp2, out_ton_cpo_imp1,	out_ton_cpo_imp2,	out_ton_pko_imp1, out_ton_pko_imp2,	
                         out_ton_rpo_imp1, out_ton_rpo_imp2, out_ton_rpko_imp1, out_ton_rpko_imp2,
                         pct_own_cent_gov_imp,	pct_own_loc_gov_imp,	pct_own_nat_priv_imp,	pct_own_for_imp)



#### BETWEEN CONFLICTS #### 

### now we want to resolve "between" conflicts, i.e. conflicts in matched company names between IBS firm_ids. 
# For this purpose, it is necessary to have all years for each ibs firm_id in the two categories 
# (resolved conflict cases and one_match cases (only those who produced CPO at least once))

## Within resolved conflicts
# import mannually done work 
wtn_cfl_resolved <- read_excel(file.path("input_data/manually_matched_ibs_uml/matching_unref/wtn_cfl_done.xlsx"))

# for each firm_id, keep only the row with the mannually deemed correct company name. 
# (because only one row per firm_id was given the resolved company name in the manual work in Excel)
wtn_cfl_resolved <- wtn_cfl_resolved[is.na(wtn_cfl_resolved$company_name)== FALSE,] 

length(unique(wtn_cfl_resolved$company_name)) 
# So there was 104 different firm_id that had within conflicting md company names. 
# For 102 of them, we could resolve the within conflict. 
# Among them, there are only 100 unique company names, meaning that there are some between conflicts. 

# rename the company name variable 
names(wtn_cfl_resolved)[names(wtn_cfl_resolved) == "company_name"] <- "within_resolved_c_name"

# merge it with the wtn_cfl data frame, 
wtn_cfl_resolved <- merge(wtn_cfl, wtn_cfl_resolved[, c("firm_id", "within_resolved_c_name")], by = c("firm_id"), all = TRUE)

# we don't keep only the records where the company_name is the one that was chosen mannually, because in some cases we might  
# need records for the same mill but with a  differently spelled name. 


## one_match cases
# one should just reproduce the procedure applied to prepare data for resolution of within conflicts. 

#keep only the one_match cases that are potential mills (they produced CPO at least once)
no_cfl <- match[match$one_match == TRUE & match$any_cpo_output == TRUE,]

# add variables on the total different matches for one ibs establishment.  
no_cfl$diff_names <- rep(NA, nrow(no_cfl))
no_cfl$n_diff_names <- rep(NA, nrow(no_cfl))
for(i in unique(no_cfl$firm_id)){
  names <- lapply(no_cfl[no_cfl$firm_id == i, "y"], function(i.elmt) i.elmt$company_name)
  no_cfl[no_cfl$firm_id == i, "diff_names"] <- paste(unique(unlist(names)), collapse = "; ")
  no_cfl[no_cfl$firm_id == i, "n_diff_names"] <- length(unique(unlist(names)))
}

#no_cfl[no_cfl$firm_id == 1763, "y"]

# extract all the matches from the list column 
l <- list()
for(i in 1:nrow(no_cfl)){
  s <- no_cfl[i, "y"][[1]] # there is always only one element anyways
  s$matched_firm_id <- rep(no_cfl[i, "firm_id"], nrow(no_cfl[i, "y"][[1]]))
  s$matched_year <- rep(no_cfl[i, "year"], nrow(no_cfl[i, "y"][[1]]))
  l[[i]]<- s
}
md_matches <- bind_rows(l)
rm(l)


# and merge them with the panel no_cfl
# rename first the year variable in md_matches
names(md_matches)[names(md_matches) == "year"] <- "md_year"
no_cfl <- merge(no_cfl, md_matches, by.x = c("firm_id","year"), by.y = c("matched_firm_id", "matched_year"), all = TRUE)

# export relevant variables to manual work
no_cfl <- dplyr::select(no_cfl,	firm_id, min_year, year,	workers_total_imp3,	n_diff_names, diff_names,  
                        md_year, company_name, no_workers, 
                        district_name, kec_name,	village_name, address,
                        main_product, 
                        in_ton_ffb_imp1,	in_ton_ffb_imp2, out_ton_cpo_imp1,	out_ton_cpo_imp2,	out_ton_pko_imp1, out_ton_pko_imp2,	
                        out_ton_rpo_imp1, out_ton_rpo_imp2, out_ton_rpko_imp1, out_ton_rpko_imp2,
                        pct_own_cent_gov_imp,	pct_own_loc_gov_imp,	pct_own_nat_priv_imp,	pct_own_for_imp)


## Spot the between conflicts
# have the within_resolved_c_name column in no_cfl data frame too.
no_cfl$within_resolved_c_name <- no_cfl$company_name

# merge them 
unique_match <- merge(wtn_cfl_resolved, no_cfl, all = TRUE)

# now identify duplicates across firm_id
btw_duplicates <- ddply(unique_match, c("within_resolved_c_name"), summarise, 
                        btw_duplicates = length(unique(firm_id)))

unique_match <- merge(unique_match, btw_duplicates, by = "within_resolved_c_name", all = TRUE)

describe(unique_match$btw_duplicates)
# 79 firm_id have no within resolved company name. 
all_na(unique_match[unique_match$btw_duplicates == 79,"within_resolved_c_name"])

# keep only between conflicts
btw_cfl <- unique_match[unique_match$btw_duplicates > 1 & is.na(unique_match$within_resolved_c_name) == FALSE,] 

setorder(btw_cfl, within_resolved_c_name, firm_id, year)






#### SPATIAL MATCHING #### 

## Indonesian CRS
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# IBS desa geometries 
ibs_desa <- readRDS(file.path("temp_data/processed_mill_geolocalization/IBSmills_desageom.Rdata"))
# - note that this is a cross section
nrow(ibs_desa[!duplicated(ibs_desa$firm_id),]) == nrow(ibs_desa)
# and it goes only to 2010 at the latest. The desa variable was not shipped for years after 2010.
# so missing geometries are not recent firms. 
ibs_desa <- dplyr::select(ibs_desa, firm_id, village_name, kec_name, district_name, geom) 
# import geo names from here, because they are the most recent valid ones. 
ibs_desa <- dplyr::filter(ibs_desa, !is.na(geom))

ibs_desa <- st_as_sf(ibs_desa, crs = 4326)
ibs_desa <- st_transform(ibs_desa, crs = indonesian_crs)

# append it to the main data 
ibs_noUML_cs <- left_join(ibs_noUML_cs, ibs_desa, by = "firm_id")

# split the ibs_noUML_cs data set into those with valid desa polygon, and those without (either bc they appeared after 2010, or because they have an invalid polygon)
ibs_noUML_cs_nodesa <- dplyr::filter(ibs_noUML_cs, st_is_empty(geom)) # 544 - these are those that cannot be matched with spatial help 
ibs_noUML_cs_wdesa <- dplyr::filter(ibs_noUML_cs, !st_is_empty(geom)) # 459

# make those with desa a spatial data frame (it's already in indonesian crs)
ibs_noUML_cs_wdesa <- st_as_sf(ibs_noUML_cs_wdesa, crs = indonesian_crs)
ibs_noUML_cs_wdesa

# remove the geom column from the other objects
ibs_noUML_cs_nodesa <- dplyr::select(ibs_noUML_cs_nodesa, -geom)
ibs_noUML_cs <- dplyr::select(ibs_noUML_cs, -geom)


# prepare Trase refineries dataset
refine$latitude <- as.numeric(refine$latitude)
refine$longitude <- as.numeric(refine$longitude)
refine$lat <- refine$latitude
refine$lon <- refine$longitude
refine <- st_as_sf(refine,	coords	=	c("longitude",	"latitude"), crs = 4326)
refine <- st_transform(refine, crs = indonesian_crs)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### MATCH IBS AND UML MILLS ####

ibs_j <- st_join(x = ibs_noUML_cs_wdesa, y = refine, join = st_contains, left = T) # is equivalent to ibs_j <- st_join(x = ibs_noUML_wdesa, y = refine) bc default is st_intersect.

#plot(st_geometry(ibs_noUML[ibs_noUML$firm_id == 1761,]))
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
ibs_unref <- left_join(x = ibs_noUML_cs_wdesa, y = st_set_geometry(oto[,c("firm_id","lat")], NULL), by = "firm_id")
ibs_unref <- filter(ibs_unref, is.na(ibs_unref$lat))

ibs_unref <- left_join(x = ibs_unref, y = st_set_geometry(noto[,c("firm_id","lon")], NULL), by = "firm_id")
ibs_unref <- filter(ibs_unref, is.na(ibs_unref$lon))

ibs_unref <-dplyr::select(ibs_unref, -lon, -lat)
# (the use of "lat" and "lon" was just an arbitrary choice of non-empty variables to flag oto and noto resp.)

# filter to those likely to be refineries
head(ibs_unref)
summary(ibs_unref$any_rpo)
summary(ibs_unref$any_rpko)
summary(ibs_unref$any_incpo)
nrow(dplyr::filter(ibs_unref, any_rpo | any_rpko)) # "or rpko" adds 8 plants in addition to the 70 that sell rpo. 
nrow(dplyr::filter(ibs_unref, any_rpo & any_incpo))
nrow(dplyr::filter(ibs_unref, any_rpo | any_incpo))

# here the criterion for being a refinery is to output at least some rpo OR some rpko. 


#### Automatic conflict resolution #### 
summary(refine$cap_mt) # this is expressed in metric tons, not in million tons, most likely 
noto <- dplyr::mutate(noto, excess_cap = max_out_ton_rpo_imp1 > cap_mt)
summary(noto$excess_cap)

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



# ***** Note ********
# Lines below are not a pb anymore, because it's ibs_noUML_cs on the left of the left join now (so keeping only final firm_id anyway)
# among the 17 firm_id that are manually changed in merging_geolocalization_works.do, 
# 13 have a desa geom and thus, are present in IBSmills_desageom.Rdata BUT NOT in IBS_UML_panel_final.dta  
# but we have identified them as being UML mills already (so we remove them here)
# nrow(ibs_noUML[is.na(ibs_noUML$uml_matched_sample),]) 
# ******   *******

unique(ibs_noUML_cs$mill_name)
# mill_name is actually useless (it has been wiped out when a firm_id was not UML)
# ibs_noUML <- dplyr::select(ibs_noUML, -mill_name)

