# SO, what we wanna do here contributes to the effort of georeferencing IBS mills. 
# More specifically; we want to georeference those mills that have an a priori valid desa_id
# but no referenced (UML) mill is within the associated village polygons.
# In many cases, the IBS observation is actually not a mill but rather a refinery or an other factory using some CPO in its inputs. 
# There are still cases where these observations are mills. They are then either 
# - UML mills actually laying slightly outside the polygon they reported in IBS. 
# - Not UML mills. 

# In the first case, we want to geolocalize them using the coordinates known in UML. 
# To do so, we need to associate a name to IBS unref mills. 
# This can be done thanks to the manufacturing directories. 


# First, screen (flag) only ibs_unref mills that produce at least once some CPO. 

# Second, MATCH 
# For each of them, there are several obs. of number of workers. 
# report all the lines in MD that match with one of these observations; 
# AND whose adresses encompasse the IBS mill's district. 

# Third, prepare for manual resolution of WITHIN CONFLICTS

# Fourth, prepare for manual resolution of BETWEEN CONFLICTS

# And finally put everything back together. 

#### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
#### OR CALLED FROM LUCFP PROJECT master.do FILE.
#### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###

# Installs all the packages required in this project, if not already installed in LUCFP/data_processing/renv/library
renv::restore()

# These are the packages needed in this particular script. 
neededPackages = c("data.table", "dplyr", "plyr", "sjmisc", "stringr","Hmisc",
                   "readxl", "writexl", "foreign", "readstata13", 
                   "sf", "rgdal")

# Load them
lapply(neededPackages, library, character.only = TRUE)


#### PREPARE DATA ####

# read in the digitized manufacturing directories.
md <- read_excel(file.path("input_data/manufacturing_directories/direktori_industri_merged_cleaned.xlsx"))
# read in the full IBS panel 
ibs <- read.dta13(file.path("temp_data/processed_IBS/IBS_PO_98_15_cleaned.dta"))

# get only firm ids of unref "mills". 
unref_cross <- readRDS(file.path("temp_data/processed_mill_geolocalization/ibs_unref.Rdata"))
unref_cross <- unref_cross[,"firm_id"]
unref_cross <- st_set_geometry(unref_cross, NULL)

# get every annual records of each unref mill.
unref <- merge(ibs, unref_cross, by = "firm_id")
setorder(unref, firm_id, year)

length(unique(unref[,"firm_id"]))  == nrow(unref_cross)

# flag those mills who never sell CPO
unref_dt <- data.table(unref)
cpo <- unref_dt[, mean(out_ton_cpo, na.rm = TRUE), by = firm_id]
names(cpo)[names(cpo) != "firm_id"] <- "any_cpo_output"
cpo[,"any_cpo_output"] <- !is.na(cpo[,"any_cpo_output"])

unref <- merge(unref, cpo, by = "firm_id", all = TRUE)
length(unique(unref[unref$any_cpo_output == TRUE,"firm_id"]))
# 268 different mills have produced CPO at least once. 


#### MATCH WITH MD #### 

## Match to each unref record, all the records of MD that have the same number of workers (i.e. be it the same year or not, same district or not)
#this adds a list column. There is one list element for each unref row. 
match <- nest_join(unref, md, by = c("workers_total_imp3" = "no_workers"), keep = TRUE)
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
describe(match[match$no_match == TRUE, "firm_id"]) # 371 establishments (1697 records)
describe(match[match$one_match == TRUE, "firm_id"]) # 111 establishments (924 records)
describe(match[match$svl_match == TRUE, "firm_id"]) # 110 establishments (1399 records)

match$i_n_match <- "never matches with anything"
match$i_n_match[match$one_match == TRUE] <- "matches always with the same company name or with nothing"
match$i_n_match[match$svl_match == TRUE] <- "matches with several company names, either the same year or across years"

ddply(match, c("i_n_match","any_cpo_output"), summarise, 
      n_mills = length(unique(firm_id)))

# la question est est-ce qu'on décide de valider systématiquement les cas où il n'y a zéro ou qu'un seul match toujours identique entre les années d'une mill ibs. 
# on pourrait dire : oui a condition qu'il y ait au moins deux occurrences de ce match. 
# ou même pas, manuellement on avait validé même quand il n'y avait qu'une obs. qui matchait.
# une partie de ces cas sont écartés ensuite pendant la phase de résolution des conflits. 

wtn_cfl <- match[match$svl_match == TRUE & match$any_cpo_output == TRUE, ]
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

write_xlsx(wtn_cfl, file.path("temp_data/processed_mill_geolocalization/wtn_cfl.xlsx"))
           
           
           
           

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
write_xlsx(btw_cfl, file.path("temp_data/processed_mill_geolocalization/btw_cfl.xlsx"))




#### PUT NAMED MILLS TOGETHER #### 
## Append the resolved within conflicts that did not have a between conflict, with the resolved between conflicts. 

# call them and order them
no_btw_cfl_panel <- unique_match[unique_match$btw_duplicates == 1,]
resolved_btw_cfl_panel <- read_excel(file.path("input_data/manually_matched_ibs_uml/matching_unref/btw_cfl_done.xlsx"))

class(resolved_btw_cfl_panel) <- class(no_btw_cfl_panel)

no_btw_cfl_panel <- setorder(no_btw_cfl_panel, firm_id, year)
resolved_btw_cfl_panel <- setorder(resolved_btw_cfl_panel, firm_id, year)


#keep only cases that could actually be resolved (between_resolved_c_name is a manually created variable) 
resolved_btw_cfl_panel <- resolved_btw_cfl_panel[is.na(resolved_btw_cfl_panel$between_resolved_c_name)==FALSE,]

# give them both a final company name column 
no_btw_cfl_panel$matched_resolved_company_name <- no_btw_cfl_panel$within_resolved_c_name
resolved_btw_cfl_panel$matched_resolved_company_name <- resolved_btw_cfl_panel$between_resolved_c_name

#and remove useless columns
no_btw_cfl_panel <- no_btw_cfl_panel[, c("firm_id", "year", "matched_resolved_company_name", "district_name", "kec_name", "village_name")]
resolved_btw_cfl_panel <- resolved_btw_cfl_panel[, c("firm_id", "year", "matched_resolved_company_name", "district_name", "kec_name", "village_name")]

# harmonize NAs (there are both "" and NA) in village name variable
no_btw_cfl_panel[no_btw_cfl_panel$village_name == "", "village_name"] <- NA
#resolved_btw_cfl_panel[is.na(resolved_btw_cfl_panel$village_name), "village_name"] <- ""
#resolved_btw_cfl_panel[resolved_btw_cfl_panel$village_name == "", "village_name"] <- NA


# keep for each data frame only the records that have the desa information 
no_btw_cfl_panel1 <- no_btw_cfl_panel[is.na(no_btw_cfl_panel$village_name) == FALSE 
                         & !str_contains(x = "*", pattern = no_btw_cfl_panel$village_name, switch = TRUE)
                         & no_btw_cfl_panel$district_name != "123456789xyz"
                         ,]
length(unique(no_btw_cfl_panel1$firm_id))

resolved_btw_cfl_panel1 <- resolved_btw_cfl_panel[is.na(resolved_btw_cfl_panel$village_name) == FALSE 
                         & !str_contains(x = "*", pattern = resolved_btw_cfl_panel$village_name, switch = TRUE)
                         & resolved_btw_cfl_panel$district_name != "123456789xyz"
                         ,]
length(unique(resolved_btw_cfl_panel1$firm_id))

# only the no_btw_cfl data frame has loss of mills (9) with the constraints on geographic info. 

# keep for each data frame only one instance for each ibs establishment (the earliest one here)
no_btw_cfl <- no_btw_cfl_panel[!duplicated(no_btw_cfl_panel$firm_id),]
resolved_btw_cfl <- resolved_btw_cfl_panel[!duplicated(resolved_btw_cfl_panel$firm_id),] 

no_btw_cfl1 <- no_btw_cfl_panel1[!duplicated(no_btw_cfl_panel1$firm_id),]
resolved_btw_cfl1 <- resolved_btw_cfl_panel1[!duplicated(resolved_btw_cfl_panel1$firm_id),] 

# and append those with a correct geographic information set
named_unref <- merge(no_btw_cfl1, resolved_btw_cfl1, all = TRUE)

## add the 9 mills with incomplete geographic information set (all in no_btw_cfl)

# mills for which something was missing
firm_id <- no_btw_cfl[!is.element(no_btw_cfl$firm_id, no_btw_cfl1$firm_id), "firm_id"]

# select those mills that only have a star desa name 
# (stars signal village names that are less trustworthy due to some checks in cleaning_IBS.do)
no_btw_cfl_panel <- no_btw_cfl_panel[no_btw_cfl_panel$firm_id %in% firm_id,]
# and add a record with their star name to the named_unref
no_btw_cfl_star <-no_btw_cfl_panel[is.na(no_btw_cfl_panel$village_name) == FALSE,]
no_btw_cfl_star <-no_btw_cfl_star[!duplicated(no_btw_cfl_star$firm_id),]
named_unref <- merge(named_unref, no_btw_cfl_star, all = TRUE)

# select those mills that never have a desa name
no_name <- ddply(no_btw_cfl_panel, "firm_id", summarise, no_name = sum(is.na(village_name))==length(village_name))
no_btw_cfl_name <- merge(no_btw_cfl_panel, no_name, by = "firm_id", all = TRUE)
no_btw_cfl_name <-no_btw_cfl_name[no_btw_cfl_name$no_name == TRUE,]
no_btw_cfl_name <- dplyr::select(no_btw_cfl_name, -no_name)
# and add any (the earliest) record from each of them 
no_btw_cfl_name <-no_btw_cfl_name[!duplicated(no_btw_cfl_name$firm_id),]
named_unref <- merge(named_unref, no_btw_cfl_name, all = TRUE)


named_unref <- setorder(named_unref, firm_id)

write_xlsx(named_unref, file.path("temp_data/processed_mill_geolocalization/named_unref.xlsx"))

length(unique(named_unref$firm_id))
length(unique(named_unref$matched_resolved_company_name))

