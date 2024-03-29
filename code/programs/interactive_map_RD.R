
##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table",
                   "foreign", "readstata13", "readxl",
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", "kableExtra",
                   "fixest", "sandwich", "lmtest", "boot", 
                   "ggplot2", "leaflet", "htmltools")
# "pglm", "multiwayvcov", "clusterSEs", "alpaca", "clubSandwich",

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
# troublePackages <- c() 
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 


# # # /!\ THIS BREAKS THE PROJECT REPRODUCIBILITY GUARANTY /!\
# troublePackages <- c("leaflet", "leaflet.providers", "png")
# # Attempt to load packages from user's default libraries.
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing

### NEW FOLDERS USED IN THIS SCRIPT 


### INDONESIAN CRS 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE
parcel_size <- 3000


#### Get the same data as those used in regressions ####

# the code is not quite identical because we want to leave some more variables, but it's important 
# that it does not change the actual analysis sample. 

d_50 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                 parcel_size/1000,"km_",
                                 "50CR.rds")))
# Split them into islands of interest
d_50_suma <- d_50[d_50$island == "Sumatra",]
d_50_kali <- d_50[d_50$island == "Kalimantan",]
#d_50_papu <- d_50[d_50$island == "Papua",]


d_30 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                 parcel_size/1000,"km_",
                                 "30CR.rds")))
d_30 <- d_30[,names(d_50)]
# Split them into islands of interest
d_30_suma <- d_30[d_30$island == "Sumatra",]
d_30_kali <- d_30[d_30$island == "Kalimantan",]
#d_30_papu <- d_30[d_30$island == "Papua",]

rm(d_30, d_50)

catchment = "CR"
outcome_variable = "lucpfip_pixelcount"
island = "Sumatra"
alt_cr = FALSE
only_sr = FALSE
commo = c("ffb", "cpo")
x_pya = 3
dynamics = TRUE
log_prices = FALSE
yoyg = FALSE
short_run = "full"
imp = 1
distribution = "quasipoisson"
fe = "lonlat + district_year"
remaining_forest = FALSE
offset = FALSE
lag_or_not = ""
controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml")#, "wa_prex_cpo_imp1""wa_pct_own_loc_gov_imp",
baseline_forest_trend = FALSE
interaction_terms = NULL #c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1")
interact_regressors = FALSE
interacted = "regressors"
spatial_control = FALSE
pya_ov = FALSE
illegal = "all"# "no_ill2"
weights = FALSE

### SPECIFICATIONS  
### BASE DATA 

if(catchment == "CR"){
  # these are preloaded, because CR is the preferred data set
  # Catchment radius
  if(island == "Sumatra" & alt_cr == FALSE){
    d <- d_30_suma
  }
  if(island == "Sumatra" & alt_cr == TRUE){
    d <- d_50_suma
  }
  if(island == "Kalimantan" & alt_cr == FALSE){
    d <- d_50_kali
  }
  if(island == "Kalimantan" & alt_cr == TRUE){
    d <- d_30_kali
  }
  if(island == "both" & alt_cr == FALSE){
    d <- rbind(d_30_suma, d_50_kali)#, d_50_papu
  }
  if(island == "both" & alt_cr == TRUE){
    d <- rbind(d_50_suma, d_30_kali)#, d_30_papu
  }
} else if(catchment == "CA"){
  # here we load it, it slows the function but we will use it only once or twice.
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                parcel_size/1000,"km_",
                                "2h_CA.rds")))
  if(island == "Sumatra"){
    d <- d[d$island == "Sumatra",]
  }
  if(island == "Kalimantan"){
    d <- d[d$island == "Kalimantan",]
  }
  
}


### SPECIFICATIONS  

if(offset & grepl("lucpf", outcome_variable)){offset_fml <- ~log(remain_pf_pixelcount)}
#if(offset & grepl("lucf", outcome_variable)){offset_fml <- ~log(remain_f30th_pixelcount)}

## VARIABLES OF INTEREST (called regressors here)
if(dynamics == FALSE){ 
  # Only the overall effect is investigated then, no short vs. medium run
  # In this case, the variables have name element _Xya_ with X the number of years over which the mean has been 
  # computed, always including the contemporaneous record. See add_parcel_variables.R
  #short_run <- ""
  
  # if we omit one commodity 
  if(length(commo) == 1){
    if(yoyg == TRUE){
      regressors <- paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not)
    }else{
      regressors <- paste0("wa_",commo,"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)
    }
  }
  # if we don't omit a commodity. 
  if(length(commo) == 2){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not))
      
    }else{
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)) 
    }
  }
}

# if, on the other hand, we want to identify the short run effects, controlling for the medium run's. 
if(dynamics == TRUE){ 
  
  if(length(commo) == 1){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                      paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
    }
    if(yoyg == FALSE){
      if(short_run == "full"){
        regressors <- c(paste0("wa_",commo,"_price_imp",imp,lag_or_not), # SR measure
                        paste0("wa_",commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
      }
      if(short_run == "dev"){
        regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not), # SR measure
                        paste0("wa_",commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
      }  
    }
  }
  
  if(length(commo) == 2){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                      paste0("wa_",commo[1],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not))    
    }
    if(yoyg == FALSE){
      if(short_run == "full"){
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,lag_or_not),# FFB SR measure
                        paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not), # FFB LR measure 
                        paste0("wa_",commo[2],"_price_imp",imp,lag_or_not),# CPO SR measure
                        paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # CPO LR measure 
      }
      if(short_run == "dev"){
        regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                        paste0("wa_",commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not))
      }
    }
  }
}

if(only_sr){
  if(length(commo) == 1){
    if(yoyg == TRUE){
      regressors <- paste0("wa_",commo,"_price_imp",imp,"_yoyg",lag_or_not)
    }else{
      regressors <- paste0("wa_",commo,"_price_imp",imp,lag_or_not)
    }
  }
  # if we don't omit a commodity. 
  if(length(commo) == 2){
    if(yoyg == TRUE){
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,"_yoyg",lag_or_not))
      
    }else{
      regressors <- c(paste0("wa_",commo[1],"_price_imp",imp,lag_or_not),
                      paste0("wa_",commo[2],"_price_imp",imp,lag_or_not)) 
    }
  }
}

# Logarithms
if(log_prices){
  for(reg in regressors){
    d[,paste0("ln_",reg)] <- log(d[,reg])
  }
  regressors <- paste0("ln_", regressors)
  rm(reg)
}

## CONTROLS
# add pya outcome variable 
if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}

# lag controls or not
if(lag_or_not=="_lag1" & length(controls)>0){controls <- paste0(controls,lag_or_not)}

# add spatial control or not
if(spatial_control){controls <- c(controls, "ngb_ov_lag4")}

# add remaining forest or not
if(remaining_forest){
  if(grepl("lucpf", outcome_variable)){
    controls <- c(controls, "remain_pf_pixelcount")
  }
  if(grepl("lucf", outcome_variable)){
    controls <- c(controls, "remain_f30th_pixelcount")
  }
  
  offset <- FALSE
}

if(baseline_forest_trend){
  if(grepl("lucpf", outcome_variable)){
    d[,"baseline_forest_trend"] <- d$pfc2000_total_pixelcount*(d$year-2000)
    controls <- c(controls, "baseline_forest_trend")
  }
  if(grepl("lucf", outcome_variable)){
    d[,"baseline_forest_trend"] <- d$fc2000_30th_pixelcount*(d$year-2000)
    controls <- c(controls, "baseline_forest_trend")
  }
}

# ### WEIGHTS
# if(weights){
#   d$sample_coverage <- d$n_reachable_ibs/d$n_reachable_uml
# }

### SELECT DATA FOR REGRESSION

## group all the variables necessary in the regression
# important to do that after outcome_variable, regressors controls etc. have been (re)defined. 
# (interactions do not need to be in there as they are fully built from the used_vars)
used_vars <- c(outcome_variable, regressors, controls,
               "lonlat", "year", "lat", "lon", "idncrs_lat", "idncrs_lon", "district", "province", "island", "district_year", "province_year", "lucpfsmp_pixelcount")
#"n_reachable_ibsuml_lag1", "sample_coverage_lag1", #"pfc2000_total_ha", 
#"remain_f30th_pixelcount","remain_pf_pixelcount"

# if we use an offset (for now only remaining forest possible) then it should be in the necessary variables, but not among the controls. 
if(offset){
  if(grepl("lucpf", outcome_variable)){
    used_vars <- c(used_vars, "remain_pf_pixelcount")
  }
  if(grepl("lucf", outcome_variable)){
    used_vars <- c(used_vars, "remain_f30th_pixelcount")
  }
} 

## Keep observations that: 

# - are covered with some of the studied forest 
# (so remove obs. from all years of grid cells that are not forested in 2000, or years after grid cells get completely deforested)
if(grepl("lucpf", outcome_variable)){
  d <- d[d$remain_pf_pixelcount > 0,]
}

# remain_f30th_pixelcount has not been computed yet  
# if(grepl("lucf", outcome_variable)){
#   d <- d[d$remain_f30th_pixelcount > 0,]
# }

# - are in years when the outcome can actually be observed
if(grepl("_slow_",outcome_variable)){
  d <- d[d$year<2011,]
}

if(grepl("fsmp_", outcome_variable) | grepl("fsp_", outcome_variable) | grepl("fmp_", outcome_variable)){
  d <- d[d$year<2015,]
}

# should be applied to d_clean rather no ? or at least d_nona...
# to speed up computation, and also because we do not want to control for values in neighboring cells that are not entering the regressions
# these values are not bringing any bias by definition as they are excluded
if(spatial_control){
  
  # restrict the data set to years where the paste outcome is available
  d <- d[d$year > 2004,]
  
  # d is a balanced panel so far so we can take any first record of a parcel to 
  # keep the most general cross section
  nrow(d) == length(unique(d$lonlat))*length(unique(d$year))
  d_cs <- d[!duplicated(d$lonlat),c("lonlat", "year", "lat", "lon", "island",outcome_variable)]
  
  # spatial
  d_cs <- st_as_sf(d_cs, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
  
  # identify neighbors
  # this definition of neighbors includes the 8 closest, surounding, grid cells. 
  d_buf <- st_buffer(d_cs, dist = parcel_size - 10)
  row.names(d_buf) <- d_buf$lonlat
  sgbp <- st_intersects(d_buf)
  
  
  # iterate on parcel ~30 minutes
  d[,"ngb_ov_lag4"] <- rep(NA, nrow(d))
  for(p_i in 1:length(sgbp)){
    # remove own index
    p_i_neighbors <- sgbp[[p_i]][sgbp[[p_i]] != p_i]
    # get the lonlat corresponding to the sgbp index
    p_i_neighbors <- as.numeric(attr(sgbp,"region.id")[p_i_neighbors])
    # compute mean of outcome variable lags in the PANEL dataset
    for(y in unique(d$year)){
      d[d$lonlat == as.numeric(attr(sgbp,"region.id")[p_i]) & d$year == y, "ngb_ov_lag4"]  <- mean(d[d$lonlat %in% p_i_neighbors & d$year == y, paste0(outcome_variable,"_lag4")],na.rm = TRUE)
    }
  } 
  
  # empty neighbor sets (those that have no neighbors but themselves) and  end up with NaN as mean(integer(0)) 
  # turn them into NA
  d[is.nan(d$ngb_ov_lag3),"ngb_ov_lag3"] <- NA
  d[is.nan(d$ngb_ov_lag4),"ngb_ov_lag4"] <- NA
  
}


# - are not in an intended RSPO-certified supply base
d <- d[d$rspo_cert==FALSE,]

# are legal, illegal, or neither
if(illegal == "no_ill1"){
  d <- d[d$illegal1 == FALSE, ]
}
if(illegal == "no_ill2"){
  d <- d[d$illegal2 == FALSE, ]
}
if(illegal == "ill1"){
  d <- d[d$illegal1 == TRUE, ]
}  
if(illegal == "ill2"){
  d <- d[d$illegal2 == TRUE, ]
}

# - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
names(usable) <- used_vars            
usable <- bind_cols(usable)
filter_vec <- base::rowSums(usable)
filter_vec <- filter_vec == length(used_vars)
d_nona <- d[filter_vec, c(used_vars)]
if(anyNA(d_nona)){stop()}
rm(filter_vec, usable)

# - sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
# probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
# just remove them if any, so that there is no bug. 
if(weights){
  d_nona <- d_nona[d_nona$sample_coverage!=0,]
  d_nona[d_nona$sample_coverage>1,"sample_coverage"] <- 1
}

pixel_area <- (27.8*27.6)/(1e4)
d_nona$lucpfip_ha_total <- d_nona$lucpfip_pixelcount*pixel_area
d_nona$lucpfsmp_ha_total <- d_nona$lucpfsmp_pixelcount*pixel_area

# - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
d_clean_ind <- d_nona[-obs2remove(fml = as.formula(paste0("lucpfip_pixelcount ~ ", fe)),
                                  d_nona, 
                                  family = "poisson"),]
d_clean_sm <- d_nona[-obs2remove(fml = as.formula(paste0("lucpfsmp_pixelcount ~ ", fe)),
                                  d_nona, 
                                  family = "poisson"),]

d_clean_ls <- list(d_clean_ind, d_clean_sm)
# 

prepare_d_clean <- function(d_clean){
  
  parcel_avg <- ddply(d_clean, "lonlat", summarise,
                      avg_lucpfip_ha_total = mean(lucpfip_ha_total),
                      avg_lucpfsmp_ha_total = mean(lucpfsmp_ha_total),
                      # avg_wa_ffb_price_imp1_4ya = mean(wa_ffb_price_imp1_4ya),
                      # avg_wa_cpo_price_imp1_4ya = mean(wa_cpo_price_imp1_4ya),
                      avg_wa_ffb_price_imp1 = mean(wa_ffb_price_imp1),
                      avg_wa_cpo_price_imp1 = mean(wa_cpo_price_imp1)
                      )
  
  d_clean <- merge(d_clean, parcel_avg, by = "lonlat")
  # rather displaying parcel averages over time, let's display a particular year (say 2008)
  # or even, the mean deviation
  d_clean$md_lucpfip_ha_total <- (d_clean$avg_lucpfip_ha_total - d_clean$lucpfip_ha_total) %>% round(1)
  d_clean$md_lucpfsmp_ha_total <- (d_clean$avg_lucpfsmp_ha_total - d_clean$lucpfsmp_ha_total) %>% round(1)
  # d_clean$md_wa_ffb_price_imp1_4ya <- (d_clean$avg_wa_ffb_price_imp1_4ya - d_clean$wa_ffb_price_imp1_4ya) %>% round(1)
  # d_clean$md_wa_cpo_price_imp1_4ya <- (d_clean$avg_wa_cpo_price_imp1_4ya - d_clean$wa_cpo_price_imp1_4ya) %>% round(1)
  d_clean$md_wa_ffb_price_imp1 <- (d_clean$avg_wa_ffb_price_imp1 - d_clean$wa_ffb_price_imp1) %>% round(1)
  d_clean$md_wa_cpo_price_imp1 <- (d_clean$avg_wa_cpo_price_imp1 - d_clean$wa_cpo_price_imp1) %>% round(1)
  
  # d_clean_cs <- d_clean[!duplicated(d_clean$lonlat),]
  
  # industrial or smallholders? 
  # d_clean[d_clean$lucpfip_ha_total>0,] %>% nrow()
  # d_clean[d_clean$lucpfsmp_ha_total>0,] %>% nrow()
  # d_clean[d_clean$lucpfip_ha_total==0 & d_clean$lucpfsmp_ha_total>0,] %>% nrow()
  # d_clean[d_clean$lucpfip_ha_total>d_clean$lucpfsmp_ha_total & d_clean$lucpfsmp_ha_total >0 ,] %>% nrow()
  # d_clean[d_clean$avg_lucpfip_ha_total>0,] %>% nrow()
  
  # prepare shapes of a cross section 
  d_clean_cs <- d_clean[d_clean$year == 2008,]
  d_clean_cs <- st_as_sf(d_clean_cs, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs)
  d_clean_cs <- st_buffer(d_clean_cs, dist = 1500)
  st_geometry(d_clean_cs) <- sapply(st_geometry(d_clean_cs), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
  d_clean_cs <- st_transform(d_clean_cs, crs = 4326)
}

d_clean_cs_ls <- lapply(d_clean_ls, FUN = prepare_d_clean)
save_d_clean_cs_ls <- d_clean_cs_ls

# both <- d_clean_cs_ls[[1]][d_clean_cs_ls[[1]]$lonlat%in%d_clean_cs_ls[[2]]$lonlat,]# in 2008, 189 parcels are in both groups (have a positive lucpfIp and lucpfSMp over all period)
# both[both$avg_lucpfip_ha_total==both$avg_lucpfsmp_ha_total,] %>% nrow() # only 3 have the same average over the period. 
# # let's just plot d_clean_cs_ls[[1]] and d_clean_cs_ls[[2]] without respect for this overlap. 

#### DISTRICTS #### 
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

district_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/district2000_sf"))

districts_sgbp <- st_intersects(district_sf, island_sf)
# all districts are covered by only one island
districts_sgbp[lengths(districts_sgbp)>1]
# or another island than the three islands of interest. Attribute them 0
districts_sgbp[lengths(districts_sgbp)==0] <- 0
# 55 districts are in one of our three islands of interest
unlist(districts_sgbp)[unlist(districts_sgbp)!=0] %>% length()
unique(unlist(districts_sgbp))
district_sf$island <- unlist(districts_sgbp) %>% as.character()
district_sf$island <- island_sf$shape_des[match(district_sf$island, row.names(island_sf))]
district_sf_intersect <- district_sf[!is.na(district_sf$island),]

#### MILLS ####

# read the sample panel of IBS geolocalized mills
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# keep only those that are matched with UML. 
ibsuml <- ibs[ibs$uml_matched_sample==1,]
# make it a cross section
# ibsuml <- ibsuml[!duplicated(ibsuml$firm_id),]
ibsuml <- ibsuml[ibsuml$year == 2008,]
ibsuml <- ibsuml[!is.na(ibsuml$lat),]
ibsuml <- st_as_sf(ibsuml, coords = c("lon", "lat"), remove = FALSE, crs = 4326)
ibsuml <- st_transform(ibsuml, crs = indonesian_crs)
mills <- st_buffer(ibsuml, dist = 500)
ibsuml[ibsuml$island_name == "Sumatra",] <- st_buffer(ibsuml[ibsuml$island_name == "Sumatra",], dist = 3e4)
ibsuml[ibsuml$island_name != "Sumatra",] <- st_buffer(ibsuml[ibsuml$island_name != "Sumatra",], dist = 5e4)
ibsuml <- st_transform(ibsuml, crs = 4326)
mills <- st_transform(mills, crs = 4326)

#### MAP #### 

cb_ind <- colorNumeric("magma", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                       domain = st_drop_geometry(d_clean_cs_ls[[1]][,"md_lucpfip_ha_total"]),
                       #bins = 4, 
                       na.color = "transparent", 
                       reverse = F)
cb_sm <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                      domain = st_drop_geometry(d_clean_cs_ls[[2]][,"md_lucpfsmp_ha_total"]),
                      #bins = 4, 
                      na.color = "transparent", 
                      reverse = F)

# popup
mills$popup <- paste0("FFB price:", round(mills$ffb_price_imp1,0), " USD/ton", "<br/>",
                      "CPO price:", round(mills$cpo_price_imp1,0), " USD/ton", "<br/>"#,mills$lat, " ", mills$lon
                      )
d_clean_cs_ls[[1]]$popup <- paste0(#"LUCFP-industrial:", d_clean_cs_ls[[1]]$lucpfip_ha_total,
                                   "LUCFP-industrial (MD):", d_clean_cs_ls[[1]]$md_lucpfip_ha_total," ha","<br/>",
                                   #"CPO price signal:", d_clean_cs_ls[[1]]$wa_cpo_price_imp1,
                                   "CPO price signal (MD):", d_clean_cs_ls[[1]]$md_wa_cpo_price_imp1, " USD/ton")

d_clean_cs_ls[[2]]$popup <- paste0(#"LUCFP-smallholders:", d_clean_cs_ls[[2]]$lucpfsmp_ha_total,
                                   "LUCFP-smallholders (MD):", d_clean_cs_ls[[2]]$md_lucpfsmp_ha_total," ha","<br/>",
                                   #"FFB price signal:", d_clean_cs_ls[[2]]$wa_ffb_price_imp1,
                                   "FFB price signal (MD):", d_clean_cs_ls[[2]]$md_wa_ffb_price_imp1, " USD/ton")
                                   

# pas de panique pour les NAs c'est juste qu'on utilise d'autres regressors pour sélectionner le sample en vrai. 
# juste il va falloir trouver un exemple qui n'est pas NA si possible
leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  setView(lat = -0.493, 
          lng = 102.367, 
          zoom = 10)

### SHOW MILLS
leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  setView(lat = -0.493, 
          lng = 102.367, 
          zoom = 10) %>% 
  addPolygons(data = mills, 
              opacity = 1, weight = 2, fill = FALSE, col = "pink", 
              popup = ~mills$popup)

#------------------------------------------------------------------------------------#
### SHOW CATCHMENT RADIUS
leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  setView(lat = -0.493, 
          lng = 102.367, 
          zoom = 10) %>% 
  addPolygons(data = mills, 
              opacity = 1, weight = 2, fill = FALSE, col = "pink", 
              popup = ~mills$popup) %>% 
  addPolygons(data = ibsuml, 
              fill = FALSE, col = "lightblue", weight = 2)

#------------------------------------------------------------------------------------#
### SHOW GRID CELLS
leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
  setView(lat = -0.493, 
          lng = 102.367, 
          zoom = 10) %>% 
  addPolygons(data = mills, 
              opacity = 1, weight = 2, fill = FALSE, col = "pink", 
              popup = ~mills$popup) %>% 
  addPolygons(data = ibsuml, 
              fill = FALSE, col = "lightblue", weight = 2) %>% 
  addPolygons(data = d_clean_cs_ls[[2]] , 
              opacity = 0, color = "black", weight = 2, 
              fill = TRUE, fillColor = ~cb_sm(d_clean_cs_ls[[2]]$md_lucpfsmp_ha_total), fillOpacity = 0.5, 
              popup = ~d_clean_cs_ls[[2]]$popup) %>% 
  addLegend(pal = cb_sm,  
            values = d_clean_cs_ls[[2]]$md_lucpfsmp_ha_total, 
            bins = 2, opacity = 0.4,
            title = "Deforestation <br/>smallholders (2008)",
            position = "bottomright") %>% 
  
  addPolygons(data = d_clean_cs_ls[[1]], 
              opacity = 0, color = "black", weight = 2, 
              fill = TRUE, fillColor = ~cb_ind(d_clean_cs_ls[[1]]$md_lucpfip_ha_total), fillOpacity = 0.5, 
              popup = ~d_clean_cs_ls[[1]]$popup) %>% 
  addLegend(pal = cb_ind,  
            values = d_clean_cs_ls[[1]]$md_lucpfip_ha_total, 
            bins = 2, opacity = 0.7,
            title = "Deforestation <br/>industrial (2008)",
            position = "bottomright")       %>% 
  ### DISTRICTS
  addPolygons(data = district_sf_intersect, 
              fill = FALSE, col = "red",
              opacity = 0.3, weight = 2)



#--------------------------------
rowid <- m.df_wide_lonlat$lonlat[m.df_wide_lonlat$lonlat == 25000]
head(m.df$lonlat)
mills$timetop1 <- dur_mat[rowid,]/60
mills <- st_transform(mills, crs = indonesian_crs)
mills <- st_buffer(mills, dist = 50000)
mills <- st_transform(mills, crs = 4326)
cb_timetop1 <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                                domain = st_drop_geometry(mills[,"timetop1"]),
                                #bins = 4, 
                                na.color = "transparent", 
                                reverse = F)
leaflet() %>% 
  addTiles()%>%
  addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%

  addPolygons(data = mills, 
              opacity = 1, weight = 2, col = ~cb_timetop1(mills$timetop1), 
              fill = TRUE, fillColor = ~cb_timetop1(mills$timetop1),
              label = mills$timetop1) %>% 
  addMarkers(data = m.df_wide_lonlat[rowid,])


