
### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table", 
                   "foreign", "readstata13", "readxl",
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "fixest", "sandwich", "lmtest", "boot", 
                   "ggplot2") 
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




#### PREPARATION
parcel_size <- 3000
catchment_radius <- 3e4
island <- c("Sumatra", "Kalimantan", "Papua")
island <- "Sumatra"
outcome_variable <- "lucpfip_pixelcount_total"
dynamics <- FALSE
commo <- c("ffb", "cpo")#
yoyg <- FALSE
short_run <- "unt level" # from c("unt_level", "dev", "yoyg")
imp <- 1
x_pya <- 3 # 2, 3 or 4
lag_or_not <- "_lag1" # from c("_lag1", "")

pya_ov <- FALSE
cluster_var <- "parcel_id"
fe_set <- "parcel_id + district^year"# needs to include year (because of the structural derivation of the first stage.
# so either "parcel_id + year" or "parcel_id + district^year" 

### BASIC DATA
d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                              parcel_size/1000,"km_",
                              catchment_radius/1000,"CR.rds")))


## Keep observations that: 

# - are in island or group of islands of interest
d <- d[d$island %in% island,]

# - were covered with some of the studied forest in 2000
if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
  d <- d[d$any_pfc2000_total,]
}
if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
   outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
   outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
  d <- d[d$any_fc2000_30th,]
}


### VARIABLES INVOLVED 
# we need a different kind of specifications as elsewhere, because we cannot have LR measures (X ya variables)
# as instrumented variables. 
endo_var <- c(paste0("wa_cpo_price_imp",imp,lag_or_not))

controls <- c(paste0("wa_ffb_price_imp",imp),
              paste0("wa_ffb_price_imp",imp,"_3pya"),
              paste0("wa_cpo_price_imp",imp,"_3pya"),
              "wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp","n_reachable_uml")
# add pya outcome variable 
if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}

# lag controls or not
if(lag_or_not=="_lag1"){controls <- paste0(controls,lag_or_not)}


#### NORMAL REGRESSION ON NORMAL SAMPLE 
### FORMULAE
direct_fml <- as.formula(paste0(outcome_variable,
                                " ~ ", paste0(endo_var, collapse = "+"), 
                                " + ",
                                paste0(controls, collapse = "+"),
                                " | ",
                                fe_set))

### SAMPLE
# - are kept by fixest, i.e. that have no NA on any of the variables used (in either stage)
filter_vec <- base::rowSums(!is.na(d[,c("parcel_id", "district", "island", "year","sample_coverage_lag1", outcome_variable, endo_var, controls)]))
filter_vec <- filter_vec == length(c("parcel_id", "district", "island", "year","sample_coverage_lag1", outcome_variable, endo_var, controls))
d_nona <- d[filter_vec, c("parcel_id", "district", "island", "year", "sample_coverage_lag1", outcome_variable, endo_var, controls)]
anyNA(d_nona)

# - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
d_nona$district_year <- paste0(d_nona$district,"_",d_nona$year)

if("parcel_id + district^year"%in%fe_set){
d_clean <- d_nona[-fixest::obs2remove(fml = as.formula(paste0(outcome_variable, " ~ parcel_id + district_year")),
                                      d_nona, 
                                      family = "poisson"),]
}else{d_clean <- d_nona[-obs2remove(direct_fml, d_nona, family = "poisson"),]}

# sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
# probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
# just remove them if any, so that there is no bug. 
d_clean <- d_clean[d_clean$sample_coverage_lag1!=0,]


# additionnally remove 2015 just to see
# d_clean <- d_clean[d_clean$year<2015,]

### REGRESSIONS
var_weights <- d_clean$sample_coverage_lag1/100 
dir_est_2nd <- fixest::feglm(direct_fml,
                             data = d_clean,
                             family = "quasipoisson", 
                             weights = var_weights)

# this reproduces the results in pdf if endo_var <- c(paste0("wa_",commo[2],"_price_imp",imp,"_4ya",lag_or_not))
fst_stages <- list()
iv_comp <- list()
for(iv in 1:6){
  instrument <- paste0("iv",iv,"_imp",imp,lag_or_not)
  
#### NORMAL ON IV SAMPLE
# now see how the same compares when the sample is trimmed to what is used in IV strategy

### FORMULAE
# it's the same as above.

### SAMPLE
filter_vec_iv <- base::rowSums(!is.na(d[,c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument)]))
filter_vec_iv <- filter_vec_iv == length(c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument))
d_nona_iv <- d[filter_vec_iv, c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument, "sample_coverage_lag1")]
anyNA(d_nona_iv)

# - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove

d_nona_iv$district_year <- paste0(d_nona_iv$district,"_",d_nona_iv$year)

if("parcel_id + district^year"%in%fe_set){
  d_clean_iv <- d_nona_iv[-fixest::obs2remove(fml = as.formula(paste0(outcome_variable, " ~ parcel_id + district_year")),
                                        d_nona_iv, 
                                        family = "poisson"),]
}else{d_clean_iv <- d_nona_iv[-obs2remove(direct_fml, d_nona_iv, family = "poisson"),]}

# sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
# probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
# just remove them if any, so that there is no bug. 
d_clean_iv <- d_clean_iv[d_clean_iv$sample_coverage_lag1!=0,]

#d_clean_iv <- d_clean_iv[d_clean_iv$year<2015,]

# while d_clean has 11615 observations (for sumatra and  
# [1] "parcel_id"                    "district"                     "year"                         "lucpfip_pixelcount_total"    
# [5] "wa_ffb_price_imp1_4ya_lag1"   "wa_cpo_price_imp1_4ya_lag1"   "wa_pct_own_loc_gov_imp_lag1"  "wa_pct_own_nat_priv_imp_lag1"
# [9] "wa_pct_own_for_imp_lag1"      "n_reachable_uml_lag1"         "iv2_imp1_lag1"                "sample_coverage_lag1" ) 
# d_clean_iv has 7444 


### REGRESSION
# exact same regression as above but on the iv-trimmed sample
var_weights_iv <- d_clean_iv$sample_coverage_lag1/100 
dir_est_2nd_iv <- fixest::feglm(direct_fml,
                             data = d_clean_iv,
                             family = "quasipoisson", 
                             weights = var_weights_iv)


#### IV (ON IV SAMPLE)

### FORMULAE
fml_1st <- as.formula(paste0(endo_var,
                             " ~ ", instrument,
                             " + ",
                             paste0(controls, collapse = "+"),
                             " | ", 
                             fe_set))

fml_2nd <- as.formula(paste0(outcome_variable,
                             " ~ ", endo_var,
                             " + est_res_1st +", 
                             paste0(controls, collapse = "+"),
                             " | ",
                             fe_set))

### SAMPLE
filter_vec_iv <- base::rowSums(!is.na(d[,c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument)]))
filter_vec_iv <- filter_vec_iv == length(c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument))
d_nona_iv <- d[filter_vec_iv, c("parcel_id", "district", "year", outcome_variable, endo_var, controls, instrument, "sample_coverage_lag1")]
anyNA(d_nona_iv)
any(filter_vec_iv)
# - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove

d_nona_iv$district_year <- paste0(d_nona_iv$district,"_",d_nona_iv$year)

if("parcel_id + district^year"%in%fe_set){
  d_clean_iv <- d_nona_iv[-fixest::obs2remove(fml = as.formula(paste0(outcome_variable, " ~ parcel_id + district_year")),
                                              d_nona_iv, 
                                              family = "poisson"),]
}else{d_clean_iv <- d_nona_iv[-obs2remove(direct_fml, d_nona_iv, family = "poisson"),]}

# sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
# probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
# just remove them if any, so that there is no bug. 
d_clean_iv <- d_clean_iv[d_clean_iv$sample_coverage_lag1!=0,]

### REGRESSIONS
var_weights_iv <- d_clean_iv$sample_coverage_lag1/100 

# # 1st stage 
est_1st <- fixest::feols(fml_1st,
                         data = d_clean_iv, 
                         weights = var_weights_iv)
# save estiamted residuals
d_clean_iv$est_res_1st <- est_1st$residuals
# 2nd stage
est_2nd <- fixest::feglm(fml_2nd,
                         data = d_clean_iv,
                         family = "quasipoisson", 
                         weights = var_weights_iv)

fst_stages[[iv]] <- etable(est_1st, coefstat = "confint")
iv_comp[[iv]] <- etable(dir_est_2nd, dir_est_2nd_iv, est_2nd, coefstat = "confint")
}
iv_comp[[1]]
iv_comp[[2]]
iv_comp[[3]]
iv_comp[[4]]
iv_comp[[5]]
iv_comp[[6]]

fst_stages[[1]]
fst_stages[[2]]
fst_stages[[3]]
fst_stages[[4]]
fst_stages[[5]]
fst_stages[[6]]

unique(d_clean_iv$year)
unique(d_clean$year)






# explore variables in time
ExPanD(df = d_clean, cs_id = "parcel_id", ts_id = "year")






# compare different parcel data sets
# island outputs of prepare_lucpfip.R
d_30km_IBS_CR_suma <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                  "Sumatra_",
                                  parcel_size/1000,"km_",
                                  catchment_radius/1000,"km_IBS_CR_",
                                  "total.rds")))
d_30km_IBS_CR_kali <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                               "Kalimantan_",
                                               parcel_size/1000,"km_",
                                               catchment_radius/1000,"km_IBS_CR_",
                                               "total.rds")))
d_30km_IBS_CR_papu <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                               "Papua_",
                                               parcel_size/1000,"km_",
                                               catchment_radius/1000,"km_IBS_CR_",
                                               "total.rds")))


d_30km_IBS_CR <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                           parcel_size/1000,"km_",
                                           catchment_radius/1000,"km_IBS_CR.rds")))

# output of wa_at_parcels.R
d_30CR <- readRDS(file.path(paste0("temp_data/processed_parcels/wa_panel_parcels_",
                   parcel_size/1000,"km_",
                   catchment_radius/1000,"CR.rds")))

nrow(d_30CR) == nrow(d_30km_IBS_CR)
nrow(d_30CR) == nrow(d_30km_IBS_CR_suma)+nrow(d_30km_IBS_CR_kali)+nrow(d_30km_IBS_CR_papu)

# they have the same number of obs.: 857664 for all islands, 30km_IBS_CR
# BUT NOT FOR THE SAME SET OF YEARS! 
unique(d_30CR$year) # takes the years of IBS
unique(d_30km_IBS_CR$year) # takes the years of GFC
unique(d_30km_IBS_CR_suma$year) # "" 
# they are both balanced panels, so although this is not the same years, this is the same amount of years so 
# sums up to the same number of obs. 
nrow(d_30CR[d_30CR$year>2000,])==nrow(d_30CR[d_30km_IBS_CR$year<2016,]) # 714720
 # and this is also the number of obs. of the final version
d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                              parcel_size/1000,"km_",
                              catchment_radius/1000,"CR.rds")))
nrow(d)

# output of add_parcel_variables.R
d_add_parc <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_final_",
                               parcel_size/1000,"km_",
                               catchment_radius/1000,"CR.rds")))
nrow(d_add_parc)
unique(d_add_parc$year)
nrow(d_add_parc[d_add_parc$year>2000,]) == nrow(d)


# output of merge_lhs_rhs_parcels.R



# there is no grid cell that never has an ibs mill reachable. 
test<- ddply(d,"parcel_id",summarise,
             n_ = sum(n_reachable_ibs))
test[test$n_==0,]
# however, some years, a grid cell may be out of reach. 
any(d$n_reachable_ibs==0)



#### RSPO #### 
# ### add the rspo data
rspo <- st_read("input_data/RSPO_supply_bases/RSPO-certified_oil_palm_supply_bases_in_Indonesia.shp")
#pnas <- st_read("unused_data/idn_rspo_supply_bases_pnas/idn_rspo_supply_bases_pnas.shp")
#unique(pnas$rspoMember)
#plot(pnas[,"rspoMember"])
names(rspo)
unique(rspo$icletdate)
unique(rspo$mill)


plot(st_geometry(rspo))
rspo[,"icletdate"] %>% plot()

# this is just Kali. 
mm <-  read.csv(file.path("unused_data/data_rspo_leakage/master_mill_data.csv"))
lk <-  read.csv(file.path("unused_data/data_rspo_leakage/long_kali.csv"))


#######################################
##### ADD OUTCOME AT NEIGHBORING GRID CELLS ##### 


outcome_variable = "lucpfip_pixelcount_total"
dynamics = FALSE
commo = c("ffb", "cpo")#
yoyg = FALSE
short_run = "unt level"
imp = 1
distribution = "quasipoisson"
fe = c("parcel_id + district_year")
x_pya = 3
lag_or_not = "_lag1"

pya_ov = FALSE 
weights = FALSE
cluster = "cluster"
spatial_control <- TRUE

island_list <- list("Sumatra", "Kalimantan")#, c("Sumatra", "Kalimantan", "Papua")

elmts <- rep(list(NA), length(island_list))
names(elmts) <- c("Sumatra", "Kalimantan")#, "All"

attributes <- c("est_results",
                "distr_lvl", "distr_cf")

i <- 1

for(i in 1:length(island_list)){
  
  controls = c("wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
               "n_reachable_uml")
  
  island <- island_list[[i]]
  
  elmts[[i]] <- list()
  length(elmts[[i]]) <- length(attributes)
  names(elmts[[i]]) <- attributes
  
  
  ### DATA FOR REGRESSIONS
  if(island == "Sumatra"){
    catchment_radius <- 3e4
  }else{
    catchment_radius <- 5e4}
  
  d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))
  
  ## Keep observations that: 
  
  # - are in island or group of islands of interest
  d <- d[d$island %in% island,]
  
  # - were covered with some of the studied forest in 2000
  if(outcome_variable == "lucpfip_pixelcount_total" | outcome_variable == "lucpfip_ha_total"){
    d <- d[d$any_pfc2000_total,]
  }
  if(outcome_variable == "lucfip_pixelcount_30th" | outcome_variable == "lucfip_ha_30th" |
     outcome_variable == "lucfip_pixelcount_60th" | outcome_variable == "lucfip_ha_60th" |
     outcome_variable == "lucfip_pixelcount_90th" | outcome_variable == "lucfip_ha_90th"){
    d <- d[d$any_fc2000_30th,]
  }
  
  if(spatial_control){
    
    # restrict the data set to years where the paste outcome is available
    d <- d[d$year > 2004,]
    
    # d is a balanced panel so far so we can take any first record of a parcel to 
    # keep the most general cross section
    nrow(d) == length(unique(d$parcel_id))*length(unique(d$year))
    d_cs <- d[!duplicated(d$parcel_id),c("parcel_id", "year", "lat", "lon", "island",outcome_variable)]
    
    # spatial
    d_cs <- st_as_sf(d_cs, coords = c("lon", "lat"), remove = FALSE, crs = indonesian_crs)
    
    # identify neighbors
    d_buf <- st_buffer(d_cs, dist = parcel_size - 10)
    row.names(d_buf) <- d_buf$parcel_id
    sgbp <- st_intersects(d_buf)
    
    # iterate on parcel ~30 minutes
    d[,"ngb_ov_lag3"] <- rep(NA, nrow(d))
    d[,"ngb_ov_lag4"] <- rep(NA, nrow(d))
    for(p_i in 1:length(sgbp)){
      # remove own index
      p_i_neighbors <- sgbp[[p_i]][sgbp[[p_i]] != p_i]
      # get the parcel_id corresponding to the sgbp index
      p_i_neighbors <- as.numeric(attr(sgbp,"region.id")[p_i_neighbors])
      # compute mean of outcome variable lags in the PANEL dataset
      for(y in unique(d$year)){
        d[d$parcel_id == as.numeric(attr(sgbp,"region.id")[p_i]) & d$year == y, "ngb_ov_lag4"]  <- mean(d[d$parcel_id %in% p_i_neighbors & d$year == y, paste0(outcome_variable,"_lag4")],na.rm = TRUE)
      }
    } 
    
  # empty neighbor sets (those that have no neighbors but themselves) and  end up with NaN as mean(integer(0)) 
  # turn them into NA
  d[is.nan(d$ngb_ov_lag3),"ngb_ov_lag3"] <- NA
  d[is.nan(d$ngb_ov_lag4),"ngb_ov_lag4"] <- NA
    
  }
  
  #d[d$parcel_id %in% p_i_neighbors & d$year == y, c("parcel_id", "year", paste0(outcome_variable,"_lag3"))]

  # - are not in an intended RSPO-certified supply base
  d_norspo <- d[d$rspo_cert==FALSE,]
  
  ### VARIABLES INVOLVED 
  
  # variable of interests (called regressors pour l'instant)
  if(dynamics == FALSE){ 
    # Only the overall effect is investigated then, no short vs. long run
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
  
  # if, on the other hand, we want to identify the short run effects, controlling for the long run's. 
  if(dynamics == TRUE){ 
    
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- c(paste0("wa_",commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                        paste0("wa_",commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
      }
      if(yoyg == FALSE){
        if(short_run == "unt level"){
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
        if(short_run == "unt level"){
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
  
  # add pya outcome variable 
  if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}
  
  # lag controls or not
  if(lag_or_not=="_lag1"){controls <- paste0(controls,lag_or_not)}
  
  if(spatial_control){controls <- c(controls, "ngb_ov_lag4")}
  
  # add rspo_cert which is not to be lagged. 
  # controls <- c(controls, "rspo_cert")
  # - are kept by fixest, i.e. that have no NA on any of the variables used (in either stage)
  used_vars <- c("parcel_id", "year", "lat", "lon", "district", "island", 
                 "n_reachable_ibsuml_lag1", "sample_coverage_lag1", 
                 outcome_variable, regressors, controls)
  
  filter_vec <- base::rowSums(!is.na(d_norspo[,used_vars]))
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d_norspo[filter_vec, used_vars]
  if(anyNA(d_nona)){stop()}
  
  # - sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
  # probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
  # just remove them if any, so that there is no bug. 
  
  if(weights == TRUE){
    d_nona <- d_nona[d_nona$sample_coverage_lag1!=0,]
  }
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  d_nona$district_year <- paste0(d_nona$district,"_",d_nona$year)
  
  d_clean <- d_nona[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
                                d_nona, 
                                family = "poisson"),]
  
  # note that obs2remove has to be done at the end, otherwise some units (along fe dimension)
  # can become "only-zero-outcome" units after other data removal.  
  
  ### FORMULA
  # list to be filled with formulae of different fixed-effect models
  
  fe_model <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                paste0(regressors, collapse = "+"),
                                " + ",
                                paste0(controls, collapse = "+"),
                                " | ",
                                fe))
  ### REGRESSIONS
  
  # run regression for each fixed-effect model 
  if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
    if(weights == TRUE){
      var_weights <- d_clean$sample_coverage_lag1/100 
      fe_reg <- fixest::feglm(fe_model,
                              data = d_clean, 
                              family = distribution, 
                              notes = TRUE, 
                              weights = var_weights)
      
    }else{
      fe_reg <- fixest::feglm(fe_model,
                              data = d_clean, 
                              family = distribution, 
                              notes = TRUE)
    }
  }else{ # no weights allowed in negative binomial
    fe_reg <- fixest::fenegbin(fe_model,
                               data = d_clean, 
                               family = distribution, 
                               notes = TRUE)
  }
  controls <- controls[-5]
  collinearity(fe_reg)
  d_clean[,c("parcel_id", "year", "ngb_ov_lag4")]
  # store it 
  elmts[[i]][["est_results"]] <- fe_reg
}



# sgbplt <- sgbpl[1:5]
# x <- sgbplt[3]
# x[[1]][x[[1]] != match(x, sgbplt)]

sgbpl2 <- lapply(sgbpl, function(x){x[x != match(list(x),sgbpl)]})

#plot(parcels_buf[parcels_buf$island == "Sumatra", "geometry"][1:20,])

max(lengths(sgbp))




# remove own index 
# there is no duplicate 

sgbp2 <- lapply(sgbp, FUN = function(x){x[-base::match(x, sgbp)]})

sgbpt <- sgbp[1:5]
lapply(sgbpt, function(x){x[-base::match(list(x),sgbpt)]})
sgbp[]
sgbp[[1]][-1]

x %in% sgbp
pa
parcels[1:500,]
plot(st_geometry(island_sf_prj))
plot(parcels$geometry[c(1:190)], add = TRUE, col = "red")




