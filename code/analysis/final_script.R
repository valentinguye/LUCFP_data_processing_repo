
# PACKAGES & OBJECTS -------------------


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("tibble", "plyr", "dplyr", "data.table", "stringr",
                   "foreign", "readstata13", "readxl",
                   "raster", "sp", "spdep", "sf",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "car",  "fixest", "sandwich", "boot", "multcomp", "urca",# 
                   "ggplot2", "viridis")#,"leaflet", "htmltools"
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


### NEW FOLDERS USED IN THIS SCRIPT 
dir.create("temp_data/reg_results")

### SET NUMBER OF THREADS USED BY {FIXEST} TO 1 (might be neccessary to avoid R session crash) 
getFixest_nthreads()

### INDONESIAN CRS necessary for some spatial operations
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE (data available only for this size currently)
parcel_size <- 3000

# PIXEL AREA
# to rescale the average partial effects and the predictions from pixel counts to hectares
pixel_area_ha <- (27.8*27.6)/(1e4)

### Prepare polygons of three Indonesian islands of interest 
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"
island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

### READ ALL POSSIBLE DATASETS HERE
# They are outputs of merge_lhs_rhs_parcels.R

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


# FUNCTIONS -------------------------------------

## REGRESSION FUNCTION ---------------------------------------------------------------
# # Commented out below are the arguments of the regression making function.
# # They may be useful to run parts of the operations within the function.
catchment = "CR"
PS = 3000
outcome_variable = "lucpfip_pixelcount"
island = "both"
start_year = 2002
end_year = 2014
alt_cr = FALSE
nearest_mill = FALSE
margin = "both"
restr_marg_def = TRUE
reduced_form_iv = FALSE
instru_share = "contemp"
commo = c("cpo")
x_pya = 3
dynamics = FALSE
annual = FALSE
price_variation = FALSE
log_prices = TRUE
yoyg = FALSE
only_sr = FALSE
short_run = "full"
imp = 1
distribution = "quasipoisson"
fe = "reachable + district_year"#
offset = FALSE
lag_or_not = "_lag1"
controls = c("n_reachable_uml", "illegal2", "pct_pfc2000_total")#, "wa_prex_cpo_imp1""wa_pct_own_loc_gov_imp",
remaining_forest = FALSE
control_lncpo = FALSE # should the CPO price treatment of interest be included in the controls (for regressions on FFB price)
interaction_terms = c("illegal2", "pct_pfc2000_total") #  NULL # "illegal2"  #c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1")
interact_regressors = TRUE
interacted = "regressors"
pya_ov = FALSE
illegal = "all"# "ill2" #
weights = FALSE
min_forest_2000 = 0
min_coverage = 0
n_iter_glm = 200
output_full = FALSE

# rm(catchment,outcome_variable,island,alt_cr,commo,x_pya,dynamics,log_prices,yoyg,short_run,imp,distribution,fe,remaining_forest,offset,lag_or_not,controls,interaction_terms ,interacted,pya_ov,illegal, nearest_mill, weights)



make_base_reg <- function(island,
                          start_year = 2002, 
                          end_year = 2014, 
                          outcome_variable = "lucpfip_pixelcount", # LHS. One of "lucfip_pixelcount", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount"p
                          PS = 3000,
                          catchment = "CR",  
                          alt_cr = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                          nearest_mill = FALSE, # whether the ibs variables should be attributed to parcels as from the nearest mill or inverse distance weighted average.  
                          margin = "both", # "both"
                          restr_marg_def = TRUE, # should the restrictive definition of extensive margin be used, if margin is either "intensive" of "extensive"
                          reduced_form_iv = FALSE, # should the regressor be the instrument, i.e. the interaction of wa cpo export shares and the spread. 
                          instru_share = "contemp", # if the above is TRUE, which export share should be used, one of "contemp", "lagged", or "avged". 
                          commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo"), commodities the price signals of which should be included in the RHS
                          x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                          dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                          annual = FALSE,
                          price_variation = FALSE, # should the regressors be price variation over the past years, or average (the default)
                          yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                          only_sr = FALSE,
                          log_prices = TRUE, # Logical, should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
                          short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                          imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                          distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                          fe = "reachable + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                          offset = FALSE, # Logical. Should the log of the remaining forest be added as an offset.  
                          lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                          controls = c("n_reachable_uml"), # , "wa_prex_cpo_imp1"character vectors of names of control variables (don't specify lags in their names)
                          control_lncpo = FALSE, # should the CPO price treatment of interest be included in the controls (for regressions on FFB price)
                          remaining_forest = FALSE, # Logical. If TRUE, the remaining forest is added as a control
                          interaction_terms = NULL, # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1"), # may be one or several of the controls specified above. 
                          interacted = "regressors",
                          interact_regressors = TRUE, # if there are two regressors (e.g. ffb and cpo), should their interaction be included in the model? 
                          pya_ov = FALSE, # logical, whether the lagged (by one year) outcome_variable should be added in controls
                          illegal = "all", # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                          min_forest_2000 = 0, 
                          min_coverage = 0, # fraction, from 0 to 1. Minimum share of reachable IBS over all reachable (UML), for an obs.to be included in sample. 
                          weights = FALSE, # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
                          n_iter_glm = 100,
                          output_full = FALSE # if TRUE, the larger dataset (not only the one prepared foranalysis within the function) is returned as a third element of the list output 
){
  
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
    
    # this is just for convenience, because this variable was not created for the CA stream, but is required to exist below 
    # (although wont be effectively used, see in specification charts). 
    # d$reachable <- rep(1, nrow(d))
  }
  
  if(PS == 1000){
    d <- readRDS("temp_data/panel_parcels_ip_final_1km_30CR.rds")
  }
  if(PS == 5000){
    d_30_5km <- readRDS("temp_data/panel_parcels_ip_final_5km_30CR.rds")
    d_50_5km <- readRDS("temp_data/panel_parcels_ip_final_5km_50CR.rds")
    
    d = 
      rbind(
        d_30_5km %>% filter(island == "Sumatra"),
        d_50_5km %>% filter(island == "Kalimantan")
      )
    rm(d_30_5km, d_50_5km)
  }
  
  ### MAKE THE lonlat 
  d$lonlat <- d$lonlat
  
  # This would also be possible if done here, once we have stacked/merged all datasets to the final one. 
  # But changes nothing and is more intricate... 
  # d <- dplyr::arrange(d, lonlat, year)
  # uni_lonlat <- unique(d$lonlat)
  # d <- mutate(d,
  #             lonlat = match(lonlat, uni_lonlat))
  # rm(uni_lonlat)
  
  
  ### SPECIFICATIONS  
  
  # Build the outcome and sample of unregulated deforestation, i.e. illegal indus or smallholder
  if(grepl("lucpfunrp", outcome_variable)){
    if(illegal %in% c("all", "illegal2")){
    d <- 
      d %>% 
      dplyr::mutate(
      lucpfunrp_pixelcount = if_else(
        # in illegal pixels, count indus + smallholders (all)
        (!is.na(illegal2) & illegal2),
        lucpfap_pixelcount,
        # otherwise (legal or unknown legal status), count only smallholders
        lucpfsmp_pixelcount
        )) %>% 
        # AND REMOVE OBS. WITH NO SMALLHOLDERS NOR ILLEGAL 
        # i.e. keep obs either illegal, or legal but with smallholders
        group_by(lonlat) %>% 
        mutate(any_sm_in_cell = any(lucpfsmp_pixelcount>0)) %>% 
        ungroup() %>% 
        filter((!is.na(illegal2) & illegal2) | (!is.na(illegal2) & !illegal2 & any_sm_in_cell)) %>% 
        as.data.frame()
      # equivalent: 
        # lucpfunrp_pixelcount = if_else(
        #   # in illegal pixels, count indus + smallholders
        #   is.na(illegal2) | !illegal2, 
        #   lucpfsmp_pixelcount,
        #   # otherwise (legal and unknown legal status), count only smallholders
        #   lucpfap_pixelcount)
      }
    # d %>% filter(illegal2) %>% pull(illegal2) %>% summary()
    # d %>% filter(is.na(illegal2) | !illegal2) %>% pull(illegal2) %>% summary()
    # d %>% filter(is.na(illegal2) | !illegal2) %>% pull(lucpfunrp_pixelcount) %>% summary()
    # d %>% filter(is.na(illegal2) | !illegal2) %>% pull(lucpfsmp_pixelcount) %>% summary()
    # 
    # d %>% filter(illegal2) %>% pull(lucpfunrp_pixelcount) %>% summary()
    # d %>% filter(illegal2) %>% pull(lucpfap_pixelcount) %>% summary()
    
    if(illegal == "illegal1"){
      d <- 
        d %>% 
        dplyr::mutate(
          lucpfunrp_pixelcount = if_else(
            # in illegal pixels, count indus + smallholders (all)
            (!is.na(illegal1) & illegal1),
            lucpfap_pixelcount,
            # otherwise (legal or unknown legal status), count only smallholders
            lucpfsmp_pixelcount
          )) %>% 
        # AND REMOVE OBS. WITH NO SMALLHOLDERS NOR ILLEGAL 
        # i.e. keep obs either illegal, or legal but with smallholders
        group_by(lonlat) %>% 
        mutate(any_sm_in_cell = any(lucpfsmp_pixelcount>0)) %>% 
        ungroup() %>% 
        filter((!is.na(illegal1) & illegal1) | (!is.na(illegal1) & !illegal1 & any_sm_in_cell)) %>% 
        as.data.frame()
    }
    if(illegal == "alt"){
      d <- 
        d %>% 
        dplyr::mutate(
          lucpfunrp_pixelcount = if_else(
            # in illegal pixels, count indus + smallholders (all)
            (!is.na(illegal2_2020) & illegal2_2020),
            lucpfap_pixelcount,
            # otherwise (legal or unknown legal status), count only smallholders
            lucpfsmp_pixelcount
          )) %>% 
        # AND REMOVE OBS. WITH NO SMALLHOLDERS NOR ILLEGAL 
        # i.e. keep obs either illegal, or legal but with smallholders
        group_by(lonlat) %>% 
        mutate(any_sm_in_cell = any(lucpfsmp_pixelcount>0)) %>% 
        ungroup() %>% 
        filter((!is.na(illegal2_2020) & illegal2_2020) | (!is.na(illegal2_2020) & !illegal2_2020 & any_sm_in_cell)) %>% 
        as.data.frame()
    }
  }
  
  if(grepl("lucfunrp", outcome_variable)){
      d <- 
        d %>% 
        dplyr::mutate(
          lucfunrp_pixelcount = if_else(
            # in illegal pixels, count indus + smallholders (all)
            (!is.na(illegal2) & illegal2),
            lucfap_pixelcount,
            # otherwise (legal or unknown legal status), count only smallholders
            lucfsmp_pixelcount
          )) %>% 
        # AND REMOVE OBS. WITH NO SMALLHOLDERS NOR ILLEGAL 
        # i.e. keep obs either illegal, or legal but with smallholders
        group_by(lonlat) %>% 
        mutate(any_sm_in_cell = any(lucfsmp_pixelcount>0)) %>% 
        ungroup() %>% 
        filter((!is.na(illegal2) & illegal2) | (!is.na(illegal2) & !illegal2 & any_sm_in_cell)) %>% 
        as.data.frame()
    }
  
  
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
        regressors <- paste0(commo,"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not)
      }else{
        regressors <- paste0(commo,"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)
      }
    }
    # if we don't omit a commodity. 
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg_",x_pya+1,"ya",lag_or_not))
        
      }else{
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)) 
      }
    }
  }
  
  # if, on the other hand, we want to identify the short run effects, controlling for the medium run's. 
  if(dynamics == TRUE){ 
    
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo,"_price_imp",imp,"_yoyg",lag_or_not), # SR measure
                        paste0(commo,"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not)) # LR measure          
      }
      if(yoyg == FALSE){
        if(short_run == "full"){
          regressors <- c(paste0(commo,"_price_imp",imp,lag_or_not), # SR measure
                          paste0(commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
        }
        if(short_run == "dev"){
          regressors <- c(paste0(commo,"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not), # SR measure
                          paste0(commo,"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # LR measure
        }  
      }
    }
    
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0(commo[1],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg_",x_pya,"pya",lag_or_not))    
      }
      if(yoyg == FALSE){
        if(short_run == "full"){
          regressors <- c(paste0(commo[1],"_price_imp",imp,lag_or_not),# FFB SR measure
                          paste0(commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not), # FFB LR measure 
                          paste0(commo[2],"_price_imp",imp,lag_or_not),# CPO SR measure
                          paste0(commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not)) # CPO LR measure 
        }
        if(short_run == "dev"){
          regressors <- c(paste0(commo[1],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                          paste0(commo[1],"_price_imp",imp,"_",x_pya,"pya",lag_or_not),
                          paste0(commo[2],"_price_imp",imp,"_dev_",x_pya,"pya",lag_or_not),
                          paste0(commo[2],"_price_imp",imp,"_",x_pya,"pya",lag_or_not))
        }
      }
    }
  }
  
  # if we want every annual price to be a regressor
  if(annual){
    regressors <- paste0(commo,"_price_imp",imp,"_lag",c(1:(x_pya+1)))
  }
  
  # if we use the IV instead of the price, in a reduced form regression
  if(reduced_form_iv){
    
    regressors <- paste0("wa_iv_",instru_share,"_imp",imp, "_lag", c(1:(x_pya+1)))
    
    if(!annual){
      # c'st un peu du bricolage pour l'instant, a terme on construit cette moyenne dans add_CR_parcel_variables.R
      newvarname <- paste0("wa_iv_",instru_share,"_imp",imp,"_",x_pya+1,"ya",lag_or_not)
      d <- dplyr::mutate(d, !!as.symbol(newvarname) := rowMeans(across(.cols = any_of(regressors)), na.rm = TRUE))
      regressors <- newvarname
    }
  }
  
  if(only_sr){
    if(length(commo) == 1){
      if(yoyg == TRUE){
        regressors <- paste0(commo,"_price_imp",imp,"_yoyg",lag_or_not)
      }else{
        regressors <- paste0(commo,"_price_imp",imp,lag_or_not)
      }
    }
    # if we don't omit a commodity. 
    if(length(commo) == 2){
      if(yoyg == TRUE){
        regressors <- c(paste0(commo[1],"_price_imp",imp,"_yoyg",lag_or_not),
                        paste0(commo[2],"_price_imp",imp,"_yoyg",lag_or_not))
        
      }else{
        regressors <- c(paste0(commo[1],"_price_imp",imp,lag_or_not),
                        paste0(commo[2],"_price_imp",imp,lag_or_not)) 
      }
    }
  }
  
  if(price_variation){
    if(length(commo) == 1){
      regressors <- paste0(commo,"_price_imp",imp,"_",x_pya+1,"yv",lag_or_not)
    }
    if(length(commo) == 2){
      regressors <- c(paste0(commo[[1]],"_price_imp",imp,"_",x_pya+1,"yv",lag_or_not), 
                      paste0(commo[[2]],"_price_imp",imp,"_",x_pya+1,"yv",lag_or_not))
    }
  }
  
  
  if(!nearest_mill & !reduced_form_iv){
    regressors <- paste0("wa_",regressors)
  }
  
  # Logarithms
  if(log_prices & !price_variation){
    for(reg in regressors){
      d[,paste0("ln_",reg)] <- log(d[,reg])
    }
    regressors <- paste0("ln_", regressors)
    rm(reg)
  }
  
  ## CONTROLS
  # Make illegal-legal dummies
  if("illegal2"%in%controls){
    d = d %>% mutate(illegal2 = if_else(illegal2, 1, 0))
  }
  if("illegal1"%in%controls){
    d = d %>% mutate(illegal1 = if_else(illegal1, 1, 0))
  }
  if("ill_or_concession" %in% controls){
    d = 
      d %>% 
      mutate(ill_or_concession = case_when(
        illegal2 & !is.na(illegal2) ~ 1, 
        !illegal2010 & !is.na(illegal2010) ~ 0, # this is parcels in concessions in 2020
        TRUE ~ NA
      ))
    # d$ill_or_concession %>% summary()
  }
  if("ill_or_leg" %in% controls){
    d = 
      d %>% 
      mutate(ill_or_leg = case_when(
        illegal2 & !is.na(illegal2) ~ 1, 
        legal2 & !is.na(legal2) ~ 0, # this is parcels in concessions in 2020
        TRUE ~ NA
      ))
    # d$ill_or_leg %>% summary()
  }


  # add lagged outcome variable 
  #if(pya_ov){controls <- c(controls, paste0(outcome_variable,"_",x_pya,"pya"))}
  # for 1 year lag
  if(any(grepl(pattern = paste0(outcome_variable,"_lag1"), x = controls))){
    d <- dplyr::arrange(d, lonlat, year)
    d <- DataCombine::slide(d,
                            Var = outcome_variable,
                            TimeVar = "year",
                            GroupVar = "lonlat",
                            NewVar = paste0(outcome_variable,"_lag1"), # the name passed in controls should correspond. 
                            slideBy = -1,
                            keepInvalid = TRUE)
    d <- dplyr::arrange(d, lonlat, year)
  }
  # for 4 year lag
  # if(any(grepl(pattern = paste0(outcome_variable,"_lag4"), x = controls))){
  #   d <- dplyr::arrange(d, lonlat, year)
  #   d <- DataCombine::slide(d,
  #                           Var = outcome_variable,
  #                           TimeVar = "year",
  #                           GroupVar = "lonlat",
  #                           NewVar = paste0(outcome_variable,"_lag4"), # the name passed in controls should correspond.
  #                           slideBy = -4,
  #                           keepInvalid = TRUE)
  #   d <- dplyr::arrange(d, lonlat, year)
  # }
  if(any(grepl(pattern = paste0("lucpfap_pixelcount_lag4"), x = controls))){ 
    for(LAG in 1:4){
    d <- 
      DataCombine::slide(d,
                         Var = "lucpfap_pixelcount_lag4",
                         TimeVar = "year",
                         GroupVar = "lonlat",
                         NewVar = paste0("lucpfap_pixelcount_lag4", "_lag", LAG), 
                         slideBy = -LAG,
                         keepInvalid = TRUE)
  }
  d <- 
  d %>% 
  mutate(lucpfap_pixelcount_lag4_4pya = rowMeans(across(.cols = c("lucpfap_pixelcount_lag4_lag1", # t-5
                                                                  "lucpfap_pixelcount_lag4_lag2", # t-6
                                                                  "lucpfap_pixelcount_lag4_lag3", # ...
                                                                  "lucpfap_pixelcount_lag4_lag4")), na.rm = FALSE))
  }
  
  if(any(grepl(pattern = paste0("ngb_ov_lag4"), x = controls))){
    for(LAG in 1:4){
    d <- 
      DataCombine::slide(d,
                         Var = "ngb_ov_lag4",
                         TimeVar = "year",
                         GroupVar = "lonlat",
                         NewVar = paste0("ngb_ov_lag4", "_lag", LAG), 
                         slideBy = -LAG,
                         keepInvalid = TRUE)
  }
  d <- 
  d %>% 
  mutate(
  ngb_ov_lag4_4pya = rowMeans(across(.cols = c("ngb_ov_lag4_lag1",
                                               "ngb_ov_lag4_lag2",
                                               "ngb_ov_lag4_lag3",
                                               "ngb_ov_lag4_lag4")), na.rm = FALSE))

  }
  
  # lag controls that are from IBS 
  select_ibs_controls <- grepl("wa", controls) & !grepl("prex", controls) # !grepl("wa_avg_", controls) & !grepl("wa_lag1_", controls) # | grepl("prex_", controls)
  if(lag_or_not=="_lag1" & length(controls)>0 & any(select_ibs_controls)){
    # find them
    controls[select_ibs_controls] <- sapply(controls[select_ibs_controls], FUN = paste0, lag_or_not)
  }
  
  # add remaining forest or not
  if("remaining_forest" %in% controls){
    if(grepl("lucpf", outcome_variable)){
      controls[controls == "remaining_forest"] <- "remain_pf_pixelcount"
    }
    if(grepl("lucf", outcome_variable)){
      controls[controls == "remaining_forest"] <- "remain_f30th_pixelcount"
    }
    
    offset <- FALSE
  }
  
  # build baseline forest trend if it is featured in the controls. 
  if("baseline_forest_trend" %in% controls){
    if(grepl("lucpf", outcome_variable)){
      d[,"baseline_forest_trend"] <- d$pfc2000_total_pixelcount*(d$year-2000)
    }
    if(grepl("lucf", outcome_variable)){
      d[,"baseline_forest_trend"] <- d$fc2000_30th_pixelcount*(d$year-2000)
    }
  }
  
  if(nearest_mill==TRUE){
    controls <- gsub("wa_", "", controls)
  }
  
  # add price level controls, in price variability regression 
  if(price_variation){
    controls <- c(controls, paste0("wa_", commo,"_price_imp",imp,"_",x_pya+1,"ya",lag_or_not))
  }
  # and log cpo price control, in regression on FFB price, if asked. 
  if(control_lncpo){
    # this is the equivalent to how the ffb price treatment is constructed
    cpo_ctrl <- paste0("wa_cpo_price_imp",imp,"_",x_pya+1,"ya",lag_or_not)
    lncpo_ctrl <- paste0("ln_",cpo_ctrl)
    d[,lncpo_ctrl] <- log(d[,cpo_ctrl])
    controls <- c(lncpo_ctrl, controls) # put it first, for make_APEs
  }
  
  # ### WEIGHTS
  # if(weights){
  #   d$sample_coverage <- d$n_reachable_ibs/d$n_reachable_uml
  # }
  
  # Make baseline forest cover share of cell, if it is featured in the interaction terms or controls
  if("pct_pfc2000_total" %in% c(interaction_terms, controls)){
    d = d %>% mutate(pct_pfc2000_total = 100 * pfc2000_total_pixelcount*pixel_area_ha/900)
  }
  
  # INTERACTING OUTCOME SHARES 
  if("share_indus" %in% controls | "share_illindus" %in% controls){
    d = 
      d %>% 
      mutate(share_indus_wna = 100 * lucpfip_pixelcount / lucpfap_pixelcount) %>% 
      group_by(reachable, district_year) %>% 
      mutate(feavg_share_indus = mean(share_indus_wna, na.rm = TRUE)) %>% 
      ungroup() %>% 
      as.data.frame() %>% 
      mutate(share_indus = if_else(lucpfap_pixelcount > 0, share_indus_wna, feavg_share_indus))  # & !is.na(illegal2)
    
    d = 
      d %>% 
      mutate(share_illindus = if_else(lucpfap_pixelcount > 0 | is.na(illegal2), illegal2*share_indus_wna, illegal2*feavg_share_indus)) # & !is.na(illegal2)
    # This means: if share_indus_wna is defined, then either illegal2 is not missing and share_illindus is defined, or illegal2 is missing and share_illindus will be NA. 
    # If share_indus_wna is not defined and illegal2 is not missing, then attribute the average of indus share interacted by illegal. 
    
    # d$share_indus_wna %>% summary()
    # d$feavg_share_indus %>% summary()
    # d$share_indus %>% summary()
    # d$share_illindus %>% summary()
  }
  # if("share_sm" %in% controls){
  #   d = 
  #     d %>% 
  #     mutate(share_sm = if_else(lucpfap_pixelcount > 0, lucpfsmp_pixelcount / lucpfap_pixelcount, 0))
  # }

  ### SELECT DATA FOR REGRESSION
  
  ## group all the variables necessary in the regression
  # important to do that after outcome_variable, regressors controls etc. have been (re)defined. 
  # (interactions do not need to be in there as they are fully built from the used_vars)
  used_vars <- c(outcome_variable, regressors, controls,
                 "lonlat",  "year", "lat", "lon", 
                 "village", "subdistrict", "district", "province", "island", "reachable", # "illegal2", IT IS IMPORTANT THAT LEGAL STATUS be not in the used_vars, because there are many missings so it changes the sample, while it is not needed as we split samples based on it before selecting only used_vars
                 "village_year", "subdistrict_year", "district_year", "province_year")#,"reachable", "reachable_year"
  #"n_reachable_ibsuml_lag1", "sample_coverage_lag1", #"pfc2000_total_ha", 
  #"remain_f30th_pixelcount","remain_pf_pixelcount"
  if(nearest_mill){
    used_vars <- c(used_vars, "nearest_firm_id")
  }
  
  # if we use an offset (for now only remaining forest possible) then it should be in the necessary variables, but not among the controls. 
  if(offset){
    if(grepl("lucpf", outcome_variable)){
      used_vars <- c(used_vars, "remain_pf_pixelcount")
    }
    if(grepl("lucf", outcome_variable)){
      used_vars <- c(used_vars, "remain_f30th_pixelcount")
    }
  } 
  
  ### KEEP OBSERVATIONS THAT: 
  
  # - are in study period 
  d <- dplyr::filter(d, year >= start_year)
  d <- dplyr::filter(d, year <= end_year)
  
  # - are covered with some of the studied forest 
  # in 2000
  if(PS == 3000){ # variables available only for 3km CR
  if(grepl("lucpf", outcome_variable)){
    d <- dplyr::filter(d, pfc2000_total_pixelcount*pixel_area_ha/900 > min_forest_2000)
    #d <- d[d$pfc2000_total_pixelcount*pixel_area_ha/900 > min_forest_2000,]
  }
  if(grepl("lucf", outcome_variable)){
    d <- dplyr::filter(d, fc2000_30th_pixelcount*pixel_area_ha/900 > min_forest_2000)
    # d <- d[d$fc2000_30th_pixelcount*pixel_area_ha/900 > min_forest_2000,]
  }

  # # in years after grid cells get completely deforested) - this is insignificant as there will almost always remain some pixels of forest
  if(grepl("lucpf", outcome_variable)){
    d <- dplyr::filter(d, remain_pf_pixelcount > 0,)
    #d <- d[d$remain_pf_pixelcount > 0,]
  }
  if(grepl("lucf", outcome_variable)){
    d <- dplyr::filter(d, remain_f30th_pixelcount > 0,)
    # d <- d[d$remain_f30th_pixelcount > 0,]
  }
  }
  # # Do not experience industrial and smallholder deforestation at the same time. 
  # if(grepl("lucpf", outcome_variable)){
  #   d <- d[d$lucp_i_and_sm_bar,]
  # }
  # if(grepl("lucf", outcome_variable)){
  #   d <- d[d$luc_i_and_sm_bar,]
  # }
  
  
  # - are in years when the outcome can actually be observed
  # if(grepl("_slow_",outcome_variable)){
  #   d <- dplyr::filter(d, year < 2011)
  #   d <- d[d$year<2011,]
  # }
  
  # remove year 2015 as it is only available for industrial plantations
  # d <- dplyr::filter(d, year < 2015)
  #  d <- d[d$year<2015,]
  # if(grepl("fsmp_", outcome_variable) | grepl("fsp_", outcome_variable) | grepl("fmp_", outcome_variable)){
  #   d <- d[d$year<2015,]
  # }
  
  # - are in intensive margin expansion areas
  if(restr_marg_def){
    if(margin == "intensive"){
      d <- dplyr::filter(d, extensive_restr == FALSE)
      # d <- d[d$extensive_restr == FALSE,]
    }
    if(margin == "extensive"){
      d <- dplyr::filter(d, extensive_restr == TRUE)
      # d <- d[d$extensive_restr == TRUE,]
    }
  } else {
    if(margin == "intensive"){
      d <- dplyr::filter(d, extensive == FALSE)
      # d <- d[d$extensive == FALSE,]
    }
    if(margin == "extensive"){
      d <- dplyr::filter(d, extensive == TRUE)
      # d <- d[d$extensive == TRUE,]
    }
  }
  
  # - have a minimum reachable mill sample coverage (and a positive amount of reachable mills)
  d <- dplyr::filter(d, n_reachable_ibs > 0 &
                       n_reachable_uml > 0 & 
                       n_reachable_ibs/n_reachable_uml >= min_coverage)
  
  # d <- d[d$n_reachable_ibs > 0 &
  #        d$n_reachable_uml > 0 & 
  #        d$n_reachable_ibs/d$n_reachable_uml >= min_coverage,]
  
  
  # - are not in an intended RSPO-certified supply base
  d <- dplyr::filter(d, rspo_cert == FALSE)
  #d <- d[d$rspo_cert==FALSE,]
  
  # are legal, illegal, or neither
  if(illegal == "no_ill1"){
    d <- d[!is.na(d$illegal1) & d$illegal1 == FALSE, ]
  }
  if(illegal == "no_ill2"){
    d <- d[!is.na(d$illegal2) & d$illegal2 == FALSE, ]
  }
  if(illegal == "ill1"){
    d <- d[!is.na(d$illegal2) & d$illegal1 == TRUE, ]
  }  
  if(illegal == "ill2"){
    d <- d[!is.na(d$illegal2) & d$illegal2 == TRUE, ]
  }
  if(illegal == "alt"){
    d <- d[!is.na(d$illegal2_2020) & d$illegal2_2020 == TRUE, ]
  }
  if(illegal == "in_concession"){
    d <- d[!is.na(d$illegal2010) & d$illegal2010 == FALSE, ]
  }
  if(illegal == "out_concession"){
    d <- d[!is.na(d$illegal2010) & d$illegal2010 == TRUE, ]
  }
  
  # make the sample comparable, by requiring that the main regressor be not missing either 
  if(reduced_form_iv){
    used_vars <- c(used_vars, paste0("wa_cpo_price_imp",imp,"_4ya",lag_or_not))
  }
  
  # - have no NA nor INF on any of the variables used (otherwise they get removed by {fixest})
  usable <- lapply(used_vars, FUN = function(var){is.finite(d[,var]) | is.character(d[,var])})
  # used_vars[!(used_vars%in%names(d))]
  names(usable) <- used_vars            
  usable <- bind_cols(usable)
  filter_vec <- base::rowSums(usable)
  filter_vec <- filter_vec == length(used_vars)
  d_nona <- d[filter_vec, c(used_vars)]
  if(anyNA(d_nona)){stop()}
  rm(filter_vec, usable)
  
  # # - sometimes there are a couple obs. that have 0 ibs reachable despite being in the sample 
  # # probably due to some small distance calculation difference between this variable computation and wa_at_parcels.R script. 
  # # just remove them if any, so that there is no bug. 
  # # hm no it's rather panel obs. that don't have any reachable mill YET. So they are important, but they should be
  # if(weights){
  #   d_nona <- d_nona[d_nona$sample_coverage!=0,]
  #   d_nona[d_nona$sample_coverage>1,"sample_coverage"] <- 1
  # }
  
  # OUTLIERS # 
  # otl <- boxplot.stats(d_nona[,regressors[1]])$out
  # d_nona <- d_nona[!(d_nona[,regressors[1]] %in% otl), ]
  
  # This is deprecated. Data are selected according to feglm criteria during feglm call, 
  # and the actual data used for estimation is extracted post estimation
  
  # - and those with not only zero outcome, i.e. that feglm would remove, see ?fixest::obs2remove
  # d_clean <- d_nona[-obs2remove(fml = as.formula(paste0(outcome_variable, " ~ ", fe)),
  #                               d_nona, 
  #                               family = "poisson"),]
  # note that obs2remove has to be the last filtering operation on data, otherwise some units (along fe dimension)
  # may become "only-zero-outcome" units after other data removal.  
  temp_est <- feglm(fml = as.formula(paste0(outcome_variable, " ~ 1 | ", fe)),
                    data = d_nona,
                    glm.iter = 50,
                    family = "poisson")
  # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  if(length(temp_est$obs_selection)>0){
    d_clean <- d_nona[unlist(temp_est$obs_selection),]
  }  else {
    d_clean <- d_nona
  }
  rm(temp_est)
  # d_clean <- d_nona
  
  
  ## INTERACTIONS
  # make the interaction variables in the data
  interaction_vars <- c()
  
  # here we produce the names of the actual interaction variables
  if(length(interaction_terms)>0){
    # unless variables of interest to be interacted are specified, they are all the presently defined regressors
    if(interacted == "regressors"){interacted <- regressors}
    # this function selects the actual controls (lagged if necessary) into the interaction terms
    make_int_term <- function(x){int_term <- controls[grepl(x,controls)] %>% paste0("X",interacted)}
    interaction_vars <- sapply(interaction_terms, FUN = make_int_term)
  }else{interacted <- NULL}
  
  # add the interaction between the regressors 
  if(length(regressors) == 2 & interact_regressors){
    interaction_vars <- c(interaction_vars, paste0(regressors[1],"X",regressors[2]))
  }  
  
  
  # and here we build the interactions 
  if(length(interaction_terms)>0){
    for(con in interaction_terms){
      actual_ <- controls[grepl(con, controls)]
      for(reg in interacted){
        d_clean[,paste0(actual_,"X",reg)] <- d_clean[,actual_]*d_clean[,reg]
        # d_clean = d_clean %>% mutate(!!as.symbol(paste0(actual_,"X",reg)) := !!as.symbol(actual_) * !!as.symbol(reg))# same as above
      }
    }
  }
  # d_clean[,paste0(actual_,"X",reg)] %>% summary()
  
  # and add the interaction between the regressors 
  if(length(regressors) == 2 & interact_regressors){
    d_clean[,paste0(regressors[1],"X",regressors[2])] <- d_clean[,regressors[1]]*d_clean[,regressors[2]]
    #d_clean[,paste0(regressors[2],"X",regressors[1])] <- d_clean[,regressors[2]]*d_clean[,regressors[1]]
  }
  
  # If this is the pooled regression with interaction by illegal, then use exactly the same samples: 
  if("illegal2Xln_wa_cpo_price_imp1_4ya_lag1" %in% interaction_vars){
    # Pooled data 
    if(outcome_variable == "lucpfip_pixelcount"){
      pdf <- rbind(
        res_data_list_full[["both_i_no_ill2"]][[2]],
        res_data_list_full[["both_i_ill2"]][[2]])
    }
    if(outcome_variable == "lucpfip_rapid_pixelcount"){
      pdf <- rbind(
        res_data_list_dyn[["both_rapid_no_ill2"]][[2]],
        res_data_list_dyn[["both_rapid_ill2"]][[2]])
    }
    d_clean = 
      d_clean %>% 
      inner_join(
      pdf %>% dplyr::select(lonlat, year),
      by = c("lonlat", "year")
    )
    stopifnot(nrow(d_clean) == nrow(pdf))
    rm(pdf)
  }  

  
  ### REGRESSIONS
  
  # Model specification
  fe_model <- as.formula(paste0(outcome_variable,
                                " ~ ",
                                paste0(regressors, collapse = "+"),
                                " + ",
                                paste0(controls, collapse = "+"),
                                " | ",
                                fe))
  
  if(length(interaction_terms)>0 | (interact_regressors & length(regressors) == 2)){
    fe_model <- as.formula(paste0(outcome_variable,
                                  " ~ ",
                                  paste0(regressors, collapse = "+"),
                                  " + ",
                                  paste0(controls, collapse = "+"),
                                  " + ",
                                  paste0(interaction_vars, collapse = "+"),
                                  " | ",
                                  fe))
    # 
    # fe_model <- as.formula(paste0(outcome_variable,
    #                               " ~ ",
    #                               paste0(regressors, collapse = "+"),
    #                               " + ",
    #                               "interact(",regressors,",",interaction_terms,")",#",","FALSE",
    #                               " + ",
    #                               paste0(controls, collapse = "+"),
    #                               " | ",
    #                               fe))
  }
  
  
  # Estimation
  if(offset == TRUE){
    if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
      if(weights == TRUE){
        var_weights <- d_clean$sample_coverage
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution,
                                 offset = offset_fml,
                                 glm.iter = n_iter_glm,
                                 notes = TRUE, 
                                 weights = var_weights)
        
      }else{
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution,
                                 offset = offset_fml,
                                 glm.iter = n_iter_glm,
                                 notes = TRUE)
      }
    }else{ # no weights allowed in negative binomial
      reg_res <- fixest::fenegbin(fe_model,
                                  data = d_clean,
                                  offset = offset_fml,
                                  notes = TRUE)
    }
  }  
  
  if(offset == FALSE){
    if(distribution != "negbin"){ # i.e. if it's poisson or quasipoisson or gaussian
      if(weights == TRUE){
        var_weights <- d_clean$sample_coverage
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution, 
                                 glm.iter = n_iter_glm,
                                 #fixef.iter = 100000,
                                 notes = TRUE, 
                                 weights = var_weights)
        
      }else{
        reg_res <- fixest::feglm(fe_model,
                                 data = d_clean, 
                                 family = distribution, 
                                 glm.iter = n_iter_glm,
                                 #fixef.iter = 100000,
                                 notes = TRUE)
      }
    }else{ # no weights allowed in negative binomial
      reg_res <- fixest::fenegbin(fe_model,
                                  data = d_clean,
                                  #fixef.iter = 100000,
                                  notes = TRUE)
    }
  }
  
  # reg_res <- fixest::feglm(lucpfip_pixelcount ~ ln_wa_cpo_price_imp1_4ya_lag1 + n_reachable_uml + illegal2 + i(illegal2, ln_wa_cpo_price_imp1_4ya_lag1, ref = 0) | reachable + district_year,
  #                          data = d_clean, 
  #                          family = distribution, 
  #                          glm.iter = n_iter_glm,
  #                          #fixef.iter = 100000,
  #                          notes = TRUE)
  
  if(output_full){# we use this only at one point, when comparing estimation sample with sample including missing values
    # so we still want to remove always zero from this latter sample
    temp_est <- feglm(fml = as.formula(paste0(outcome_variable, " ~ 1 | ", fe)),
                      data = d,
                      family = "poisson")
    # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
    if(length(temp_est$obs_selection)>0){
      d <- d[unlist(temp_est$obs_selection),]
    }  else { 
      d <- d
    }
    rm(temp_est)
    
    toreturn <- list(reg_res, d_clean, d)
  }else{
    toreturn <- list(reg_res, d_clean)
  }
  
  rm(d, d_nona, d_clean)
  return(toreturn)
  
}

## APE FUNCTIONS -------------------------------------------------------
# helper function that transforms the list of results into a data frame of average partial effects (APEs) and their standard errors (SEs), 
# for the K first regressors in the models fitted by make_base_reg 
# If there are interactions in the models, the APEs (and SEs) of the interaction effects are computed (may not work if K > 1 then)

# res_data <- res_data_list_byplantation[[3]]
# reg_res <- res_data[[1]]
# d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
# rm(res_data)
k = 1
K=1
cumulative <- TRUE
controls_pe = FALSE
# SE = "cluster"
CLUSTER = "reachable"
stddev = FALSE
rel_price_change = 0.01 # sd/m #
abs_price_change = 1
rounding = 2

make_APEs <- function(#res_data, 
                      reg_elm = 1,
                      K=1, 
                      controls_pe = FALSE, 
                      #SE = "cluster", 
                      CLUSTER = "reachable",#"subdistrict",
                      stddev = FALSE, # if TRUE, the PEs are computed for a one standard deviation (after removing variation in the fixed-effect dimensions)
                      rel_price_change = 0.01, 
                      abs_price_change = 1, 
                      rounding = 2){
  
  # store APEs and their deltaMethod statistics in this list 
  dM_ape_roi_list <- list()
  
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working  
  # res_data <- res_data_list_full[[reg_elm]]
  # reg_res <- res_data[[1]]
  # d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation. 
  # rm(res_data)
  
  ## Redefine changes in regressor to one standard deviation if asked 
  if(stddev){
    # remove fixed effect variations from the regressor
    reg_sd <- fixest::feols(fml = as.formula(paste0(
      names(reg_res$coefficients)[1],
      " ~ 1 | ", 
      paste0(reg_res$fixef_vars, collapse = "+"))),
      data = d_clean)
    
    # and take the remaining standard deviation
    rel_price_change <- sd(reg_sd$residuals)
    abs_price_change <- sd(reg_sd$residuals)
  }
  
  
  # identify the nature of different variables
  coeff_names <- names(coef(reg_res))
  # define actual interaction terms (possibly lagged)
  interaction_effects <- coeff_names[grepl(pattern = "X", names(coef(reg_res)))]
  others <- coeff_names[!(grepl(pattern = "X", coeff_names))]
  interaction_terms <- others[paste0(others,"X",coeff_names[1]) %in% interaction_effects]
  
  
  # data POINTS that deltaMethod will grab 
  
  # averages of the interaction terms -the controls that we interacted with the regressor of interest) 
  int_term_avg <- list()
  if(length(interaction_terms) >0){
    for(i in 1:length(interaction_terms)){
      int_term_avg[[i]] <- mean(d_clean[,interaction_terms[i]])
    }
  }
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values) 
  
  
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  
  ## repeat the following for the K regressors of interest 
  linear_ape_fml_list <- list() # to store some results, see at end of loop
  for(k in 1:K){
    
    # add interaction term and average of the other regressor of interest if there are several (K == 2 only for now)
    if(K == 1){
      interaction_terms_k <- interaction_terms
      int_term_avg_k <- int_term_avg
    }
    
    if(K == 2 & k == 1){
      interaction_terms_k <- c(interaction_terms, coeff_names[2])
      int_term_avg_k <- c(int_term_avg, mean(d_clean[,coeff_names[2]]))
    }
    if(K == 2 & k == 2){
      interaction_terms_k <- c(interaction_terms, coeff_names[1])
      int_term_avg_k <- c(int_term_avg, mean(d_clean[,coeff_names[1]]))
    }
    
    
    linear_ape_fml <- paste0(coeff_names[k])
    
    # if there are interactions, the partial effect formula needs to be appended.
    i_t <- 1
    while(i_t<=length(interaction_terms_k)){
      linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[grepl(coeff_names[k],interaction_effects)][i_t],"*",int_term_avg_k[[i_t]])
      i_t <- i_t +1
    }
    rm(i_t)
    
    # the final formula is different depending on the regressor of interest being in the log scale or not. 
    if(grepl("ln_",coeff_names[k])){
      ape_fml_roi <- paste0("((",1+rel_price_change,")^(",linear_ape_fml,") - 1)*100")#*fv_bar*",pixel_area_ha)
    } else{
      ape_fml_roi <- paste0("(exp(",linear_ape_fml,"*",abs_price_change,") - 1)*100")#*fv_bar*",pixel_area_ha)
    } 
    
    dM_ape_roi <- deltaMethod(object = coef(reg_res),
                              vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE,
                              g. = ape_fml_roi,
                              rhs = 0)
    
    # dM_ape_roi <- matrix(c(1,1,3), 1)
    # colnames(dM_ape_roi) = c("Estimate","2.5 %","97.5 %")
    
    row.names(dM_ape_roi) <- NULL
    dM_ape_roi <- as.matrix(dM_ape_roi)
    dM_ape_roi <- dM_ape_roi[,c("Estimate","2.5 %","97.5 %")]
    # dM_ape_roi_list[[k]]

    dM_ape_roi_list[[k]] <- dM_ape_roi
    

    
    ### COMPUTE PARTIAL EFFECTS OF INTERACTION TERMS 
    dM_ape_roi_list[[K + k]] <- list()
    
    if(length(interaction_terms) >0){
      # average of the regressor of interest
      reg_bar <- mean(d_clean[,coeff_names[k]]) 
      
      for(i in 1:length(interaction_terms)){
        
        # there are several parts 
        # 1. the coeff of the interaction effect of interest
        iei <- paste0(interaction_terms[i], "X", coeff_names[k])
        ape_fml1 <- paste0(iei," + ") 
        
        # 2. the linear expression derivative wrt. the regressor of interest 
        ape_fml2 <- paste0("(",linear_ape_fml,")")
        
        # 3. the derivative of the linear expression wrt. the interaction term; that multiplies the previous linear expression
        ape_fml3 <- paste0("*(",
                           coeff_names[coeff_names==interaction_terms[i]]," + ",
                           iei,"*",reg_bar,")")#
        
        # paste everything together; divide by 100 to approximate the log transformation of regressor of interest if it's in the log scale
        if(grepl("ln_", coeff_names[k])){
          ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")*0.01")#,"*",fv_bar,"*",pixel_area_ha)
          
        } else{
          ape_fml_it <- paste0("(",ape_fml1,ape_fml2, ape_fml3,  ")")#,"*",fv_bar,"*",pixel_area_ha)
        }
        
        dM <- deltaMethod(object = coef(reg_res), 
                          vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE,
                          g. = ape_fml_it, 
                          rhs = 0)
        
        row.names(dM) <- NULL
        dM <- as.matrix(dM)
        dM <- dM[,c("Estimate","2.5 %","97.5 %")]
        
        dM_ape_roi_list[[K + k]][[i]] <- dM
        
      } # closes loop on i (interaction terms)
    } # closes condition on length(interaction_terms) > 0
    
    # save the linear derivative formula wrt. the regressor of interest
    linear_ape_fml_list[[k]] <- linear_ape_fml
  }# closes loop over k. 
  
  
  ### COMPUTE THE PARTIAL EFFECTS OF CONTROLS
  
  if(controls_pe & K == 1){ # we provide controls' partial effects only in the simplest case 
   
    for(coi in others[-K]){ # !grepl("price",others)
      
      # make the linear formula 
      ape_fml_coi <- paste0("(exp(",coi,"*",1,") - 1)*100")
      
      if(grepl("ln_",coi)){
        ape_fml_coi <- paste0("((",1+rel_price_change,")^(",coi,") - 1)*100")#*fv_bar*",pixel_area_ha)
      } else{
        ape_fml_coi <- paste0("(exp(",coi,"*",abs_price_change,") - 1)*100")#*fv_bar*",pixel_area_ha)
      } 
      
      dM_ape_coi <- deltaMethod(object = coef(reg_res), 
                                vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                                g. = ape_fml_coi, 
                                rhs = 0)   
      
      row.names(dM_ape_coi) <- NULL
      dM_ape_coi <- as.matrix(dM_ape_coi)
      dM_ape_coi <- dM_ape_coi[,c("Estimate","2.5 %","97.5 %")]
      
      dM_ape_roi_list[[length(dM_ape_roi_list) + 1]] <- dM_ape_coi    
    }
  }
  
  # FINALLY, COMPUTE THE APE OF THE INTERACTION BETWEEN TWO REGRESSORS OF INTEREST
  if(K == 2){
    # there are several parts 
    # 1. the coeff of the interaction effect of interest
    iei <- paste0(coeff_names[1], "X", coeff_names[2]) # the order is important, as in the regression function, the interaction name follows this order with the two first regressors
    ape_fml1 <- paste0(iei," + ") 
    
    # 2. the linear expression derivative wrt. the regressor of interest 
    # linear_ape_fml will be the last computed one in the loop over k above, 
    # so if K = 2, it is computed as coeff_names[2] being the regressor with respect to which the expression is first derivate 
    ape_fml2 <- paste0("(",linear_ape_fml_list[[1]],")")
    
    # 3. and therefore, the derivative in this part should be wrt. coeff_names[1] that multiplies the previous linear expression
    ape_fml3 <- paste0("(",linear_ape_fml_list[[2]],")")
    
    # paste everything together; divide by 100 to approximate the log transformation of regressor of interest if it's in the log scale
    if(any(grepl("ln_", coeff_names[c(1,2)]))){
      ape_fml_it <- paste0("0.01*(",ape_fml1,ape_fml2, "*", ape_fml3, ")")#,"*",fv_bar,"*",pixel_area_ha)
      
    } else{
      ape_fml_it <- paste0("(",ape_fml1,ape_fml2, "*", ape_fml3, ")")#,"*",fv_bar,"*",pixel_area_ha)
    }
    
    dM <- deltaMethod(object = coef(reg_res), 
                      vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                      g. = ape_fml_it, 
                      rhs = 0)
    
    row.names(dM) <- NULL
    dM <- as.matrix(dM)
    dM <- dM[,c("Estimate","2.5 %","97.5 %")]
    dM_ape_roi_list[[length(dM_ape_roi_list) + 1]] <- dM
  }
  
  
  # make a one column matrix with all computed APEs' estimates, LB and HB values. 
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(dM_ape_roi_list)), 
                data = unlist(dM_ape_roi_list))  
  
  
  ## Shape the matrix of APE stats
  
  # round and coerce to character
  # mat[1,] <- mat[1,] %>% round(digits = 2)
  # mat[c(2,3),] <- mat[c(2,3),] %>% round(digits = 2)
  # mat[4,] <- mat[4,] %>% formatC(digits = 0, format = "f")
  # mat[5,] <- mat[5,] %>% formatC(digits = 0, format = "f")
  # 
  # mat[2,] <- paste0("[",mat[2,],"; ",mat[3,],"]")
  # mat <- mat[-3,] %>% as.matrix()
  
  mat <- round(mat, digits = rounding)
  
  k  <- 1
  while(k < nrow(mat)){
    mat[k+1,] <- paste0("[",mat[k+1,],"; ",mat[k+2,],"]")
    k <- k + 3
  } 
  row.names(mat) <- c(rep(c("Estimate","CI","delete"), nrow(mat)/3))
  mat <- mat[row.names(mat)!="delete",] %>% as.matrix()
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Observations"
  mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of clusters
  mat <- rbind(mat, length(unique(d_clean[,CLUSTER])))
  row.names(mat)[nrow(mat)] <- "Clusters"
  mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  # # add a row with the size of the variation 
  # if(grepl("ln_", coeff_names[k])){
  #   mat <- rbind(mat, rel_price_change)
  # } else{
  #   mat <- rbind(mat, abs_price_change)
  # }  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  
  
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, ape_fml_it, dM_ape_roi_list, linear_ape_fml_list)
  return(mat)
}

# note that this function yields a warning that is not important: 
# In rm(coeff_names, interaction_effects, others, interaction_terms,  :
# objet 'ape_fml_it' introuvable

# This function computes the cumulative APE of annual price signals. The formula is the sum of the individual APEs (not the APE of the sum). 
# but yields very similar results btw. 
make_cumulative_APE <- function(#res_data, 
                                reg_elm = 1,
                                cumulative = TRUE, # should the cumulative APE be returned, or the annual ones.
                                CLUSTER = "reachable",# "subdistrict", 
                                stddev = FALSE,
                                rel_price_change = 0.01, 
                                abs_price_change = 1, 
                                rounding = 2){
  # # get estimation results and data 
  # reg_res <- res_data[[1]]
  # d_clean <- res_data[[2]]
  
  # store APEs and their deltaMethod statistics in this list
  dM_ape_roi_list <- list()
  
  ## Redefine changes in regressor to one standard deviation if asked 
  if(stddev){
    # remove fixed effect variations from the regressor
    reg_sd <- fixest::feols(fml = as.formula(paste0(
      names(reg_res$coefficients)[1],
      " ~ 1 | ", 
      paste0(reg_res$fixef_vars, collapse = "+"))),
      data = d_clean)
    
    # and take the remaining standard deviation
    rel_price_change <- sd(reg_sd$residuals)
    abs_price_change <- sd(reg_sd$residuals)
  }
  
  # identify the nature of different variables
  coeff_names <- names(coef(reg_res))
  
  
  ## FORMULA FOR SUM OF APE OF REGRESSORS OF INTEREST 
  
  # select regressors of interest
  roi <- grep(pattern = "price", x = coeff_names, value = TRUE)
  if(length(roi)==0){
    roi <- grep(pattern = "_iv_", x = coeff_names, value = TRUE)
  }
  
  ape_formulas <- c()
  i <- 1
  for(r in roi){
    # the final formula is different depending on the regressor of interest being in the log scale or not. 
    if(grepl("ln_",coeff_names[1])){
      ape_formulas[i] <- paste0("((",1+rel_price_change,")^(",r,") - 1)*100")#*fv_bar*",pixel_area_ha)
    } else{
      ape_formulas[i] <- paste0("(exp(",r,"*",abs_price_change,") - 1)*100")#*fv_bar*",pixel_area_ha)
    }  
    i <- i+1
  }  
  
  if(cumulative){
    
    final_formula <- paste0(ape_formulas, collapse = " + ")
    
    dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                              vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                              g. = final_formula, 
                              rhs = 0)
    
    row.names(dM_ape_roi) <- NULL
    dM_ape_roi <- as.matrix(dM_ape_roi)
    dM_ape_roi_list[[1]] <- dM_ape_roi[,c("Estimate","2.5 %","97.5 %")]
    
  } else { 
    ape_list <- list()
    for(i in 1:length(ape_formulas)){
      dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                                vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                                g. = ape_formulas[i], 
                                rhs = 0) 
      
      row.names(dM_ape_roi) <- NULL
      dM_ape_roi <- as.matrix(dM_ape_roi)
      dM_ape_roi_list[[i]] <- dM_ape_roi[,c("Estimate","2.5 %","97.5 %")]
    }  
    
  }
  
  
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(dM_ape_roi_list)), 
                data = unlist(dM_ape_roi_list))  
  
  mat <- round(mat, digits = rounding)
  
  k  <- 1
  while(k < nrow(mat)){
    mat[k+1,] <- paste0("[",mat[k+1,],"; ",mat[k+2,],"]")
    k <- k + 3
  } 
  row.names(mat) <- c(rep(c("Estimate","CI","delete"), nrow(mat)/3))
  mat <- mat[row.names(mat)!="delete",] %>% as.matrix()
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Observations"
  mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of clusters
  mat <- rbind(mat, length(unique(d_clean[,CLUSTER])))
  row.names(mat)[nrow(mat)] <- "Clusters"
  mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  rm(coeff_names, reg_res, d_clean, 
     ape_formulas, roi, final_formula, dM_ape_roi)
  return(mat)
}


# this is a simpler version of make_APEs function above, that does not handle interaction terms, 
# and has purposedly res_data as an argument (not like the other APE functions)
make_APEs_1regr <- function(res_data, 
                            # reg_elm = 1,
                            #SE = "cluster", 
                            stddev = TRUE,
                            CLUSTER = "reachable",# "subdistrict", 
                            rel_price_change = 0.01, 
                            abs_price_change = 1){
  # Get estimation results and data 
  # For some reasons, how this function is used (within another function) requires to create objects here, unlike for the main APE function. 
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]]
  
  ## Redefine changes in regressor to one standard deviation if asked 
  if(stddev){
    # remove fixed effect variations from the regressor
    reg_sd <- fixest::feols(fml = as.formula(paste0(
      names(reg_res$coefficients)[1],
      " ~ 1 | ", 
      paste0(reg_res$fixef_vars, collapse = "+"))),
      data = d_clean)
    
    # and take the remaining standard deviation
    rel_price_change <- sd(reg_sd$residuals)
    abs_price_change <- sd(reg_sd$residuals)
  }
  
  # identify the nature of different variables
  coeff_names <- names(coef(reg_res))
  # define actual interaction terms (possibly lagged)
  interaction_effects <- coeff_names[grepl(pattern = "X", names(coef(reg_res)))]
  others <- coeff_names[!(grepl(pattern = "X", coeff_names))]
  interaction_terms <- others[paste0(others,"X",coeff_names[1]) %in% interaction_effects]
  
  
  # data POINTS that deltaMethod will grab 
  
  # averages of the interaction terms -the controls that we interacted with the regressor of interest) 
  int_term_avg <- list()
  if(length(interaction_terms) >0){
    for(i in 1:length(interaction_terms)){
      int_term_avg[[i]] <- mean(d_clean[,interaction_terms[i]])
    }
  }
  
  # average fitted values
  fv_bar <- mean(reg_res$fitted.values) 
  
  ## FORMULA FOR APE OF REGRESSOR OF INTEREST 
  
  linear_ape_fml <- paste0(coeff_names[1])
  i_t <- 1
  while(i_t<=length(interaction_terms)){
    linear_ape_fml <- paste0(linear_ape_fml," + ", interaction_effects[grepl(coeff_names[1],interaction_effects)][i_t],"*",int_term_avg[[i_t]])
    i_t <- i_t +1
  }
  
  # the final formula is different depending on the regressor of interest being in the log scale or not. 
  if(grepl("ln_",coeff_names[1])){
    ape_fml_roi <- paste0("((",1+rel_price_change,")^(",linear_ape_fml,") - 1)*100")#*fv_bar*",pixel_area_ha)
  } else{
    ape_fml_roi <- paste0("(exp(",linear_ape_fml,"*",abs_price_change,") - 1)*100")#*fv_bar*",pixel_area_ha)
  }   
  
  dM_ape_roi <- deltaMethod(object = coef(reg_res), 
                            vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE, 
                            g. = ape_fml_roi, 
                            rhs = 0)
  
  
  row.names(dM_ape_roi) <- NULL
  dM_ape_roi <- as.matrix(dM_ape_roi)
  dM_ape_roi <- dM_ape_roi[,c("Estimate","SE","Pr(>|z|)")]
  
  dM_ape_list <- list()
  dM_ape_list[[1]] <- dM_ape_roi
  
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(dM_ape_list)), 
                data = unlist(dM_ape_list))  
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  
  #row.names(mat) <- rep(c("Estimate","SE","p-value"),1+length(interaction_terms))
  
  rm(coeff_names, interaction_effects, others, interaction_terms, reg_res, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, dM_ape_list)
  return(mat)
}

## DESCRIPTIVE FUNCTIONS -------------------------------------------- 
make_desstats <- function(sample_1, sample_2){
  
  # get estimation results and data 
  # this is the exact set of observations used for estimation
  # #reg_res <- res_data[[1]]
  # parcels <- sample_1
  # d <- res_data[[3]]
  # sample_1 <- res_data_list_des[["both_i"]][[2]]
  # sample_2 <- res_data_list_des[["both_sm"]][[2]]
  
  
  # add variables that were not used in the regressions
  d_all <- rbind(d_30_suma, d_50_kali)
  # determine dependent variable
  dep_var_1 <- names(sample_1)[grepl("pixelcount",names(sample_1))]
  dep_var_2 <- names(sample_2)[grepl("pixelcount",names(sample_2))]
  
  # it is important that both have the same name
  # dep_var_1 is always length 1, but dep_var_2 is not if sample_2 is d
  if(length(dep_var_2)>1){
    dep_var_2 <- dep_var_1
  }
  
  
  
  
  
  # group public ownership
  #d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  
  # process a bit sample_2, in case it's the with-NA sample. - this is done in make_base_reg directly now
  # # identify grid cells with null lucfp every year
  # temp_est <- feglm(fml = as.formula(paste0(dep_var_2, " ~ 1 | reachable + district_year")),
  #                   data = sample_2[!is.na(sample_2[,dep_var_2]),], 
  #                   family = "poisson")
  # # it's possible that the removal of always zero dep.var in some FE dimensions above is equal to with the FE currently implemented
  # if(length(temp_est$obs_selection)>0){
  #   sample_2 <- sample_2[unlist(temp_est$obs_selection),]
  # }  else { 
  #   sample_2 <- sample_2
  # }
  # rm(temp_est)
  
  # deprecated code with obs2remove
  # sample_2 <- sample_2[-obs2remove(fml = as.formula(paste0(dep_var_2, " ~ lonlat + district_year")),
  #                                  sample_2[!is.na(sample_2[,dep_var_2]),], 
  #                                  family = "poisson"),]
  
  variables_1 <- c(dep_var_1,
                   "wa_cpo_price_imp1_4ya_lag1",
                   # "wa_pct_own_gov_imp_lag1", "wa_pct_own_nat_priv_imp_lag1", "wa_pct_own_for_imp_lag1",
                   "n_reachable_uml") 
  variables_2 <- c(dep_var_2,
                   "wa_cpo_price_imp1_4ya_lag1",
                   # "wa_pct_own_gov_imp_lag1", "wa_pct_own_nat_priv_imp_lag1", "wa_pct_own_for_imp_lag1",
                   "n_reachable_uml") 
  
  sample_1 <- left_join(sample_1[,c("lonlat","year")], d_all[,c("lonlat","year",variables_1)], by = c("lonlat","year"))#
  sample_2 <- left_join(sample_2[,c("lonlat","year")], d_all[,c("lonlat","year",variables_2)], by = c("lonlat","year"))
  
  # pixelcount to hectares
  sample_1[,dep_var_1] <- sample_1[,dep_var_1]*pixel_area_ha
  sample_2[,dep_var_2] <- sample_2[,dep_var_2]*pixel_area_ha
  
  
  statistics <- c("mean", "std.dev", "median", "min", "max")
  
  ## Des. stats. for sample_1
  des_sample_1 <- matrix(NA, nrow = length(variables_1), ncol = length(statistics))
  row.names(des_sample_1) <- variables_1
  colnames(des_sample_1) <- c(statistics)
  
  for(var in variables_1){
    des_sample_1[var,statistics] <- summarise(sample_1,
                                              mean = mean(get(var), na.rm=TRUE),
                                              std.dev = sd(get(var), na.rm= TRUE),
                                              median = median(get(var), na.rm= TRUE), 
                                              min = min(get(var), na.rm= TRUE),
                                              max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    des_sample_1[var, "median"] <- paste0(des_sample_1[var,"median"]," [",des_sample_1[var, "min"],"; ",des_sample_1[var,"max"],"]")
  }
  
  des_sample_1 <- des_sample_1[,c("mean", "std.dev", "median")]
  
  length(unique(sample_1$lonlat)) %>% paste0(" number of grid cells in sample") %>%  print()
  nrow(sample_1) %>% paste0(" number of observations in sample") %>%  print()
  
  
  ## Des. stats for sample_2
  des_sample_2 <- matrix(NA, nrow = length(variables_2), ncol = length(statistics))
  row.names(des_sample_2) <- variables_2
  colnames(des_sample_2) <- c(statistics)
  
  for(var in variables_2){
    des_sample_2[var,statistics] <- summarise(sample_2,
                                              mean = mean(get(var), na.rm=TRUE),
                                              std.dev = sd(get(var), na.rm= TRUE),
                                              median = median(get(var), na.rm= TRUE), 
                                              min = min(get(var), na.rm= TRUE),
                                              max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    des_sample_2[var, "median"] <- paste0(des_sample_2[var,"median"]," [",des_sample_2[var, "min"],"; ",des_sample_2[var,"max"],"]")
  }
  
  des_sample_2 <- des_sample_2[,c("mean", "std.dev", "median")]
  
  length(unique(sample_2$lonlat)) %>% paste0(" total number of grid cells") %>%  print()
  nrow(sample_2) %>% paste0(" total number of observations") %>%  print()
  
  
  ### add the t-test column
  t_test <- matrix(NA, nrow = length(variables_1), ncol = 1)
  row.names(t_test) <- variables_1
  colnames(t_test) <- "t test"
  
  for(var in variables_1[-1]){ # do not include the first one as it would not be the right one for sample_2
    test <- t.test(x = sample_1[,var],
                   y = sample_2[,var], 
                   alternative = "two.sided", 
                   mu = 0, 
                   paired = F, 
                   var.equal = FALSE)
    
    t_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  # t test on dependent variables
  test <- t.test(x = sample_1[,dep_var_1],
                 y = sample_2[,dep_var_2], 
                 alternative = "two.sided", 
                 mu = 0, 
                 paired = F, 
                 var.equal = FALSE)
  
  t_test[dep_var_1,] <- test$p.value %>% formatC(digits = 3, format = "f")
  # interpretation: if the 95% CI includes 0, then the difference in means between two samples 
  # is not statistically different from 0. Hence, the two samples are "similar" in means 
  # with respect to the variable tested. 
  # (In other words, we cannot reject the null hypothesis that the difference is null
  # -i.e. the two samples are alike - with 95% confidence) 
  
  ### add the Smirnov test
  ks_test <- matrix(NA, nrow = length(variables_1), ncol = 1)
  row.names(ks_test) <- variables_1
  colnames(ks_test) <- "KS test"
  
  for(var in variables_1[-1]){
    test <- ks.test(x = sample_1[,var],
                    y = sample_2[,var], 
                    alternative = "two.sided", 
                    exact = FALSE)
    
    ks_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  test <- ks.test(x = sample_1[,dep_var_1],
                  y = sample_2[,dep_var_2], 
                  alternative = "two.sided", 
                  exact = FALSE)
  
  ks_test[dep_var_1,] <- test$p.value %>% formatC(digits = 3, format = "f")
  # intepretation: if p-value < 0.05 we cannot reject with 95% confidence that 
  # two distributions are equal, implying that they are different. 
  
  
  ### Bind all 
  des_table <- cbind(des_sample_1, des_sample_2, t_test, ks_test)
  # row names
  row.names(des_table) <- c("Deforestation (ha)",
                            "Price signal ($/tCPO)", 
                            #"Public ownership (%)", 
                            #"Domestic private ownership (%)", 
                            #"Foreign ownership (%)", 
                            #"Competition", 
                            "# reachable mills")
  
  return(des_table)
}



make_desstats_simple <- function(sample_1){
  
  # determine dependent variable
  dep_var_1 <- names(sample_1)[grepl("pixelcount",names(sample_1))]
  
  # Other variables to describe
  other_vars <- c("wa_cpo_price_imp1_4ya_lag1",
                  "wa_ffb_price_imp1_4ya_lag1", 
                   # do NOT include ownership actually, because not a main control variable, and not very relevant at pixel level. 
                   # "wa_pct_own_gov_imp_lag1", "wa_pct_own_nat_priv_imp_lag1", "wa_pct_own_for_imp_lag1",
                  "n_reachable_uml",
                  "pfc2000_total_pixelcount", 
                  "fc2000_30th_pixelcount", 
                  "lucfip_pixelcount",
                  "lucfsmp_pixelcount"
                   ) 
  
  # add variables that were not used in the regressions
  d_all <- rbind(d_30_suma, d_50_kali) %>% dplyr::select(lonlat, year, all_of(other_vars))
  # group public ownership
  #d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  # d_all$wa_pct_own_gov_imp_lag1 <- d_all$wa_pct_own_loc_gov_imp_lag1 + d_all$wa_pct_own_cent_gov_imp_lag1
  
  sample_1 <- left_join(sample_1[,c("lonlat","year", dep_var_1)], 
                        d_all, 
                        by = c("lonlat","year"))#
  
  variables_1 <- c(dep_var_1, other_vars)
  
  # pixelcount to hectares
  sample_1[,dep_var_1] <- sample_1[,dep_var_1]*pixel_area_ha
  sample_1[,"pfc2000_total_pixelcount"] <- sample_1[,"pfc2000_total_pixelcount"]*pixel_area_ha
  sample_1[,"fc2000_30th_pixelcount"]   <- sample_1[,"fc2000_30th_pixelcount"]*pixel_area_ha
  sample_1[,"lucfip_pixelcount"]   <- sample_1[,"lucfip_pixelcount"]*pixel_area_ha
  sample_1[,"lucfsmp_pixelcount"]   <- sample_1[,"lucfsmp_pixelcount"]*pixel_area_ha
  
  statistics <- c("mean", "std.dev", "median", "min", "max")
  
  ## Des. stats. for sample_1
  des_sample_1 <- matrix(NA, nrow = length(variables_1), ncol = length(statistics))
  row.names(des_sample_1) <- variables_1
  colnames(des_sample_1) <- c(statistics)
  
  for(var in variables_1){
    des_sample_1[var,statistics] <- summarise(sample_1,
                                              mean = mean(get(var), na.rm=TRUE),
                                              std.dev = sd(get(var), na.rm= TRUE),
                                              median = median(get(var), na.rm= TRUE), 
                                              min = min(get(var), na.rm= TRUE),
                                              max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix() 
    
    if(grepl("price", var)){
      des_sample_1[var,statistics] =
        des_sample_1[var,statistics] %>% 
        as.numeric() %>% 
        round(digits = 0) %>% 
        formatC(drop0trailing = TRUE, 
                format = "fg", flag = "-", zero.print = TRUE)
    } else{
      des_sample_1[var,statistics] =
        des_sample_1[var,statistics] %>% 
        as.numeric() %>% 
        round(digits = 1) %>% 
        formatC(drop0trailing = TRUE, 
                format = "fg", flag = "-", zero.print = TRUE)
    }
    
    # group meand and SD in one single string, for displying issues
    des_sample_1[var, "mean"] <- paste0(des_sample_1[var,"mean"]," (",des_sample_1[var, "std.dev"],")")
    
    # group median min and max in one single string, for displying issues
    des_sample_1[var, "median"] <- paste0(des_sample_1[var,"median"]," [",des_sample_1[var, "min"],"; ",des_sample_1[var,"max"],"]")
  }
  
  des_sample_1 <- des_sample_1[,c("mean", "median")]
  
  length(unique(sample_1$lonlat)) %>% paste0(" number of grid cells in sample") %>%  print()
  nrow(sample_1) %>% paste0(" number of observations in sample") %>%  print()
  
  return(des_sample_1)
}

# We want to produce, for a set of ibs variables, the mean, median, std.dev. min and max statistics, 
# for the geo_sample that we are going to use, and the is_mill sample. 

make_des_table_ibs <- function(ibs_isl){
  
  variables <- c("min_year", 
                 "ffb_price_imp1", "in_ton_ffb_imp1",
                 "cpo_price_imp1", "out_ton_cpo_imp1",
                 "pko_price_imp1", "out_ton_pko_imp1",
                 "prex_cpo_imp1", 
                 "pct_own_cent_gov_imp", "pct_own_loc_gov_imp", "pct_own_nat_priv_imp", "pct_own_for_imp")
  
  statistics <- c("mean", "std.dev", "median", "min", "max")
  
  
  ## Matrix for analysis_sample
  rhs_des_sample <- matrix(NA, nrow = length(variables), ncol = length(statistics))
  row.names(rhs_des_sample) <- variables
  colnames(rhs_des_sample) <- c(statistics)
  
  for(var in variables){
    rhs_des_sample[var,statistics] <- summarise(dplyr::filter(ibs_isl, analysis_sample == TRUE),
                                                mean = mean(get(var), na.rm=T),
                                                std.dev = sd(get(var), na.rm= TRUE),
                                                median = median(get(var), na.rm= TRUE), 
                                                min = min(get(var), na.rm= TRUE),
                                                max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    rhs_des_sample[var, "median"] <- paste0(rhs_des_sample[var,"median"]," [",rhs_des_sample[var, "min"],"; ",rhs_des_sample[var,"max"],"]")
  }
  
  rhs_des_sample <- rhs_des_sample[,c("mean", "std.dev", "median")]
  
  # number of establishments in sub-group
  # N_sample_row <- c(length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"])), 
  #                           rep("", length(statistics))) 
  # rhs_des_sample <- rbind(N_sample_row, rhs_des_sample)
  
  ## Matrix for all ibs mills
  rhs_des_pop <- matrix(NA, nrow = length(variables), ncol = length(statistics))
  row.names(rhs_des_pop) <- variables
  colnames(rhs_des_pop) <- c(statistics)
  
  for(var in variables){
    rhs_des_pop[var,statistics] <- summarise(dplyr::filter(ibs_isl, is_mill == TRUE),
                                             mean = mean(get(var), na.rm=T),
                                             std.dev = sd(get(var), na.rm= TRUE),
                                             median = median(get(var), na.rm= TRUE), 
                                             min = min(get(var), na.rm= TRUE),
                                             max = max(get(var), na.rm= TRUE)) %>% 
      as.matrix()  %>% 
      round(digits = 2) %>% 
      formatC(drop0trailing = TRUE, 
              format = "fg", flag = "-", zero.print = TRUE)
    
    # group median min and max in one single string, for displying issues
    rhs_des_pop[var, "median"] <- paste0(rhs_des_pop[var,"median"]," [",rhs_des_pop[var, "min"],"; ",rhs_des_pop[var,"max"],"]")
  }
  
  rhs_des_pop <- rhs_des_pop[,c("mean", "std.dev", "median")]
  # # number of establishments in sub-group
  # N_pop_row <- c(length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"])), 
  #                   rep("", length(statistics))) 
  # rhs_des_pop <- rbind(N_pop_row, rhs_des_pop)
  
  # bind two groups together
  rhs_des <- cbind(rhs_des_sample, rhs_des_pop)    
  
  # add the t-test column
  t_test <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(t_test) <- variables
  colnames(t_test) <- "t test"
  
  for(var in variables){
    test <- t.test(x = ibs_isl[ibs_isl$analysis_sample==TRUE,var],
                   y = ibs_isl[ibs_isl$is_mill==TRUE,var], 
                   alternative = "two.sided", 
                   mu = 0, 
                   paired = F, 
                   var.equal = FALSE)
    
    t_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  # interpretation: if the 95% CI includes 0, then the difference in means between two samples 
  # is not statistically different from 0. Hence, the two samples are "similar" in means 
  # with respect to the variable tested. 
  # (In other words, we cannot reject the null hypothesis that the difference is null
  # -i.e. the two samples are alike - with 95% confidence) 
  
  # add the Smirnov test
  ks_test <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(ks_test) <- variables
  colnames(ks_test) <- "KS test"
  
  for(var in variables){
    test <- ks.test(x = ibs_isl[ibs_isl$analysis_sample==TRUE,var],
                    y = ibs_isl[ibs_isl$is_mill==TRUE,var], 
                    alternative = "two.sided", 
                    exact = FALSE)
    
    ks_test[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  
  # intepretation: if p-value < 0.05 we cannot reject with 95% confidence that 
  # two distributions are equal, implying that they are different. 
  
  rhs_des <- cbind(rhs_des, t_test, ks_test)
  
  
  
  # row names
  row.names(rhs_des) <- c("First year in IBS", "FFB MUV (USD/ton)", "FFB input (ton)", 
                          "CPO MUV (USD/ton", "CPO output (ton)", 
                          "PKO MUV (USD/ton)", "PKO output (ton)", 
                          "CPO export share (%)", 
                          "Central government ownership (%)", 
                          "Local government ownership (%)", 
                          "National private ownership (%)", 
                          "Foreign ownership (%)")
  
  return(rhs_des)
}


# MAIN REGRESSIONS ------------------------------------------------------------------------

## Table 3. CPO PRICE ELASTICITY ACROSS PLANTATION TYPES -----------------
### REGRESSIONS ### 
# infrastructure to store results
res_data_list_full <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

# size_list <- list("i","sm", "a")
size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")
# ill_status <- c("in_concession", paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      if(ILL == "no_ill2"){# this is necessary to handle some convergence issue
        higher_iter_glm <- 2000 
      } else {
        higher_iter_glm <- 200
      } 
      res_data_list_full[[elm]] <- make_base_reg(island = ISL,
                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                 illegal = ILL,
                                                 n_iter_glm = higher_iter_glm,
                                                 offset = FALSE)
      names(res_data_list_full)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    higher_iter_glm <- 200
    res_data_list_full[[elm]] <- make_base_reg(island = ISL,
                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                               illegal = ILL,
                                               n_iter_glm = higher_iter_glm,
                                               offset = FALSE)
    names(res_data_list_full)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_full[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs() # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
# ape_mat <- lapply(res_data_list_full, FUN = make_APEs) # this was the way to do it until it stopped working, on R 4.5 

ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticities of deforestation across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, # indent = FALSE,
            bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              underline = FALSE, # useless - line under "Deforestation for:" will need to be removed manually in Overleaf 
              width = "10em",
              latex_valign = "m") 


rm(ape_mat)


# ### PARTIAL EFFECTS WITH ALL CONTROLS 
# rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
# ape_mat <- lapply(res_data_list_full, FUN = make_APEs, controls_pe = TRUE) # and for the same reason, this cannot be wrapped in other functions (an environment problem)
# ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
# row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
# ape_mat
# colnames(ape_mat) <- NULL
# 
# options(knitr.table.format = "latex")
# kable(ape_mat, booktabs = T, align = "c",
#       caption = "Price elasticity  and partial effects of control variables on deforestation across Indonesian oil palm sectors") %>% #of 1 percentage change in medium-run price signal
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c(" " = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1),
#                    bold = F,
#                    align = "c") %>%
#   add_header_above(c(" " = 1,
#                      "Industrial plantations" = 3,
#                      "Smallholder plantations" = 3, 
#                      "All" = 3),
#                    align = "c",
#                    strikeout = F) %>%
#   pack_rows("Price elasticity", 1, 2, 
#             italic = TRUE, bold = TRUE)  %>%
#   pack_rows("Partial effects of:", 3, 8,
#             italic = FALSE, bold = TRUE) %>%
#   pack_rows("Domestic private \n mill ownership", 3, 4,
#             italic = TRUE, bold = TRUE)  %>%
#   pack_rows("Foreign mill \n ownership", 5, 6,
#             italic = TRUE, bold = TRUE)  %>%
#   pack_rows("# reachable mills", 7, 8,
#             italic = TRUE, bold = TRUE)  %>%
#   # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "1em", hline_before = TRUE) %>% 
#   column_spec(column = 1,
#               width = "7em",
#               latex_valign = "m") %>% 
#   column_spec(column = c(2:(ncol(ape_mat))),
#               width = "7em",
#               latex_valign = "m") %>% 
#   column_spec(column = ncol(ape_mat)+1,
#               bold = TRUE)

## Table 4. INDUSTRIAL LUC DYNAMICS --------------------------------------------------------------------------------------

## Regressions
# infrastructure to store results
res_data_list_dyn <- list()
elm <- 1

dyn_list <- list("rapid", "slow")

# rapid and slow in comparable times: i.e. up to 2010
for(DYN in dyn_list){
  for(ILL in ill_status){
    res_data_list_dyn[[elm]] <- make_base_reg(island = "both", 
                                              outcome_variable = paste0("lucpfip_",DYN,"_pixelcount"),
                                              illegal = ILL,
                                              end_year = 2010,
                                              offset = FALSE)
    names(res_data_list_dyn)[elm] <- paste0("both_",DYN,"_",ILL)
    elm <- elm + 1
  }
}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_dyn)){
  # get estimation results and data; this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_dyn[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticity of immediate and transitional deforestation in Indonesian industrial plantations") %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "Immediate conversion" = 3,
                     "Transitional conversion" = 3),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 6),
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE,
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "10em",
              latex_valign = "m") #%>%
  # column_spec(column = c(2:(ncol(ape_mat))),
  #             width = "4em",
  #             latex_valign = "c")


rm(res_data_list_dyn)


## Table 5. FFB SMALLHOLDERS -----------------------------------------------------------------------------------------
# infrastructure to store results
res_data_list_ffbsm <- list()
elm <- 1
SIZE <- "sm"
ILL <- "all"

## FFB price, NOT controlling for CPO prices
res_data_list_ffbsm[[elm]] <- make_base_reg(island = "both",
                                            outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                            illegal = ILL,
                                            commo = c("ffb"), 
                                            offset = FALSE)
names(res_data_list_ffbsm)[elm] <- paste0("both_",SIZE,"_",ILL,"_ffb")
elm <- elm + 1


## FFB price, controlling for CPO prices
res_data_list_ffbsm[[elm]] <- make_base_reg(island = "both",
                                            outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                            illegal = ILL,
                                            commo = c("ffb"), # important that ffb is indeed in first position of the vector. 
                                            control_lncpo = TRUE,
                                            offset = FALSE)
names(res_data_list_ffbsm)[elm] <- paste0("both_",SIZE,"_",ILL,"_ffbcpo")
elm <- elm + 1


## FFB price, controlling for CPO prices AND LUCPFIP
res_data_list_ffbsm[[elm]] <- make_base_reg(island = "both",
                                            outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                            illegal = ILL,
                                            commo = c("ffb"), # important that ffb is indeed in first position of the vector. 
                                            control_lncpo = TRUE,
                                            controls = c("n_reachable_uml", "lucpfip_rapid_pixelcount"),
                                            offset = FALSE)
names(res_data_list_ffbsm)[elm] <- paste0("both_",SIZE,"_",ILL,"_ffbcpo")
elm <- elm + 1


## FFB price, controlling for CPO prices, and interacting FFB AND CPO with LUCPFIP  
res_data_list_ffbsm[[elm]] <- make_base_reg(island = "both",
                                            outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                            illegal = ILL,
                                            commo = c("ffb"), # important that ffb is indeed in first position of the vector. 
                                            control_lncpo = TRUE,
                                            controls = c("n_reachable_uml", "lucpfip_rapid_pixelcount"),
                                            interaction_terms = c("lucpfip_rapid_pixelcount"),
                                            offset = FALSE)
names(res_data_list_ffbsm)[elm] <- paste0("both_",SIZE,"_",ILL,"_ffbcpo_interactreachable")
elm <- elm + 1

## PARTIAL EFFECTS
# Compute APEs and select only FFB and CPO estimates, as well as # obs. and clusters 
# Do not loop because each case is particular 
SM_REG <- 1
rm(ape_mat1, d_clean)
res_data <- res_data_list_ffbsm[[SM_REG]]
reg_res <- res_data[[1]]
d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
rm(res_data)
ape_mat1 <- make_APEs(reg_elm = SM_REG, K = 1) 
ape_mat1 <- ape_mat1 %>% as.matrix()
# add 2 rows to ape_mat1, equivalent to the estimate and CI for CPO control missing here.
ape_mat1 <- rbind(ape_mat1, matrix(ncol = ncol(ape_mat1), nrow = 2, data = ""))
# add 2 rows to indicate whether we control
ape_mat1 <- rbind(ape_mat1, matrix(ncol = ncol(ape_mat1), nrow = 2, data = ""))
# sort rows
ape_mat1 <- ape_mat1[c(1,2,5:8,3,4),]
ape_mat1

SM_REG <- 2
rm(ape_mat2, d_clean)
res_data <- res_data_list_ffbsm[[SM_REG]]
reg_res <- res_data[[1]]
d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
rm(res_data)
ape_mat2 <- make_APEs(reg_elm = SM_REG, K = 1, controls_pe=TRUE) 
ape_mat2 <- ape_mat2[c(1,2,3,4,7,8)] %>% as.matrix()
# add 2 rows to indicate whether we control
ape_mat2 <- rbind(ape_mat2, matrix(ncol = ncol(ape_mat2), nrow = 2, data = ""))
# sort rows
ape_mat2 <- ape_mat2[c(1:4,7:8,5:6),]
ape_mat2

SM_REG <- 3
rm(ape_mat3, d_clean)
res_data <- res_data_list_ffbsm[[SM_REG]]
reg_res <- res_data[[1]]
d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
rm(res_data)
ape_mat3 <- make_APEs(reg_elm = SM_REG, K = 1, controls_pe=TRUE) 
ape_mat3 <- ape_mat3[c(1,2,3,4,9,10)] %>% as.matrix()
# add 2 rows to indicate whether we control
ape_mat3 <- rbind(ape_mat3, matrix(ncol = ncol(ape_mat3), nrow = 2, data = c("", "Yes")))
# sort rows
ape_mat3 <- ape_mat3[c(1:4,7:8,5:6),]
ape_mat3

SM_REG <- 4
# order of variables differ in this case, due to where make_APEs computes interaction effects
rm(ape_mat4, d_clean)
res_data <- res_data_list_ffbsm[[SM_REG]]
reg_res <- res_data[[1]]
d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
rm(res_data)
ape_mat4 <- make_APEs(reg_elm = SM_REG, K = 1, controls_pe=TRUE) 
ape_mat4 <- ape_mat4[c(1,2,5,6,11,12)] %>% as.matrix()
# add 2 rows to indicate whether we control
ape_mat4 <- rbind(ape_mat4, matrix(ncol = ncol(ape_mat4), nrow = 2, data = c("Yes", "Yes")))
# sort rows
ape_mat4 <- ape_mat4[c(1:4,7:8,5:6),]
ape_mat4

ape_mat <- cbind(ape_mat1, ape_mat2, ape_mat3, ape_mat4)
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"),2),
                        "Industrial expansion proxy",
                        "... interacted with FFB price",
                        "Observations", 
                        "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Partial effects of FFB prices on deforestation for smallholder oil palm plantations") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Smallholder plantations" = 4),
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 4, 
            bold = TRUE)  %>%
  pack_rows("FFB price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("CPO price signal", 3, 4, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Additional controls", 5, 6, indent = FALSE, latex_gap_space = "1em",
            bold = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "15em",
              latex_valign = "m")# %>% 
  # column_spec(column = c(2:(ncol(ape_mat))),
  #             width = "7em",
  #             latex_valign = "c")




rm(res_data_list_ffbsm)

## Table 6. PRICE DYNAMICS ---------------------------------------------------------------------------------

## Regressions

res_data_list_prdyn <- list()
elm <- 1

size_list <- list("i", "sm", "unr", "a")
# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")
for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      # With SR only 
      # res_data_list_prdyn[[elm]] <- make_base_reg(island = ISL,
      #                                             outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
      #                                             illegal = ILL,
      #                                             only_sr = TRUE,
      #                                             offset = FALSE)
      # names(res_data_list_prdyn)[elm] <- paste0(ISL,"_",SIZE, "_",ILL,"_onlySR")
      # elm <- elm + 1
      # With SR AND MDR
      res_data_list_prdyn[[elm]] <- make_base_reg(island = "both",
                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                  illegal = ILL,
                                                  dynamics = TRUE,
                                                  interact_regressors = TRUE,
                                                  offset = FALSE)
      names(res_data_list_prdyn)[elm] <- paste0(ISL,"_",SIZE, "_",ILL,"_prdyn")
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    # With SR only 
    # res_data_list_prdyn[[elm]] <- make_base_reg(island = ISL,
    #                                             outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
    #                                             illegal = ILL,
    #                                             only_sr = TRUE,
    #                                             offset = FALSE)
    # names(res_data_list_prdyn)[elm] <- paste0(ISL,"_",SIZE, "_",ILL,"_onlySR")
    # elm <- elm + 1
    # With SR AND MDR
    res_data_list_prdyn[[elm]] <- make_base_reg(island = "both",
                                                outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                illegal = ILL,
                                                dynamics = TRUE,
                                                interact_regressors = TRUE,
                                                offset = FALSE)
    names(res_data_list_prdyn)[elm] <- paste0(ISL,"_",SIZE, "_",ILL,"_prdyn")
    elm <- elm + 1
  }}
names(res_data_list_prdyn)

## PARTIAL EFFECTS
# With SR only 
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_prdyn)){ #c(1,3,5,7,9,11)
  # # get estimation results and data; this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  # res_data <- res_data_list_prdyn[[REGELM]]
  # reg_res <- res_data[[1]]
  # d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  # rm(res_data)
  # temp_mat <- make_APEs(reg_elm = REGELM, K = 1) # and for the same reason, this cannot be wrapped in other functions
  # # add 4 lines to temp_mat (equivalent to estimate and CI of MR and interaction)
  # temp_mat <- rbind(temp_mat, matrix(ncol = ncol(temp_mat), nrow = 4))
  # temp_mat <- temp_mat[c(1,2,5,6,7,8,3,4),] # rearrange rows
  # temp_mat[is.na(temp_mat)] <- "" 
  # ape_mat[[REGELM]] <- temp_mat
  # reg_res_list[[REGELM]] <- reg_res
  # rm(d_clean, reg_res)
  # REGELM <- REGELM + 1
  
  # With SR AND MDR
  res_data <- res_data_list_prdyn[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM, K = 2, rounding = 3) 
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Partial effects of short- and medium-run price signals on deforestation across Indonesian oil palm plantations") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1#,
                     # "(7)" = 1,
                     # "(8)" = 1, 
                     # "(9)" = 1, 
                     # "(10)" = 1,
                     # "(11)" = 1, 
                     # "(12)" = 1
                     ),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 4, 
            bold = TRUE)  %>%
  pack_rows("Short-run price signal", 1, 2, indent = FALSE, # short run is indeed always before MR in the regressor vector construction in make_base_reg
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Medium-run price signal", 3, 4, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Interaction", 5, 6, 
            bold = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "12em",
              latex_valign = "m") 


rm(res_data_list_prdyn)


# DESCRIPTIVE STATISTICS ------------------------------------------------------------------------
## Table 1. Des.stats. by plantation type ---------------------------------

## Descriptive statistics for different types of plantations
# template to store 
list_desstat_all <- list()
i <- 1
for(ILL in ill_status){
  list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[[paste0("both_i_",ILL)]][[2]])
  i <- i +1
}
# Add smallholders
list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[["both_sm_all"]][[2]])
i <- i +1

# Add unregulated
list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[["both_unr_all"]][[2]])
i <- i +1

# Add all
list_desstat_all[[i]] <- make_desstats_simple(sample_1 = res_data_list_full[["both_a_all"]][[2]])
i <- i +1


des_table <- bind_cols(list_desstat_all) %>% as.matrix()
# row names
row.names(des_table) <- c("Deforestation (ha)",
                          "CPO price signal (USD/ton)", 
                          "FFB price signal (USD/ton)", 
                          "# reachable UML mills", 
                          "Primary forest cover 2000 (ha)",
                          "Secondary forest cover 2000 (ha)",
                          "Secondary forest loss for industrial",
                          "Secondary forest loss for smallholders"
                          # "Public ownership (%)", 
                          # "Domestic private ownership (%)", 
                          # "Foreign ownership (%)", 
                          #"Competition", 
                          # "# reachable mills"
                          )
# "n_reachable_uml",
# "pfc2000_total_pixelcount", 
# "fc2000_30th_pixelcount", 
# "lucfip_pixelcount",
# "lucfsmp_pixelcount"
des_table
colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Estimation samples - descriptive statistics") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean (SD)" = 1, "med. [min;max]" = 1,
                     "mean (SD)" = 1, "med. [min;max]" = 1,
                     "mean (SD)" = 1, "med. [min;max]" = 1,
                     "mean (SD)" = 1, "med. [min;max]" = 1,
                     "mean (SD)" = 1, "med. [min;max]" = 1,
                     "mean (SD)" = 1, "med. [min;max]" = 1),
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "# grid cells = 3983 \n # grid cell-year = 24131" = 2,
                     "# grid cells = 3189 \n # grid cell-year = 17091" = 2, 
                     "# grid cells = 11782 \n # grid cell-year = 65368" = 2, 
                     "# grid cells = 3211 \n # grid cell-year = 20721" = 2, 
                     "# grid cells = 4266 \n # grid cell-year = 24596" = 2, 
                     "# grid cells = 12687 \n # grid cell-year = 71926" = 2),
                   align = "c",
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Legal" = 2,
                     "Illegal" = 2, 
                     "All" = 2,
                     " " = 6), 
                   bold = FALSE,
                   align = "c",
                   strikeout = F) %>% 
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 6,
                     "Smallholder plantations" = 2, 
                     "Unregulated plantations" = 2, 
                     "All" = 2), 
                   bold = FALSE,
                   align = "c",
                   strikeout = F) %>% 
  row_spec(1:2, extra_latex_after = "\\\\") 


## Table A.1. Accumulated deforestation in different catchment radius / sample ------------------------------------------------------------------------

### TOTAL DEFORESTATION 
rows <- c("Sumatra", "Kalimantan", "both")
cols <- c("sample", "30km", "50km", "total")

accu_lucpfp <- matrix(ncol = length(cols), nrow = length(rows)) # 4 cols for sample, 30km, 50km and total. 3 rows for Sumatra, Kalimantan, and both. 
row.names(accu_lucpfp) <- rows
colnames(accu_lucpfp) <- cols



#island <- "Sumatra"
for(island in c("Sumatra", "Kalimantan")){
  # LUCPFIP
  defo_dyna <- list()
  for(DYNA in c("rapid", "slow")){
    brick_lucpfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",DYNA,"_",island,"_",parcel_size/1000,"km_total.tif")))
    
    # select 2002-2014 layers 
    layer_names <- paste0("parcel_lucpfip_",DYNA,"_",island,"_",parcel_size/1000,"km_total.",c(2:14))
    
    brick_lucpfip <- raster::subset(brick_lucpfip, layer_names)
    
    # Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
    r_accu_lucfp <- calc(brick_lucpfip, fun = sum, na.rm = TRUE)
    
    defo_dyna[[match(DYNA, c("rapid", "slow"))]] <- raster::extract(x = r_accu_lucfp, 
                                                                    y = island_sf_prj[island_sf_prj$shape_des==island,] %>% st_geometry() %>% as("Spatial"), 
                                                                    fun = sum, 
                                                                    na.rm = TRUE) 
  }
  accu_lucpfp[island, "total"] <- sum(unlist(defo_dyna))
  
  # LUCPSMP
  defo_sm <- list()
  for(size in c("s", "m")){
    brick_lucpfsmp <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpf",size,"p_",island,"_",parcel_size/1000,"km_total.tif")))
    
    # select 2002-2014 layers 
    layer_names <- paste0("parcel_lucpf",size,"p_",island,"_",parcel_size/1000,"km_total.",c(2:14))
    
    brick_lucpfsmp <- raster::subset(brick_lucpfsmp, layer_names)
    
    # Add up annual aggregated LUCFP (result is a single layer with cell values = the sum of annual cell values over the selected time period)
    r_accu_lucfp <- calc(brick_lucpfsmp, fun = sum, na.rm = TRUE)
    
    defo_sm[[match(size, c("s", "m"))]] <- raster::extract(x = r_accu_lucfp, 
                                                           y = island_sf_prj[island_sf_prj$shape_des==island,] %>% st_geometry() %>% as("Spatial"), 
                                                           fun = sum, 
                                                           na.rm = TRUE) 
  }
  
  # Add indus and smallholder deforestation
  accu_lucpfp[island, "total"] <- sum(unlist(defo_dyna)) + sum(unlist(defo_sm)) 
}

### Now restrict to full CR. 

# sum pixelcounts over parcels and years, under the condition that n_reachable is not null
accu_lucpfp["Sumatra", "30km"] <- dplyr::filter(d_30_suma, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Sumatra", "50km"] <- dplyr::filter(d_50_suma, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Kalimantan", "30km"] <- dplyr::filter(d_30_kali, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()
accu_lucpfp["Kalimantan", "50km"] <- dplyr::filter(d_50_kali, n_reachable_ibs > 0 & year > 2001 & year < 2015)$lucpfap_pixelcount %>% sum()

### FOr analysis sample
for(ISL in c("Sumatra", "Kalimantan")){
  res_data_list  <- make_base_reg(island = ISL,
                                  outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                  output_full = FALSE)
  
  d_clean <- res_data_list[[2]]
  accu_lucpfp[ISL, "sample"] <- sum(d_clean$lucpfap_pixelcount)
}


# Fill the row of total over islands
for(col in c(1:ncol(accu_lucpfp))){
  accu_lucpfp["both", col] <- sum(accu_lucpfp["Sumatra", col], 
                                  accu_lucpfp["Kalimantan", col])
}

accu_lucpfp <- accu_lucpfp*pixel_area_ha/1000 # convert to *kilo* hectares
accu_lucpfp <- accu_lucpfp %>% round(digits = 2)

row.names(accu_lucpfp)[3] <- "Both"

accu_lucpfp_save <- accu_lucpfp

# Look at how much land is concerned with immediate and transitional deforestation
res_data_list  <- make_base_reg(island = "both",
                                #end_year = 2010,
                                outcome_variable = paste0("lucpfip_rapid_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                output_full = FALSE)

d_clean <- res_data_list[[2]]
sum(d_clean$lucpfip_rapid_pixelcount)*pixel_area_ha/1000 %>% round(digits = 2)

res_data_list  <- make_base_reg(island = "both",
                                end_year = 2010,
                                outcome_variable = paste0("lucpfip_slow_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                output_full = FALSE)

d_clean <- res_data_list[[2]]
sum(d_clean$lucpfip_slow_pixelcount)*pixel_area_ha/1000 %>% round(digits = 2)
rm(d_clean)
rm(res_data_list)
# Finally, add figures from Austin et al. 2017 (SI) - no because not comparable enough because we add smallholders without breaking down here. 
# austin <- matrix(nrow = 2, ncol = ncol(accu_lucpfp))
# austin["Sumatra", 4] <- "(448)"
# austin["Kalimantan", 4] <- "(1021)"

### Print the LateX table code ###

accu_lucpfp
colnames(accu_lucpfp) <- NULL

options(knitr.table.format = "latex") 
kable(accu_lucpfp, booktabs = T, align = "c", 
      caption = "Deforestation accumulated over 2002-2014, in kha.") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "Sample" = 1, 
                     "30km from sample mill" = 1, 
                     "50km from sample mill" = 1,
                     "Total" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  column_spec(column = 1,
              width = "5em") 

# %>% column_spec(column = c(2:ncol(accu_lucpfp)),
#             width = "5em") 

rm(res_data_list_des)

### INDUSTRIAL VS SMALLHOLDERS
# res_data_list_des <- list()
# for(SIZE in c("i", "sm")){
#   res_data_list_des[[paste0("both_",SIZE)]] <- make_base_reg(island = "both",
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                                                  output_full = FALSE)
# }
# # export for analysis on other plateform
# d_clean_sf <- st_as_sf(d_clean, coords = c("lon","lat"), crs = 4326)
# st_write(d_clean_sf, file.path("temp_data/actual_analysis_data_Kalimantan_a"), driver = "ESRI Shapefile", delete_dsn = TRUE)



# des_table <- make_desstats(sample_1 = res_data_list_des[["both_i"]][[2]],
#                            sample_2 = res_data_list_des[["both_sm"]][[2]])
# colnames(des_table) <- NULL
# 
# options(knitr.table.format = "latex") 
# kable(des_table, booktabs = T, align = "c", 
#       caption = "Estimation sample - descriptive statistics") %>% 
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
#   add_header_above(c(" " = 1, 
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
#                      "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
#                      "p-value" = 1, "p-value" = 1),
#                    align = "c", 
#                    strikeout = F) %>% 
#   add_header_above(c(" " = 1, 
#                      "# grid cells = 4757 \n # grid cell-year = 31721" = 3,
#                      "# grid cells = 8368 \n # grid cell-year = 91620" = 3, 
#                      " " = 2),
#                    align = "c",
#                    strikeout = F) %>% 
#   add_header_above(c(" " = 1,
#                      "Without missing values" = 3,
#                      "With missing values" = 3, 
#                      "t test" = 1,
#                      "KS test" = 1), 
#                    bold = TRUE,
#                    align = "c",
#                    strikeout = F) %>% 
#   # pack_rows("Ownership signals (%)", 3, 5, 
#   #           italic = TRUE)  %>%
#   # column_spec(column = c(2,3,5,6,8,9),
#   #             width = "3em") %>%
#   # column_spec(column = c(4,7,10),
#   #             width = "9em") %>%
#   
#   footnote(general = c("Note"),
#            threeparttable = TRUE, 
#            escape = TRUE) 




## Table A.2 IBS -----------------------------------------------------------

# all IBS oil palm related establishments. 
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))

length(unique(ibs$firm_id))
length(unique(ibs$firm_id[ibs$uml_matched_sample==1])) 
length(unique(ibs$firm_id[ibs$is_mill==1])) 
length(unique(ibs$firm_id[ibs$analysis_sample==1]))
# ibs[ibs$geo_sample & !ibs$is_mill,1:40]
# out of 1473 establishments initially extracted from IBS,
# 1004 are involved at least one year in some FFB processing or CPO or PKO producing, but 930 are not in Java nor Bali 
# out of which 468 have been geolocalized (and 2 more are in Java)
# 2 additional IBS firms have been matched with UML and geolocalized but have no sign of FFB, PKO or CPO activity. 
# (firm_id 2036 and 55630)


# Now we also remove those that have some FFB-CPO-PKO activity but have been identified as refineries. 
# But why would we remove refineries? They may be different but if they input FFB, they have an impact on proximate LUC. 
# Because we excluded refineries from geo_sample, and here the purpose is to compare this sample to the population of mills. 
# let's first see the comparative stat des without removing refineries. 


#### Print the LateX table code 
des_table <- make_des_table_ibs(ibs) 
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs[ibs$analysis_sample==TRUE, "firm_id"]))
length(unique(ibs[ibs$is_mill==TRUE, "firm_id"]))
des_table
colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, all islands") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "p-value" = 1,
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 587 mills" = 3,
                     "All IBS palm oil mills \n n = 930 mills" = 3, 
                     "t-test" = 1,
                     "KS test" = 1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "4em") %>% 
  column_spec(column = c(4,7),
              width = "9em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 


## Table A.3 Des.stats. with/without NAs  ------------------------------------------------------------------------

### OVERALL 
res_data_list_des <- list()
res_data_list_des[["both_a"]] <- make_base_reg(island = "both",
                                               outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                               output_full = TRUE)

# # export for analysis on other plateform
# d_clean_sf <- st_as_sf(d_clean, coords = c("lon","lat"), crs = 4326)
# st_write(d_clean_sf, file.path("temp_data/actual_analysis_data_Kalimantan_a"), driver = "ESRI Shapefile", delete_dsn = TRUE)



des_table <- make_desstats(sample_1 = res_data_list_des[["both_a"]][[2]],
                           sample_2 = res_data_list_des[["both_a"]][[3]])
des_table
colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Estimation sample - descriptive statistics") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "p-value" = 1, "p-value" = 1),
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1, 
                     "# grid cells = 12687 \n # grid cell-year = 71926" = 3,
                     "# grid cells = 22570 \n # grid cell-year = 215667" = 3, 
                     " " = 2),
                   align = "c",
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Without missing values" = 3,
                     "With missing values" = 3, 
                     "t test" = 1,
                     "KS test" = 1), 
                   bold = TRUE,
                   align = "c",
                   strikeout = F) 
# pack_rows("Ownership signals (%)", 3, 5, 
#           italic = TRUE)  %>%
# column_spec(column = c(2,3,5,6,8,9),
#             width = "3em") %>%
# column_spec(column = c(4,7,10),
#             width = "9em") %>%




## Figure B.3 Price & deforestation time series -------------------------------------------------------

# Price time series
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# Restrict to analysis sample 
as_cs = 
  ibs %>% 
  filter(analysis_sample==1) 
price_ts = 
  rbind(
    summarise(as_cs, .by = "year", 
              value = mean(cif_rtdm_cpo, na.rm = TRUE)) %>%
      mutate(variable = "CIF Rotterdam CPO"), 
    summarise(as_cs, .by = "year", 
              value = mean(cpo_price_imp1, na.rm = TRUE)) %>%
    mutate(variable = "Mill-gate CPO"), 
    summarise(as_cs, .by = "year", 
              value = mean(ffb_price_imp1, na.rm = TRUE)) %>%
    mutate(variable = "Mill-gate FFB")) %>% 
  mutate(type = "Price (USD/ton)")

  
# Deforestation (kha) time series
defo_ts_list <- list()
for(i in 1:length(res_data_list_full)){
  sample_1 = res_data_list_full[[i]][[2]]
  dep_var_1 <- names(sample_1)[grepl("pixelcount",names(sample_1))]
  
  defo_ts_list[[i]] = 
    sample_1 %>% 
    summarise(.by = "year", 
              value := sum(!!as.symbol(dep_var_1), na.rm = TRUE) / 1000) %>% # (convert to kha)
    mutate(
      variable = names(res_data_list_full[i]),
      type = "Deforestation (kha)")
}
# Stack and rename
df = 
  defo_ts_list %>% 
  bind_rows() %>% 
  rbind(price_ts) %>% 
  filter(year < 2015) %>% 
  mutate(variable = case_when(
    variable == "both_i_no_ill2" ~ "Legal industrial plantations",
    variable == "both_i_ill2" ~ "Illegal industrial plantations",
    variable == "both_i_all" ~ "All industrial plantations",
    variable == "both_sm_all" ~ "Smallholder plantations",
    variable == "both_unr_all" ~ "Unregulated plantations",
    variable == "both_a_all" ~ "All plantations",
    TRUE ~ variable
  ), 
  variable = factor(variable, levels=c("CIF Rotterdam CPO", 
                                       "Mill-gate CPO", 
                                       "Mill-gate FFB",
                                       "Legal industrial plantations",
                                       "Illegal industrial plantations", 
                                       "All industrial plantations",
                                       "Smallholder plantations",
                                       "Unregulated plantations", 
                                       "All plantations")))
  

# Rescale factor for Deforestation (kha) to match Price (USD/ton) scale (adjust as needed)
scale_factor <- max(df[df$type=="Price (USD/ton)", "value"], na.rm = T) / max(df[df$type=="Deforestation (kha)", "value"])

df <- df %>% mutate(value_rescaled = ifelse(type == "Deforestation (kha)", value * scale_factor, value))

# Plot
ts_plot <-
  ggplot() +
  # Prices
  geom_line(
    data = df %>% filter(type == "Price (USD/ton)"),
    aes(x = year, y = value_rescaled, color = variable, linetype = variable), linewidth = 1.1
  ) +
  # Deforestation (kha) (rescaled)
  geom_line(
    data = df %>% filter(type == "Deforestation (kha)"),
    aes(x = year, y = value_rescaled, color = variable, linetype = variable), linewidth = 1.1
  ) +
  scale_y_continuous(
    name = "Price (USD/ton)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Deforestation (kha)")
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "CIF Rotterdam CPO" = "black",
      "Mill-gate CPO" = "pink", 
      "Mill-gate FFB" = "orange",
      "Legal industrial plantations" = "red",
      "Illegal industrial plantations" = "green",
      "All industrial plantations" = "purple",
      "Smallholder plantations" = "brown",
      "Unregulated plantations" = "blue",
      "All plantations" = "darkgreen"
    )
  ) +
  scale_linetype_manual(
    name = "",
    values = c(
      "CIF Rotterdam CPO" = "solid",
      "Mill-gate CPO" = "solid", 
      "Mill-gate FFB" = "solid",
      "Legal industrial plantations" = "dotted",
      "Illegal industrial plantations" = "dotdash",
      "All industrial plantations" = "longdash",
      "Smallholder plantations" = "twodash",
      "Unregulated plantations" = "dashed",
      "All plantations" = "dashed"
      )
  ) +
  guides(
    color = guide_legend(title = ""),
    linetype = guide_legend(title = "")
  ) +
  scale_x_continuous(breaks = seq(1998, 2014, by = 2)) +
  theme_minimal() +
  labs(x = "Year") +
  theme(
    axis.text =  element_text(size = 14),
    axis.title =  element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )

ts_plot

ggsave(ts_plot, 
       filename = "ts_plot.png", 
       width=20, height=10)


## Price variation within district-year --------------------------------------------------------------
# Restrict to analysis sample 
as <- ibs[ibs$analysis_sample==TRUE,]

# remove no variation 
rm_fevar <- fixest::feols(fml = as.formula("cpo_price_imp1 ~ 1"),
                          data = as)
# and take the standard deviation of remaining variation
sd(rm_fevar$residuals)
# which is exactly equal to:
sd(as$cpo_price_imp1, na.rm= T)

# remove year variation 
rm_fevar <- fixest::feols(fml = as.formula("cpo_price_imp1 ~ 1 | year"),
                          data = as)
# and take the standard deviation of remaining variation (which is also the RMSE of this reg)
sd(rm_fevar$residuals) # That is 20% relative to average price. 

# remove district-year variation
rm_fevar <- fixest::feols(fml = as.formula("cpo_price_imp1 ~ 1 | district^year"),
                          data = as)
# and take the standard deviation of remaining variation (which is also the RMSE of this reg)
sd(rm_fevar$residuals) # That is 20% relative to average price. 


# Repeat in FFB prices 
# Removing no variation, 
sd(as$ffb_price_imp1, na.rm = T)

# Removing year variation, 
rm_fevar <- fixest::feols(fml = as.formula("ffb_price_imp1 ~ 1 | year"),
                                     data = as)
# and take the standard deviation of remaining variation
sd(rm_fevar$residuals)

# Removing province level variation, 
# which is the level at which prices are supposed to be collectively determined 
rm_fevar <- fixest::feols(fml = as.formula("ffb_price_imp1 ~ 1 | province^year"),
                              data = as)
# and take the standard deviation of remaining variation
sd(rm_fevar$residuals)

# check district-year
rm_fevar_ffb <- fixest::feols(fml = as.formula("ffb_price_imp1 ~ 1 | district^year"),
                              data = as)
# same order of magnitude relative to average price (22%)
sd(rm_fevar_ffb$residuals)

rm(as)



## Sets of reachable mills -----------------
d_clean <-  res_data_list_full[["both_a_all"]][[2]]

# How many they are
d_clean$reachable %>% unique() %>% length() # 1441

# What is the average area and duration of a set of reachable mills 
reachset = 
  d_clean %>% 
  summarise(.by = "reachable", 
            n_obs = n(),
            n_cells = length(unique(lonlat)), 
            n_years = length(unique(year)))

reachset$n_obs %>% summary()
reachset$n_cells %>% summary()
reachset$n_years %>% summary()


# MAPS -----------------------------------------------------------

## Figure B.1 Sample cells and mills 
# prepare backgroud layers with other countries
countries <- st_read(file.path("input_data/Global_LSIB_Polygons_Detailed"))
countries <- countries[countries$COUNTRY_NA == "Indonesia" | 
                         countries$COUNTRY_NA == "Malaysia" | 
                         countries$COUNTRY_NA == "Thailand" | 
                         countries$COUNTRY_NA == "Singapore" | 
                         countries$COUNTRY_NA == "Brunei", ]
# these two lines to speed up mapping
countries <- st_transform(countries, crs = indonesian_crs) %>% st_simplify(dTolerance = 1000)
countries <- st_transform(countries, crs = 4326)
indonesia_sf <- countries %>% filter(COUNTRY_NA == "Indonesia")

res_data_both_a_all <- make_base_reg(island = "both",
                                     outcome_variable = paste0("lucpfap_pixelcount"), 
                                     illegal = "all",
                                     output_full = FALSE)

d_clean <- res_data_both_a_all[[2]]
d_clean_cs <- d_clean[!duplicated(d_clean$lonlat),]

d_cs <- st_as_sf(d_clean_cs, coords = c("lon", "lat"), crs = 4326)
d_cs <- st_transform(d_cs, crs = indonesian_crs)

d_cs <- st_buffer(d_cs, dist = 1600)
st_geometry(d_cs) <- sapply(st_geometry(d_cs), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
d_geo <- st_union(st_geometry(d_cs))
d_geo <- st_transform(d_geo, crs = 4326)
d_geo <- d_geo %>% st_as_sf() %>% mutate(label = "Plantation sites")


# d_cs <- ddply(d_clean, "lonlat", summarise, 
#               accu_lucfp = sum(lucpfap_pixelcount))
# 
# d_cs <- left_join(d_cs, d_clean_cs[,c("lonlat", "lon", "lat")], by = "lonlat")

### MILLS
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# keep only those that are used in analysis. 
ibs <- ibs[ibs$analysis_sample ==1,]
# make it a cross section
ibs <- ibs[!duplicated(ibs$firm_id),]
ibs <- ibs[!is.na(ibs$lat),]
ibs <- st_as_sf(ibs, coords = c("lon", "lat"), remove = FALSE, crs = 4326)

uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
# count the number of UML mills
uml$uml_id %>% unique() %>% length()
uml$trase_code %>% unique() %>% length()

# we do not plot uml actually, it's too much. 
# uml <- uml[!is.na(uml$lat),]
# uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = FALSE, crs = 4326)
# # remove those matched with ibs
# uml <- uml[!(uml$trase_code %in% ibs$trase_code),]

legend_df <- data.frame(x = 0, y = 0, label = "Mills") %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326)

map1 =
  ggplot() +
    geom_sf(data = countries, fill = "grey", col = "black") +
    geom_sf(data = indonesia_sf, fill = "white", col = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "lightgrey"), 
          legend.position = "bottom",
          legend.text = element_text(size = 18, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.key = element_blank(), 
          axis.text = element_text(size = 10),
          # legend.key.spacing = unit(2, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.key.width = unit(1.5, "cm")
          # plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
          ) +
    # geom_sf(data=st_geometry(d_geo), color=alpha("grey",0.2))+
    geom_sf(data= d_geo, aes(fill = label, alpha=label), color = "transparent") +
    scale_fill_manual(
      name = " ",
      values = c("Plantation sites" = "blue")
    ) +
    scale_alpha_manual(
      name = " ",
      values = c("Plantation sites" = 0.2)
    ) +
    geom_sf(data = legend_df, aes(shape = label, size = label, color = label), fill = NA) +
    scale_size_manual(
      name = " ",
      values = c("Mills" = 4)
    ) +
    scale_shape_manual(
      name = " ",
      values = c("Mills" = 16)
    ) +
    scale_color_manual(
      name = " ",
      values = c("Mills" = "red")
    ) +
    geom_sf(data = ibs, color = "red", size = 0.15) + 
    # geom_sf(data = st_geometry(countries), fill = "transparent") +
    coord_sf(xlim = c(94, 119), ylim = c(-7, 7), expand = FALSE)
    

ggsave(map1, 
       filename = "map_mills_sites.png", 
       width=18, height=11)

# leaflet() %>% 
#   # we do not add these tiles if the figure is purposed for print in white and black. 
#   # addTiles()%>%
#   # addProviderTiles(providers$Esri.WorldGrayCanvas) %>% # CartoDB.PositronNoLabels 
#   addPolygons(data = countries, stroke = TRUE, color = "black", weight = 1,
#               fill = FALSE) %>%
#   setView(lat = -0.493, 
#           lng = 107.367, 
#           zoom = 4.5) %>% 
#   #addPolygons(data = island_sf, stroke = FALSE, fill = FALSE, fillColor = "grey", fillOpacity = 0) %>% 
#   addPolygons(data = st_geometry(d_geo), 
#               stroke = TRUE, opacity = 0.25,  weight = 2, color = "green") %>%   #~cb_ind(d_cs$accu_lucfp)
#   addCircleMarkers(data = ibs, radius = 0.001, fillOpacity = 1, fillColor = "black", stroke = FALSE, weight = 0) 
# # exported in width = 1150 and height = 560 and zoom once 

rm(d_clean_cs, d_cs)

##### Figure B.2 Accumulated deforestation #####  

# Display outcome vars in grid cells within 82km of a known (UML) mill. 
# This one for industrial has legal/illegal info (dirty naming but its UML CR, see add_CR_parcel_variables.R)
lucpfip <- readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_land_des_",
                                    parcel_size/1000,"km_",
                                    "82CR.rds"))) 
# Use lucpfip_pixelcount_total, not a problem since this is just for plotting and ~the same as rapid + slow. 
all_indus = 
  lucpfip %>% 
  filter(year >= 2002 & year <= 2014) %>% 
  mutate(lucpfip_pixelcount = lucpfip_pixelcount_total) %>% 
  summarise(.by = lonlat, 
            lon = unique(lon),
            lat = unique(lat),
            # Keep legal var here
            illegal2 = unique(illegal2),
            # the rounding to 0 decimales on percentage points means that cells with less 
            # than 0.5% (4.5ha) of deforestation are rounded to 0 and then converted to NA (transparent).  
            accu_LT_cellpct = round(100*sum(lucpfip_pixelcount)*pixel_area_ha/900, 0)) %>% 
  mutate(loss_type = "all_indus")

# Legal industrial
leg_indus = 
  all_indus %>% 
  filter(!is.na(illegal2) & !illegal2) %>% 
  mutate(loss_type = "leg_indus")

# Illegal industrial
ill_indus = 
  all_indus %>% 
  filter(!is.na(illegal2) & illegal2) %>% 
  mutate(loss_type = "ill_indus")

# Smallholders
lucpfsmp <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfsmp_panel_",
                                     parcel_size/1000,"km_",
                                     "82km_UML_CR.rds")))
all_sm = 
  lucpfsmp %>% 
  filter(year >= 2002 & year <= 2014) %>% 
  mutate(lucpfsmp_pixelcount = lucpfsp_pixelcount_total + lucpfmp_pixelcount_total) %>% 
  summarise(.by = lonlat, 
            lon = unique(lon),
            lat = unique(lat), 
            accu_LT_cellpct = round(100*sum(lucpfsmp_pixelcount)*pixel_area_ha/900, 0)) %>% 
  mutate(loss_type = "all_sm")

# unregulated
unreg = 
  full_join(ill_indus, 
            all_sm, 
            by = c("lonlat", "lon", "lat")) %>% 
  mutate(accu_LT_cellpct = accu_LT_cellpct.x + accu_LT_cellpct.y) %>% 
  dplyr::select(!contains(".")) %>% 
  mutate(loss_type = "unreg")

# all 
all = 
  full_join(all_indus, 
            all_sm, 
            by = c("lonlat", "lon", "lat")) %>% 
  mutate(accu_LT_cellpct = accu_LT_cellpct.x + accu_LT_cellpct.y) %>% 
  dplyr::select(!contains(".")) %>% 
  mutate(loss_type = "all")

d_stack = 
  rbind(leg_indus %>% dplyr::select(names(all_sm)), 
        ill_indus %>% dplyr::select(names(all_sm)), 
        all_indus %>% dplyr::select(names(all_sm)), 
        all_sm,
        unreg %>% dplyr::select(names(all_sm)),
        all %>% dplyr::select(names(all_sm))) 

d_stack$accu_LT_cellpct <- na_if(d_stack$accu_LT_cellpct, 0)

# Spatialize to 3km grid cells
d_stack <- st_as_sf(d_stack, coords = c("lon", "lat"), crs = 4326)
d_stack <- st_transform(d_stack, crs = indonesian_crs)
d_stack <- st_buffer(d_stack, dist = 1600) # half the size of a cell + TAKING SOME MARGIN TO PREVENT WHITE LINES BETWEEN GRIDS nearer to the equator
st_geometry(d_stack) <- sapply(st_geometry(d_stack), FUN = function(x){st_as_sfc(st_bbox(x))}) %>% st_sfc(crs = indonesian_crs)
d_stack <- st_transform(d_stack, crs = 4326)

d_stack %>% filter(loss_type == "all_indus") %>% pull(accu_LT_cellpct) %>% summary()
d_stack %>% filter(loss_type == "all_sm") %>% pull(accu_LT_cellpct) %>% summary()
d_stack %>% filter(loss_type == "all") %>% pull(accu_LT_cellpct) %>% summary()

d_stack = 
  d_stack %>% 
  mutate(
    loss_type = case_when(
      loss_type=="leg_indus"   ~ "(1) Legal industrial plantations", 
      loss_type=="ill_indus"   ~ "(2) Illegal industrial plantations", 
      loss_type=="all_indus"   ~ "(3) All industrial plantations", 
      loss_type=="all_sm"      ~ "(4) Smallholder plantations", 
      loss_type=="unreg"       ~ "(5) Unregulated plantations",
      loss_type=="all"         ~ "(6) All plantations"
    ), 
    loss_type = factor(loss_type, levels=c("(1) Legal industrial plantations",
                                           "(2) Illegal industrial plantations", 
                                           "(3) All industrial plantations",
                                           "(4) Smallholder plantations",
                                           "(5) Unregulated plantations", 
                                           "(6) All plantations"))
  ) 

ibs_cr_sf = 
  rbind(
    ibs %>%
      filter(island_name == "Sumatra") %>% 
      st_transform(crs = indonesian_crs) %>%
      st_buffer(dist = 30*1e3) %>% 
      st_transform(crs = 4326),
    ibs %>%
      filter(island_name == "Kalimantan") %>% 
      st_transform(crs = indonesian_crs) %>%
      st_buffer(dist = 50*1e3) %>% 
      st_transform(crs = 4326)
  ) %>% 
  st_geometry() %>% 
  st_union()

# ibs_cr_sf <- ibs_cr_sf %>% st_as_sf() %>% mutate(label = label_cr)
legend_df <- data.frame(x = 0, y = 0, label = "Catchment radius \nof mills in analysis") %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326)


map <- 
  ggplot() +
  geom_sf(data = countries, fill = "grey", col = "black") +
  geom_sf(data = indonesia_sf, fill = "white", col = "black") +
  geom_sf(data = d_stack, aes(fill = accu_LT_cellpct), col = "transparent") + # , lwd = NA this would prevent the key in the legend to appear in case of factors 
  scale_fill_viridis(name = "% of 900 ha cell",
                     discrete = FALSE, # TRUE,
                     option="A",
                     direction = -1, 
                     na.value = "transparent") +
  facet_wrap(facets = ~loss_type, ncol = 2, nrow = 3, 
             strip.position = "top") +
  
  geom_sf(data=ibs_cr_sf , color = "darkgrey", fill = NA, linewidth = 0.5, show.legend = FALSE) +
  # Dummy point layer to force a circular legend key
  geom_sf(data = legend_df, aes(shape = label, size = label, color = label), fill = NA) +
  scale_size_manual(
    name = " ",
    values = c("Catchment radius \nof mills in analysis" = 9)
  ) +
  scale_shape_manual(
    name = " ",
    values = c("Catchment radius \nof mills in analysis" = 21)
    ) +
  # scale_discrete_manual(aesthetic = "linewidth", values = c("Catchment radius \nof mills in analysis" = 3)) +
  scale_color_manual(
    name = " ",
    values = c("Catchment radius \nof mills in analysis" = "darkgrey")
    ) +
  # Zoom in 
  coord_sf(xlim = c(94, 118.5), ylim = c(-4, 5), expand = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "lightgrey"), 
        legend.position = "bottom",
        legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key = element_blank(), 
        axis.text = element_text(size = 10),
        # legend.key.spacing = unit(2, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(2.2, "cm"),
        strip.text = element_text(size = 15, face = "bold"), 
        strip.background = element_rect(fill = "white", colour = "lightgrey")) 


ggsave(map, 
       filename = "map_accu_inCR_bytype.png", 
       width=20, height=13)

# map


#### DESCRIPTIVE PARTIAL AUTOCORRELATION FUNCTION OF PRICES 
# # read this panel, as it still features the annual price observations from 1998 (final data only from 2001)
# # 50km CR because it is the most general
# RHS_50 <-  readRDS(file.path(paste0("temp_data/processed_parcels/parcels_panel_w_dyn_",
#                                     parcel_size/1000,"km_",
#                                     "50CR.rds")))
# 
# prices <- RHS_50[,c("lonlat", "year", "cpo_price_imp1")]
# 
# nrow(prices)/length(unique(prices$year))
# 
# # # keep only series with non-missing values
# prices_nona <-
#   prices %>%
#   filter(year >= 2000 & year <= 2015) %>%
#   group_by(lonlat) %>%
#   mutate(never_missing = !anyNA(cpo_price_imp1)) %>%
#   ungroup() %>%
#   filter(never_missing) %>%
#   dplyr::select(-never_missing)
# prices_nona$lonlat %>% unique() %>% length()
# 
# # keep only one
# sample_ids <- prices_nona$lonlat %>% sample(size = 1) # 700 this is roughly 1% of the grid cells. 
# prices_nona <- prices_nona[prices_nona$lonlat %in% sample_ids,]
# 
# prices_ts <- 
#   prices_nona %>% 
#   dplyr::arrange(lonlat, year) %>% 
#   pull(cpo_price_imp1) %>% 
#   ts()
# 
# length_panel <- length(prices_ts)
# ur.kpss(prices_ts) %>% summary() # --> data need to be differenced
# ur.kpss(diff(prices_ts, differences = 1)) %>% summary() # --> first differencing suffices
# ur.kpss(diff(prices_ts, differences = 5)) %>% summary() # --> first differencing suffices
# store <- list()
# for(ar in 1:11){
#   res <- stats::arima(prices_ts, order = c(ar,1,0), method = "ML") 
#   store[[ar]] <- BIC(res)
# }
# 
# # opt_ar <- which(max(unlist(store)))
# BIC(unlist(store))
# 
# arima_res <- 
# arima_res
# stargazer(arima_res, font.size = "footnotesize")
# 
# 
# 
# # Apply descriptive statistics to a sample of plantations, because neighboring plantations can have spatially correlated price time series, 
# # that would conflate the inference performed here.
# sample_ids <- prices$lonlat %>% sample(size = 700, replace = FALSE) # 700 this is roughly 1% of the grid cells. 
# prices <- prices[prices$lonlat %in% sample_ids,]
# 
# # we need it to be balanced
# nrow(prices) == length(unique(prices$lonlat)) * length(unique(prices$year))
# length(unique(prices$lonlat))
# 
# length_panel <- length(unique(prices$year))
# 
# prices_extd <- dplyr::arrange(prices, lonlat, year)
# # prices_ts <- matrix(prices_extd[,c("cpo_price_imp1")], 
# #                     nrow = length_panel, 
# #                     byrow = FALSE) # this fills by colum, i.e. 1 lonlat id by one. 
# # prices_ts <- ts(prices_ts, frequency = length_panel)
# 
# prices_ts <- ts(prices_extd[,c("cpo_price_imp1")], frequency = length_panel)
# 
# 
# 
# # fp <- prices[,c("lonlat", "year")]
# 
# # separate each time series (of each plantation) by an "empty" (NA) time series, to isolate them from each others. 
# # length_panel <- length(unique(prices$year))
# # fp <- mutate(fp, year = year + length_panel)
# # fp$cpo_price_imp1 <- NA
# # prices_extd <- rbind(prices, fp)
# # prices_extd <- dplyr::arrange(prices_extd, lonlat, year)
# # prices_ts <- prices_extd[,c("cpo_price_imp1")]
# # prices_ts <- ts(prices_ts, deltat = 1/length_panel*2)#frequency = length_panel
# 
# ur.kpss(prices_ts) %>% summary() # --> data need to be differenced
# ur.kpss(diff(prices_ts, differences = 1)) %>% summary() # --> first differencing suffices
# # ur.kpss(diff(prices_ts, differences = 2)) %>% summary() # --> first differencing suffices
# 
# arima_res <- stats::arima(prices_ts, order = c(length_panel - 1,1,0), method = "ML") 
# arima_res
# stargazer(arima_res, font.size = "footnotesize")
# 
# # ur.df(prices_ts)
# # this does not handle properly the specific structure of our time series (with NA sequences separating grid cells' respective time series)
# # pacf_14 <- pacf(diff(prices_ts, differences = 1), na.action = na.exclude, lag.max = length_panel - 1)
# 
# rm(RHS_50)




# ADDITIONAL REGRESSIONS --------------------------------------------------------



## Table A.2 ANNUAL EFFECTS ---------------------------------------------------------------------------------

# infrastructure to store results
res_data_list_full_annual <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      res_data_list_full_annual[[elm]] <- make_base_reg(island = ISL,
                                                        outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                        illegal = ILL,
                                                        annual = TRUE,
                                                        commo = "cpo",
                                                        n_iter_glm = 400,
                                                        offset = FALSE)
      names(res_data_list_full_annual)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    res_data_list_full_annual[[elm]] <- make_base_reg(island = ISL,
                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                      illegal = ILL,
                                                      annual = TRUE,
                                                      commo = "cpo",
                                                      n_iter_glm = 400,
                                                      offset = FALSE)
    names(res_data_list_full_annual)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full_annual)){
  # get estimation results and data; this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_full_annual[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_cumulative_APE(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

# Repeat without cumulative, i.e. annual APEs
rm(ape_mat_annual, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat_annual = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full_annual)){
  # get estimation results and data; this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_full_annual[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat_annual[[REGELM]] <- make_cumulative_APE(reg_elm = REGELM, cumulative = FALSE) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
ape_mat_annual <- bind_cols(ape_mat_annual)  %>% as.matrix()
row.names(ape_mat_annual) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat_annual)/2)-1)), "Observations", "Clusters") 
ape_mat_annual

# for(SIZE in size_list){
#   for(ILL in ill_status){
#     if(SIZE == "sm" & ILL == "ill2"){ # this is necessary to handle some convergence issue, but even with 10000 iterations algo does not converge. 
#       higher_n_iter <- 2000 # so we set it at 2000 to spare time
#     }else{higher_n_iter <- 400}
#     res_data_list_full_annual[[elm]] <- make_base_reg(island = ISL,
#                                                       outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                                                       illegal = ILL,
#                                                       annual = TRUE,
#                                                       commo = "cpo",
#                                                       n_iter_glm = higher_n_iter,
#                                                       offset = FALSE)
#     names(res_data_list_full_annual)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#     elm <- elm + 1
#   }
# }

ape_mat <- rbind(ape_mat[1:2,], ape_mat_annual)
colnames(ape_mat) <- NULL

# lapply(res_data_list_full_annual, function(x){x[[1]]})

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Cumulative and annual price elasticities of deforestation across Indonesian oil palm sectors") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Cumulative elasticity to:", 1, 2, 
            bold = TRUE)  %>%
  pack_rows("Annual CPO price signals", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Annual elasticity to:", 3, 10,
            bold = TRUE) %>%
  pack_rows("CPO price signal in t", 3, 4,
            italic = TRUE, bold = FALSE, indent = FALSE)  %>%
  pack_rows("CPO price signal in t-1", 5, 6,
            italic = TRUE, bold = FALSE, indent = FALSE)  %>%
  pack_rows("CPO price signal in t-2", 7, 8,
            italic = TRUE, bold = FALSE, indent = FALSE) %>%
  pack_rows("CPO price signal in t-3", 9, 10,
            italic = TRUE, bold = FALSE, indent = FALSE) %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>% 
  column_spec(column = 1,
              width = "12em",
              latex_valign = "m") 

## Table A.4 LAGGED OUTCOME ROBUSTNESS CHECKS -------------------------

res_data_list_lagov <- list()
SIZE <- "a"
ISL <- "both"
i <- 1

# Control for 5-year lagged deforestation 
res_data_list_lagov[[i]] <- make_base_reg(island = ISL,
                                          outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                          controls = c("n_reachable_uml",
                                                       "lucpfap_pixelcount_lag4_lag1"))


# Control for 5-year lagged deforestation in the neighbor grid cells
i <- i + 1
res_data_list_lagov[[i]] <- make_base_reg(island = ISL,
                                          outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                          controls = c("n_reachable_uml",
                                                       "ngb_ov_lag4_lag1"))

min_year_with_lags <- res_data_list_lagov[[i]][[2]] %>% pull(year) %>% min()
min_year_with_lags

# The main specification, on the same time period as with lags
i <- i + 1
res_data_list_lagov[[i]] <- make_base_reg(island = ISL,
                                          start_year = min_year_with_lags,
                                          outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                          controls = c("n_reachable_uml"))


# control for t-5 to t-8 average deforestation 
i <- i + 1
res_data_list_lagov[[i]] <- make_base_reg(island = ISL,
                                          outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                          controls = c("n_reachable_uml",
                                                       "lucpfap_pixelcount_lag4_4pya"))

# control for t-5 to t-8 average deforestation in the neighbor grid cells
i <- i + 1
res_data_list_lagov[[i]] <- make_base_reg(island = ISL,
                                          outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                          controls = c("n_reachable_uml",
                                                       "ngb_ov_lag4_4pya"))

min_year_with_lags <- res_data_list_lagov[[i]][[2]] %>% pull(year) %>% min()
min_year_with_lags

## The main specification, on the same time period as with lags
i <- i + 1
res_data_list_lagov[[i]] <- make_base_reg(island = ISL,
                                          start_year = min_year_with_lags,
                                          outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                          controls = c("n_reachable_uml"))
i

rm(i)

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_lagov)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_lagov[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticities of deforestation conditional on past deforestation") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Sample period: 2006-2014" = 3,
                     "Sample period: 2009-2014" = 3),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "In plantation site i" = 1,
                     "In i's 8 neighbors" = 1,
                     " " = 1,
                     "In plantation site i" = 1,
                     "In i's 8 neighbors" = 1,
                     " " = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c("Conditional on deforestation" = 1,
                     "the year preceding price signal" = 2,
                     "Not" = 1,
                     "the four years preceding price signal" = 2,
                     "Not" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c("Deforestation for:" = 1,
                     "All plantations" = 6),
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "10em",
              latex_valign = "m") 

rm(ape_mat)


## FALSIFICATION TEST (MILL-LEVEL) ####

ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))

# Restrict to analysis sample 
as <- ibs[ibs$analysis_sample==TRUE,] 

as <- mutate(as, in_kton_ffb_imp1 = in_ton_ffb_imp1/1000)

# make lagged var at mill level 
for(LAG in 1:3){
  as <- 
    DataCombine::slide(as,
                       Var = "in_kton_ffb_imp1",
                       TimeVar = "year",
                       GroupVar = "firm_id",
                       NewVar = paste0("in_kton_ffb_imp1_lag",LAG), 
                       slideBy = -LAG,
                       keepInvalid = TRUE)
  
  as <- 
    DataCombine::slide(as,
                       Var = "cpo_price_imp1",
                       TimeVar = "year",
                       GroupVar = "firm_id",
                       NewVar = paste0("cpo_price_imp1_lag",LAG), 
                       slideBy = -LAG,
                       keepInvalid = TRUE)
}

as <- 
  as %>% 
  mutate(
    # Make PYA for both dep and indep vars. 
    cpo_price_imp1_4pya   = rowMeans(across(.cols = starts_with("cpo_price_imp1")), na.rm = FALSE), 
    in_kton_ffb_imp1_4pya = rowMeans(across(.cols = starts_with("in_kton_ffb_imp1")), na.rm = FALSE), 
    # This would be the closest equivalent to the price signal used in used in analysis. 
    ln_cpo_price_imp1 = log(cpo_price_imp1), 
    ln_cpo_price_imp1_4pya = log(cpo_price_imp1_4pya)
  )

# as %>% dplyr::select(firm_id, year, starts_with("in_kton_ffb_imp1")) %>% View()

### Table A.5 CPO price in level ---------------------- 
# Regressions Price CPO ~ FFB vol 
fals_1 <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_1_wth <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 | district^year + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_2 <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_2_wth <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_3 <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_3_wth <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_4 <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 + in_kton_ffb_imp1_lag3 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_4_wth <-
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 + in_kton_ffb_imp1_lag3 | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_5 <- 
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1_4pya | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_5_wth <- 
  fixest::feols(fml = as.formula("cpo_price_imp1 ~ in_kton_ffb_imp1_4pya | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

# collect 
res_col <- list(fals_1,
                fals_1_wth,
                fals_2,
                fals_2_wth,
                fals_3,
                fals_3_wth,
                fals_4,
                fals_4_wth,
                fals_5,
                fals_5_wth)

df <- 
  etable(res_col,
    coefstat = "confint",
    tex = F, 
    se.below = TRUE,
    digits = "s3",
    signif.code=NA,
    depvar = FALSE,
    dict = c(cpo_price_imp1 = "CPO price",
             in_kton_ffb_imp1 = "FFB input",
             in_kton_ffb_imp1_lag1 = "FFB input t-1",
             in_kton_ffb_imp1_lag2 = "FFB input t-2",
             in_kton_ffb_imp1_lag3 = "FFB input t-3",
             in_kton_ffb_imp1_4pya = "FFB input 4 past year average", 
             firm_id = "mill")
  )

names(df)[1] <- "colnm"
df

df <- df %>% filter(!colnm %in% c("Fixed-Effects:", "_____________________", "_____________________________", 
                                  "S.E.: Clustered", "VCOV: Clustered", "R2", "Within R2"))
row.names(df) <- df$colnm
df <- df %>% dplyr::select(-colnm) %>% as.matrix()
# View(df)

df <- rbind(df, sapply(res_col, FUN = function(x) x$fixef_sizes[1]))
row.names(df) <- c(rep(c("Estimate","95% CI"), (nrow(df)-4)/2), "District-year", "Mill", "Observations", "Clusters") 
# df
colnames(df) <- NULL

options(knitr.table.format = "latex")
kable(df, booktabs = T, align = "c",
      caption = "Mill-level regressions of CPO output prices on FFB input volumes") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  
  add_header_above(c(" " = 1,
                     "CPO output price" = 10),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1, 
                     "(7)" = 1,
                     "(8)" = 1, 
                     "(9)" = 1, 
                     "(10)" = 1),
                   align = "c") %>%
  pack_rows("Regression coefficient of:", 1, 10, 
            bold = TRUE)  %>% 
  pack_rows(" FFB input kton in year t", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton in year t-1", 3, 4, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton in year t-2", 5, 6, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton in year t-3", 7, 8, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton, 4 past year average", 9, 10, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  
  pack_rows("Fixed-effects:", 11, 12,
            bold = TRUE)  %>% 
  pack_rows(" ", nrow(df)-1, nrow(df), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              underline = FALSE, # useless - line under "Deforestation for:" will need to be removed manually in Overleaf 
              width = "15em",
              latex_valign = "m") 



### Table A.6 CPO price in log ---------------------- 
# Regressions log(Price CPO) ~ FFB vol 
fals_1 <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_1_wth <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 | district^year + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_2 <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_2_wth <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_3 <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_3_wth <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_4 <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 + in_kton_ffb_imp1_lag3 | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_4_wth <-
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1 + in_kton_ffb_imp1_lag1 + in_kton_ffb_imp1_lag2 + in_kton_ffb_imp1_lag3 | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

fals_5 <- 
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1_4pya | district^year"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))
fals_5_wth <- 
  fixest::feols(fml = as.formula("ln_cpo_price_imp1 ~ in_kton_ffb_imp1_4pya | district^year  + firm_id"),
                data = as) # %>% filter(!is.na(in_kton_ffb_imp1_4pya)))

# collect 
res_col <- list(fals_1,
                fals_1_wth,
                fals_2,
                fals_2_wth,
                fals_3,
                fals_3_wth,
                fals_4,
                fals_4_wth,
                fals_5,
                fals_5_wth)

df <- 
  etable(res_col,
         coefstat = "confint",
         tex = F, 
         se.below = TRUE,
         digits = "s3",
         signif.code=NA,
         depvar = FALSE,
         dict = c(ln_cpo_price_imp1 = "CPO price",
                  in_kton_ffb_imp1 = "FFB input",
                  in_kton_ffb_imp1_lag1 = "FFB input t-1",
                  in_kton_ffb_imp1_lag2 = "FFB input t-2",
                  in_kton_ffb_imp1_lag3 = "FFB input t-3",
                  in_kton_ffb_imp1_4pya = "FFB input 4 past year average", 
                  firm_id = "mill")
  )

names(df)[1] <- "colnm"
df

df <- df %>% filter(!colnm %in% c("Fixed-Effects:", "_____________________", "_____________________________", 
                                  "S.E.: Clustered", "VCOV: Clustered", "R2", "Within R2"))
row.names(df) <- df$colnm
df <- df %>% dplyr::select(-colnm) %>% as.matrix()
View(df)

df <- rbind(df, sapply(res_col, FUN = function(x) x$fixef_sizes[1]))
row.names(df) <- c(rep(c("Estimate","95% CI"), (nrow(df)-4)/2), "District-year", "Mill", "Observations", "Clusters") 
df
colnames(df) <- NULL

options(knitr.table.format = "latex")
kable(df, booktabs = T, align = "c",
      caption = "Mill-level regressions of CPO output prices on FFB input volumes") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  
  add_header_above(c(" " = 1,
                     "Log of CPO output price" = 10),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1, 
                     "(7)" = 1,
                     "(8)" = 1, 
                     "(9)" = 1, 
                     "(10)" = 1),
                   align = "c") %>%
  pack_rows("Regression coefficient of:", 1, 10, 
            bold = TRUE)  %>% 
  pack_rows(" FFB input kton in year t", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton in year t-1", 3, 4, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton in year t-2", 5, 6, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton in year t-3", 7, 8, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" FFB input kton, 4 past year average", 9, 10, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  
  pack_rows("Fixed-effects:", 11, 12,
            bold = TRUE)  %>% 
  pack_rows(" ", nrow(df)-1, nrow(df), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              underline = FALSE, # useless - line under "Deforestation for:" will need to be removed manually in Overleaf 
              width = "15em",
              latex_valign = "m") 

## Table A.6 INTERACTION BY TYPE OF PLANTATION & INITIAL FOREST COVER -------------------------------------------------------------

# To do pooled regression with interaction, we need the residuals' variance to be the same across groups made by the interaction term. 
# Check whether this is the case with the level of clustering used in the main analysis (reachable): 

# First, get the pooled data and the residuals from the pooled regression
# For the three levels of heterogeneity, 
# Not for immediate conversion only, because in this case the equal variance assumption is not satisfied. 
res_data_list_byplantation <- list()
elm <- 1
ISL <- "both"
for(illegal_var in c("illegal2", "ill_or_concession", "ill_or_leg")){
  # if(Y == "lucpfip_rapid_pixelcount"){ends = 2010}else{ends = 2014}
  res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = "all",
                                                     outcome_variable = "lucpfip_pixelcount",
                                                     interaction_terms = c(illegal_var),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                     controls          = c(illegal_var, "n_reachable_uml"), #
                                                     n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                     offset = FALSE)
  names(res_data_list_byplantation)[elm] <- paste0(ISL,"_",illegal_var,"_i_byILL")
  elm <- elm + 1
  
  res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = "all",
                                                     outcome_variable = "lucpfip_pixelcount",
                                                     interaction_terms = c(illegal_var, "n_reachable_uml", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                     controls =          c(illegal_var, "n_reachable_uml", "pct_pfc2000_total"),
                                                     n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                     offset = FALSE)
  names(res_data_list_byplantation)[elm] <- paste0(ISL,"_",illegal_var,"_i_byILLandIFCandNREACH")
  elm <- elm + 1
  
  res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = "all",
                                                     outcome_variable = "lucpfip_pixelcount",
                                                     interaction_terms = c(illegal_var, "n_reachable_uml", "pct_pfc2000_total", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                     controls =          c(illegal_var, "n_reachable_uml", "pct_pfc2000_total", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp"),
                                                     n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                     offset = FALSE)
  names(res_data_list_byplantation)[elm] <- paste0(ISL,"_",illegal_var,"_i_byILLandIFCandNREACHandOWN")
  elm <- elm + 1
  
}

table_names = c("Estimate", 
                "95% CI", 
                "Estimate Illegality",
                "95% CI Illegality",
                "Estimate # reachable mills",
                "95% CI # reachable mills",
                "Estimate Initial forest cover (%)",
                "95% CI Initial forest cover (%)",
                "Estimate Mill ownership dom. private (%)",
                "95% CI Mill ownership dom. private (%)",
                "Estimate Mill ownership foreign (%)",
                "95% CI Mill ownership foreign (%)",
                "p-value",
                "Observations", 
                "Clusters")

tab_df = matrix(ncol = length(res_data_list_byplantation), 
                nrow = length(table_names)) 

row.names(tab_df) <- table_names
REGELM <- elm -1 
for(REGELM in 1:length(res_data_list_byplantation)){
  reg_name <- res_data_list_byplantation[REGELM] %>% names() %>% print()
  reg_res = res_data_list_byplantation[[REGELM]][[1]]
  d_clean = res_data_list_byplantation[[REGELM]][[2]]
  
  coeff_table =  etable(reg_res,
                        coefstat = "confint",
                        tex = F, 
                        se.below = TRUE,
                        digits = "s3",
                        signif.code=NA,
                        depvar = FALSE)
  
  row.names(coeff_table) = coeff_table[,1]
  # APE elasticity 
  ape_mat1 <- make_APEs(rounding = 3)
  # Works but useless.. 
  # if(nrow(ape_mat1)==6){row.names(ape_mat1)  = table_names[c(1:4, 10,11)]}
  # if(nrow(ape_mat1)==8){row.names(ape_mat1)  = table_names[c(1:6, 10,11)]}
  # if(nrow(ape_mat1)==10){row.names(ape_mat1) = table_names[c(1:8, 10,11)]}
  
  # Main price elasticity estimate and sample sizes
  tab_df[1:2,REGELM] <- ape_mat1[1:2]
  tab_df[(nrow(tab_df)-1):nrow(tab_df),REGELM] <- ape_mat1[(nrow(ape_mat1)-1):nrow(ape_mat1)] 
  
  # Coefficients 
  illegal_vars <- grep("ill", row.names(reg_res$coeftable), value = TRUE)
  illegal_var <- illegal_vars[1]
  illint_var <- illegal_vars[2]
  tab_df["Estimate Illegality", REGELM] <- reg_res$coeftable[illint_var,"Estimate"] %>% round(3)
  tab_df["95% CI Illegality",   REGELM] <- coeff_table[match(illint_var, row.names(coeff_table))+1,"reg_res"] %>% str_squish()

  tab_df["Estimate Illegality", REGELM] <- reg_res$coeftable[illint_var,"Estimate"] %>% round(3)
  tab_df["95% CI Illegality",   REGELM] <- coeff_table[match(illint_var, row.names(coeff_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate # reachable mills", REGELM] <- reg_res$coeftable["n_reachable_umlXln_wa_cpo_price_imp1_4ya_lag1" ,"Estimate"] %>% round(3)
  tab_df["95% CI # reachable mills",   REGELM] <- coeff_table[match("n_reachable_umlXln_wa_cpo_price_imp1_4ya_lag1", row.names(coeff_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Initial forest cover (%)", REGELM] <- reg_res$coeftable["pct_pfc2000_totalXln_wa_cpo_price_imp1_4ya_lag1" ,"Estimate"] %>% round(3)
  tab_df["95% CI Initial forest cover (%)",   REGELM] <- coeff_table[match("pct_pfc2000_totalXln_wa_cpo_price_imp1_4ya_lag1", row.names(coeff_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Mill ownership dom. private (%)", REGELM] <- reg_res$coeftable["wa_pct_own_nat_priv_impXln_wa_cpo_price_imp1_4ya_lag1" ,"Estimate"] %>% round(3)
  tab_df["95% CI Mill ownership dom. private (%)",   REGELM] <- coeff_table[match("wa_pct_own_nat_priv_impXln_wa_cpo_price_imp1_4ya_lag1", row.names(coeff_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Mill ownership foreign (%)", REGELM] <- reg_res$coeftable["wa_pct_own_for_impXln_wa_cpo_price_imp1_4ya_lag1" ,"Estimate"] %>% round(3)
  tab_df["95% CI Mill ownership foreign (%)",   REGELM] <- coeff_table[match("wa_pct_own_for_impXln_wa_cpo_price_imp1_4ya_lag1", row.names(coeff_table))+1,"reg_res"] %>% str_squish()
  
  # Equal variance tests
  d_clean$resid = reg_res$residuals
  # unclustered test
  res_variance <- d_clean %>%
    group_by(lonlat, year, !!as.symbol(illegal_var)) %>% 
    summarise(var_resid = mean(resid^2), .groups = "drop")
  t.test(var_resid ~ as.symbol(illegal_var), data = res_variance) %>% print()

  # clustered test
  res_variance <- d_clean %>%
    group_by(reachable, as.symbol(illegal_var)) %>% # reachable, 
    summarise(var_resid = mean(resid^2), .groups = "drop")
  eqvars_clust_test <- t.test(var_resid ~ as.symbol(illegal_var), data = res_variance) %>% print()
  
  # save p-value of clustered test 
  tab_df["p-value", REGELM] =  eqvars_clust_test$p.value %>% round(3)
}






REGELM <- elm - 1
res_data <- res_data_list_byplantation[[REGELM]]
reg_res <- res_data[[1]]
d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
rm(res_data)
ape_mat1 <- make_APEs(rounding = 3) # and for the same reason, this cannot be wrapped in other functions
reg_res %>% summary(cluster = "reachable")

res_data_list_byplantation[[REGELM]][[1]]$coeftable
## PARTIAL EFFECTS
# rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_byplantation)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_byplantation[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat1 <- make_APEs(rounding = 3) # and for the same reason, this cannot be wrapped in other functions
  reg_res %>% summary(cluster = "reachable")
  # arrange depending on interaction set - must make 12 rows in sum 
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUS"){ 
    # This has 6 rows. Add 6 rows after 4th one. 
    ape_mat1 <- rbind(matrix(ape_mat1[1:4,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 6, data = ""), 
                      matrix(ape_mat1[5:6,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUSandIFC"){ 
    # This has 8 rows. Add 4 rows after 4th one (4 last rows are IFC and N)
    ape_mat1 <- rbind(matrix(ape_mat1[1:4,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[5:8,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUSandILL"){ 
    # This has 8 rows. Add 4 rows after 6th one (2 last rows are N)
    ape_mat1 <- rbind(matrix(ape_mat1[1:6,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[7:8,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUSandILLandIFC"){ 
    # This has 10 rows. Add 2 rows after 6th one (4 last rows are IFC and N)
    ape_mat1 <- rbind(matrix(ape_mat1[1:6,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 2, data = ""), 
                      matrix(ape_mat1[7:10,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byILLINDUS"){ 
    # This has 6 rows. Add 4 rows after 2nd one, and 2 rows after the 4th one.  
    ape_mat1 <- rbind(matrix(ape_mat1[1:2,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[3:4,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 2, data = ""),
                      matrix(ape_mat1[5:6,])
                      )
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byILLINDUSandIFC"){ 
    # This has 8 rows. Add 4 rows after 2nd one
    ape_mat1 <- rbind(matrix(ape_mat1[1:2,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[3:8,])
    )
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byIFC"){ 
    # This has 6 rows. Add 6 rows after 2nd one
    ape_mat1 <- rbind(matrix(ape_mat1[1:2,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 6, data = ""), 
                      matrix(ape_mat1[3:6,])
    )
  }
  ape_mat[[REGELM]] <- ape_mat1
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)

ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL


options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticity heterogeneity by plantation types") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Deforestation for:" = 1,
                     "All plantations" = 6),
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  
  pack_rows("Moderation effect of:", 3, 10, 
            bold = TRUE, italic = FALSE)  %>%
  pack_rows("Share of industrial", 3, 4, indent = FALSE,
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Illegal", 5, 6, indent = FALSE,
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Share of industrial illegal", 7, 8, indent = FALSE,
            bold = FALSE, italic = TRUE)  %>%
  pack_rows("Share of initial forest cover", 9, 10, indent = FALSE,
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "12em",
              latex_valign = "m") 


rm(res_data_list_byplantation)


# 1. by INDUS 
res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
                                                   outcome_variable = paste0("lucpfap_pixelcount"),
                                                   interaction_terms =           c("share_indus"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                   controls = c("n_reachable_uml", "share_indus"),
                                                   n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                   offset = FALSE)
names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUS")
elm <- elm + 1

# 2. by INDUS and IFC
res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
                                                   outcome_variable = paste0("lucpfap_pixelcount"),
                                                   interaction_terms =           c("share_indus", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                   controls = c("n_reachable_uml", "share_indus", "pct_pfc2000_total"),
                                                   n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                   offset = FALSE)
names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUSandIFC")
elm <- elm + 1


# # 3. by INDUSxILL
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_pixelcount"),
#                                                    interaction_terms =           c("share_illindus"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUS_inINDUS")
# elm <- elm + 1
# # 4. by INDUSxILL and IFC, controlling for indus
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_pixelcount"),
#                                                    interaction_terms =           c("share_illindus", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus", "pct_pfc2000_total", "share_indus"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUSandIFC_ctrlINDUS")
# elm <- elm + 1

# # 3. by INDUS and ILL
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_indus", "illegal2"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_indus", "illegal2"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUSandILL")
# elm <- elm + 1
# 
# # 4. by INDUS and ILL and IFC
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_indus", "illegal2", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_indus", "illegal2", "pct_pfc2000_total"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUSandILLandIFC")
# elm <- elm + 1
# 
# # 5. by INDUSxILL
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_illindus"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUS")
# elm <- elm + 1
# 
# # 6. by INDUSxILL and IFC
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_illindus", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus", "pct_pfc2000_total"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUSandIFC")
# elm <- elm + 1


# # OR 
# ### For each plantation type 
# res_data_list_interact_ifc <- list()
# elm <- 1
# 
# size_list <- list("i","sm", "unr", "a")
# 
# # legality definition
# ill_def <- 2
# ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")
# 
# for(SIZE in size_list){
#   if(SIZE == "i"){
#     # Industrial by illegal status
#     for(ILL in ill_status){
#       res_data_list_interact_ifc[[elm]] <- make_base_reg(island = ISL,
#                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                      interaction_terms = c("pfc2000_total_pixelcount"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                      controls = c("n_reachable_uml", "pfc2000_total_pixelcount"),
#                                                      illegal = ILL,
#                                                      n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                      offset = FALSE)
#       names(res_data_list_interact_ifc)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#       elm <- elm + 1
#     }
#   } else {
#     # If size == sm or a, we do want to include all legal statuses.
#     # If size == unr, the selection of illegal indus is handled in make_base_reg 
#     ILL <- "all"
#     res_data_list_interact_ifc[[elm]] <- make_base_reg(island = ISL,
#                                                    outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                    interaction_terms = c("pfc2000_total_pixelcount"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "pfc2000_total_pixelcount"),
#                                                    illegal = ILL,
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
#     names(res_data_list_interact_ifc)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#     elm <- elm + 1
#   }}
# 
# ## PARTIAL EFFECTS
# rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
# ape_mat = list()
# reg_res_list = list()
# for(REGELM in 1:length(res_data_list_interact_ifc)){
#   # get estimation results and data
#   # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
#   res_data <- res_data_list_interact_ifc[[REGELM]]
#   reg_res <- res_data[[1]]
#   d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
#   rm(res_data)
#   ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM, rounding = 7) # and for the same reason, this cannot be wrapped in other functions
#   reg_res_list[[REGELM]] <- reg_res
#   rm(d_clean, reg_res)
# }
# # ape_mat <- lapply(res_data_list_interact_ifc, FUN = make_APEs) # this was the way to do it until it stopped working, on R 4.5 
# lapply(reg_res_list, FUN = function(x) x$convStatus)
# 
# ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
# row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
# ape_mat
# colnames(ape_mat) <- NULL


### Coeff difference tests --------------------
# This below is not used anymore. 
# compare_APEs_across_groups <- function(group1, group2, m0 = 0, alternative = "two.sided") { 
#   # important that it's selecting the "1" (resp. 2) row, and not the "Estimate" (resp. "SE") row in cases where there are several 
#   # rows with the same names, i.e. if APEs of interaction effects have been computed as well. 
#   ape1 <- ape_mat[1,group1] 
#   ape2 <- ape_mat[1,group2]
#   # # n1 <- ape_mat["Observations","Sumatra industrial"]
#   # # n2 <- ape_mat["Observations","Sumatra smallholders"]
#   sigma1 <- ape_mat[2,group1]^2
#   sigma2 <- ape_mat[2,group2]^2
#   
#   statistic <- (ape1 - ape2 - m0) / sqrt(sigma1 + sigma2)
#   
#   pval <- if (alternative == "two.sided") { 
#     2 * pnorm(abs(statistic), lower.tail = FALSE) 
#   } else if (alternative == "less") { 
#     pnorm(statistic, lower.tail = TRUE) 
#   } else { 
#     pnorm(statistic, lower.tail = FALSE) 
#   } 
#   # LCL <- (M1 - M2 - S * qnorm(1 - alpha / 2)) UCL <- (M1 - M2 + S * qnorm(1 - alpha / 2)) value <- list(mean1 = M1, mean2 = M2, m0 = m0, sigma1 = sigma1, sigma2 = sigma2, S = S, statistic = statistic, p.value = p, LCL = LCL, UCL = UCL, alternative = alternative) 
#   # print(sprintf("P-value = %g",p)) # print(sprintf("Lower %.2f%% Confidence Limit = %g", 
#   # alpha, LCL)) # print(sprintf("Upper %.2f%% Confidence Limit = %g", # alpha, UCL)) return(value) } test <- t.test_knownvar(dat1$sample1, dat1$sample2, V1 = 1, V2 = 1 )
#   return(pval)
# }

compare_coeffs_across_groups <- function(group1, group2, m0 = 0, alternative = "two.sided") { 
  # important that it's selecting the "1" (resp. 2) row, and not the "Estimate" (resp. "SE") row in cases where there are several 
  # rows with the same names, i.e. if APEs of interaction effects have been computed as well. 
  reg_sum1 <- res_data_list_full[[group1]][[1]] %>% summary()
  reg_sum2 <- res_data_list_full[[group2]][[1]] %>% summary()
  
  coeff1 <- reg_sum1$coefficients[1]
  coeff2 <- reg_sum2$coefficients[1]

  sigma1 <- reg_sum1$se[1]^2
  sigma2 <- reg_sum2$se[1]^2
  
  statistic <- (coeff1 - coeff2 - m0) / sqrt(sigma1 + sigma2)
  
  pval <- if (alternative == "two.sided") { 
    2 * pnorm(abs(statistic), lower.tail = FALSE) 
  } else if (alternative == "less") { 
    pnorm(statistic, lower.tail = TRUE) 
  } else { 
    pnorm(statistic, lower.tail = FALSE) 
  } 
  # LCL <- (M1 - M2 - S * qnorm(1 - alpha / 2)) UCL <- (M1 - M2 + S * qnorm(1 - alpha / 2)) value <- list(mean1 = M1, mean2 = M2, m0 = m0, sigma1 = sigma1, sigma2 = sigma2, S = S, statistic = statistic, p.value = p, LCL = LCL, UCL = UCL, alternative = alternative) 
  # print(sprintf("P-value = %g",p)) # print(sprintf("Lower %.2f%% Confidence Limit = %g", 
  # alpha, LCL)) # print(sprintf("Upper %.2f%% Confidence Limit = %g", # alpha, UCL)) return(value) } test <- t.test_knownvar(dat1$sample1, dat1$sample2, V1 = 1, V2 = 1 )
  return(pval)
}




## infrastructure to store results
ape_mat_list <- list()
elm <- 1

# make the APE. Note the difference: Standard errors are returnd with make_APEs_1regr, not CI95. 
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full)){
  res_data <- res_data_list_full[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs_1regr(res_data = res_data_list_full[[REGELM]]) # and for the same reason, this cannot be wrapped in other functions
  rm(d_clean, reg_res, res_data)
}
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()

rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full)){
  res_data <- res_data_list_full[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs() # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
colnames(ape_mat) <- names(res_data_list_full)
# keep ape_mat like this for later comparisons between APEs

# fill the matrix of p values of equality tests BY ROW

# comp_ape_mat["industrial = smallholders",paste0("both_a_all")] 
i_vs_sm <- compare_APEs_across_groups(group1 = paste0("both_i_all"), 
                                      group2 = paste0("both_sm_all")) 

notill_vs_ill <- compare_APEs_across_groups(group1 = paste0("both_i_no_ill2"), 
                                            group2 = paste0("both_i_ill2")) 

leg_vs_ill <- compare_APEs_across_groups(group1 = paste0("both_i_in_concession"), 
                                         group2 = paste0("both_i_ill2")) 


illindus_vs_sm <- compare_APEs_across_groups(group1 = paste0("both_i_ill2"), 
                                             group2 = paste0("both_sm_all")) 

# Old code to produce former table of p-values. 
# groups <- c("both_i_no_ill2", "both_i_ill2",
#             "both_i_all",
#             "both_sm_all")
# 
# comparisons <- c("industrial = smallholders",
#                  "legal indus. = illegal indus.")#,"primary forest = broadly def. forest"
# 
# comp_ape_mat <- matrix(ncol = length(groups), nrow = length(comparisons), data = NA)
# colnames(comp_ape_mat) <- groups
# row.names(comp_ape_mat) <- comparisons
# 
# for(SIZE in size_list){
#   comp_ape_mat["legal = illegal",paste0("both_",SIZE,"_all")] <- compare_APEs_across_groups(group1 = paste0("both_",SIZE,"_no_ill2"), 
#                                                                                             group2 = paste0("both_",SIZE,"_ill2")) 
# }
# 
# comp_ape_mat <- comp_ape_mat %>% formatC(digits = 4, format = "f")
# comp_ape_mat[comp_ape_mat=="   NA"] <- "" 
# comp_ape_mat
# colnames(comp_ape_mat) <- NULL
# 
# options(knitr.table.format = "latex")
# kable(comp_ape_mat, booktabs = T, align = "c",
#       caption = "p-values from equality tests of price elasticities") %>% #of 1 percentage change in medium-run price signal
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c("Ho" = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1,
#                      " " = 1,
#                      " " = 1),
#                    bold = F,
#                    align = "c") %>%  
#   add_header_above(c("Deforestation for:" = 1,
#                      "All plantations" = 3,
#                      "Industrial plantations" = 1,
#                      "Smallholder plantations" = 1),
#                    bold = F,
#                    align = "c") %>%
#   column_spec(column = 1,
#               width = "13em",
#               latex_valign = "m") %>% 
#   column_spec(column = c(2:(ncol(comp_ape_mat))),
#               width = "5em",
#               latex_valign = "m")



## Table A.7 PRICE VARIABILITY ----------------------------------------------------
# infrastructure to store results
res_data_list_variability <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      res_data_list_variability[[elm]] <- make_base_reg(island = ISL,
                                                        outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                        illegal = ILL,
                                                        price_variation = TRUE,
                                                        n_iter_glm = 500,
                                                        offset = FALSE)
      names(res_data_list_variability)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    res_data_list_variability[[elm]] <- make_base_reg(island = ISL,
                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                      illegal = ILL,
                                                      price_variation = TRUE,
                                                      n_iter_glm = 500,
                                                      offset = FALSE)
    names(res_data_list_variability)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_variability)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_variability[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Effects of price variability on deforestation across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Average Partial Effect of:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal variability", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "14em",
              latex_valign = "m") 


rm(ape_mat)

## Table A.8 ISLAND BREAKDOWN ----------------------------------------------------
# For temporal breakdown, check github before 30/06/2025 but even there it is old code that was not used at that time. 
ape_mat_list <- list()
isl_list <- list("Sumatra", "Kalimantan") 
for(ISL in isl_list){
  res_data_list_bd <- list()
  elm <- 1
  
  size_list <- list("i","sm", "unr", "a")
  
  # legality definition
  ill_def <- 2
  ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")
  
  for(SIZE in size_list){
    if(SIZE == "i"){
      # Industrial by illegal status
      for(ILL in ill_status){
        res_data_list_bd[[elm]] <- make_base_reg(island = ISL,
                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                 illegal = ILL,
                                                 n_iter_glm = 2000,
                                                 offset = FALSE)
        names(res_data_list_bd)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
        elm <- elm + 1
      }
    } else {
      # If size == sm or a, we do want to include all legal statuses.
      # If size == unr, the selection of illegal indus is handled in make_base_reg 
      ILL <- "all"
      res_data_list_bd[[elm]] <- make_base_reg(island = ISL,
                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                               illegal = ILL,
                                               n_iter_glm = 2000,
                                               offset = FALSE)
      names(res_data_list_bd)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }}
  
  ## PARTIAL EFFECTS
  # Collect partial effects for this island
  rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
  ape_mat = list()
  reg_res_list = list()
  for(REGELM in 1:length(res_data_list_bd)){
    # get estimation results and data
    # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
    res_data <- res_data_list_bd[[REGELM]]
    reg_res <- res_data[[1]]
    d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
    rm(res_data)
    ape_mat[[REGELM]] <- make_APEs() # and for the same reason, this cannot be wrapped in other functions
    reg_res_list[[REGELM]] <- reg_res
    rm(d_clean, reg_res)
  }
  lapply(reg_res_list, FUN = function(x) x$convStatus) %>% print()
  
  ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
  row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
  ape_mat
  colnames(ape_mat) <- NULL
  # Store the island specific results
  ape_mat_list[[ISL]] <- ape_mat
}

stacked_ape_mat <- rbind(ape_mat_list[[1]], ape_mat_list[[2]])

options(knitr.table.format = "latex")
kable(stacked_ape_mat, booktabs = T, align = "c",
      caption = "Price elasticities of deforestation across the oil palm sector, by island") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Sumatra elasticity to:", 1, 4, bold = TRUE)  %>% 
  pack_rows("CPO price signal", 1, 4, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", 3, 4, indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  
  pack_rows("Kalimantan elasticity to:", 5, 8, bold = TRUE)  %>%
  pack_rows("CPO price signal", 5, 8, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", 7, 8, indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "14em",
              latex_valign = "m") 




## Table A.9 SECONDARY FORESTS -----------------------

# infrastructure to store results
res_data_list_full_2ndry <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      res_data_list_full_2ndry[[elm]] <- make_base_reg(island = ISL,
                                                       outcome_variable = paste0("lucf",SIZE,"p_pixelcount"),
                                                       illegal = ILL,
                                                       offset = FALSE)
      names(res_data_list_full_2ndry)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    res_data_list_full_2ndry[[elm]] <- make_base_reg(island = ISL,
                                                     outcome_variable = paste0("lucf",SIZE,"p_pixelcount"),
                                                     illegal = ILL,
                                                     offset = FALSE)
    names(res_data_list_full_2ndry)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full_2ndry)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_full_2ndry[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
# ape_mat <- lapply(res_data_list_full_2ndry, FUN = make_APEs) # this was the way to do it until it stopped working, on R 4.5 
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

# for(SIZE in size_list){
#   for(ILL in ill_status){
#     res_data_list_full_2ndry[[elm]] <- make_base_reg(island = ISL,
#                                                      outcome_variable = paste0("lucf",SIZE,"p_pixelcount"), # note the absence of "p": it is not primary forest data. 
#                                                      illegal = ILL,
#                                                      offset = FALSE)
#     names(res_data_list_full_2ndry)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#     elm <- elm + 1
#   }
# }

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticities of secondary forest loss, across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Secondary forest loss for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, # indent = FALSE,
            bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "12em",
              latex_valign = "m")


rm(ape_mat)

## Table A.10 INTERACTION BY N REACHABLE -------------------------------------------------------------
res_data_list_interact <- list()
elm <- 1

size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      res_data_list_interact[[elm]] <- make_base_reg(island = ISL,
                                                     outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                     interaction_terms = c("n_reachable_uml"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                     illegal = ILL,
                                                     n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                     offset = FALSE)
      names(res_data_list_interact)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    res_data_list_interact[[elm]] <- make_base_reg(island = ISL,
                                                   outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                   interaction_terms = c("n_reachable_uml"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                   illegal = ILL,
                                                   n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                   offset = FALSE)
    names(res_data_list_interact)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_interact)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_interact[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM, rounding = 4) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
# ape_mat <- lapply(res_data_list_interact, FUN = make_APEs) # this was the way to do it until it stopped working, on R 4.5 
lapply(reg_res_list, FUN = function(x) x$convStatus)

ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

# for(SIZE in size_list){
#   for(ILL in ill_status){
#     res_data_list_interact[[elm]] <- make_base_reg(island = "both",
#                                                    outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                                                    interaction_terms = c("n_reachable_uml"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    illegal = ILL,
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
#     names(res_data_list_interact)[elm] <- paste0("both_",SIZE,"_",ILL)
#     elm <- elm + 1
#   }
# }

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticity heterogeneity across ownership and local market development") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  
  pack_rows("Moderation effect of:", 3, 4, 
            bold = TRUE, italic = FALSE)  %>%
  pack_rows("# reachable mills", 3, 4, indent = FALSE, 
            bold = FALSE, italic = TRUE)  %>%
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "10em",
              latex_valign = "m")


rm(res_data_list_interact)

## Table A.11 ILLEGAL MEASUREMENT ROBUSTNESS #### 
res_data_list_illrob <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","unr", "a")

# legality definition
ILL <- "alt"

for(SIZE in size_list){
  res_data_list_illrob[[elm]] <- make_base_reg(island = ISL,
                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                                               illegal = ILL,
                                               offset = FALSE)
  names(res_data_list_illrob)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
  elm <- elm + 1
}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_illrob)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_illrob[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price elasticities of illegal deforestation according to 2020 concession map") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 1,
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1),
                   align = "c") %>%pack_rows("Elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "10em",
              latex_valign = "m")

rm(ape_mat)


## Table A. 12 INTENSIVE MARGIN ROBUSTNESS -------------------------------------

# infrastructure to store results
res_data_list_intensive <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      if(ILL == "no_ill2"){# this is necessary to handle some convergence issue
        higher_iter_glm <- 2000 
      } else {
        higher_iter_glm <- 200
      } 
      res_data_list_intensive[[elm]] <- make_base_reg(island = ISL,
                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                      illegal = ILL,
                                                      margin = "intensive",
                                                      n_iter_glm = higher_iter_glm,
                                                      offset = FALSE)
      names(res_data_list_intensive)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    higher_iter_glm <- 200
    res_data_list_intensive[[elm]] <- make_base_reg(island = ISL,
                                                    outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                    illegal = ILL,
                                                    margin = "intensive",
                                                    n_iter_glm = higher_iter_glm,
                                                    offset = FALSE)
    names(res_data_list_intensive)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_intensive)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_intensive[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL

options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "CPO price elasticity of deforestation across Indonesian oil palm plantations, away from direct mill vicinity") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "10em",
              latex_valign = "m") 

rm(ape_mat)







## Table A.13 IN-LEVEL PRICES ROBUSTNESS ####
# infrastructure to store results
res_data_list_level <- list()
elm <- 1

isl_list <- list("both")#"Sumatra", "Kalimantan", 
ISL <- "both"

size_list <- list("i","sm", "unr", "a")

# legality definition
ill_def <- 2
ill_status <- c(paste0("no_ill",ill_def), paste0("ill",ill_def), "all")

for(SIZE in size_list){
  if(SIZE == "i"){
    # Industrial by illegal status
    for(ILL in ill_status){
      res_data_list_level[[elm]] <- make_base_reg(island = ISL,
                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                 illegal = ILL,
                                                 log_prices = FALSE,
                                                 n_iter_glm = 2000,
                                                 offset = FALSE)
      names(res_data_list_level)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
      elm <- elm + 1
    }
  } else {
    # If size == sm or a, we do want to include all legal statuses.
    # If size == unr, the selection of illegal indus is handled in make_base_reg 
    ILL <- "all"
    res_data_list_level[[elm]] <- make_base_reg(island = ISL,
                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                               illegal = ILL,
                                               log_prices = FALSE,
                                               n_iter_glm = 2000,
                                               offset = FALSE)
    names(res_data_list_level)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
    elm <- elm + 1
  }}

## PARTIAL EFFECTS
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_level)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_level[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM) # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)
# ape_mat <- lapply(res_data_list_level, FUN = make_APEs) # this was the way to do it until it stopped working, on R 4.5 

ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
ape_mat
colnames(ape_mat) <- NULL


options(knitr.table.format = "latex")
kable(ape_mat, booktabs = T, align = "c",
      caption = "Price semi-elasticities of deforestation across the Indonesian oil palm sector") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "Legal" = 1,
                     "Illegal" = 1,
                     "All" = 1, 
                     " " = 3),
                   bold = F,
                   align = "c") %>%
  add_header_above(c("Deforestation for:" = 1,
                     "Industrial plantations" = 3,
                     "Smallholder plantations" = 1, 
                     "Unregulated plantations" = 1, 
                     "All" = 1),
                   align = "c",
                   strikeout = F) %>%
  add_header_above(c(" " = 1,
                     "(1)" = 1,
                     "(2)" = 1, 
                     "(3)" = 1, 
                     "(4)" = 1,
                     "(5)" = 1, 
                     "(6)" = 1),
                   align = "c") %>%
  pack_rows("Semi-elasticity to:", 1, 2, bold = TRUE)  %>% # pack_rows(start_row =  nrow(ape_mat)-1, end_row = nrow(ape_mat),  latex_gap_space = "0.5em", hline_before = FALSE) %>% 
  pack_rows("CPO price signal", 1, 2, indent = FALSE, 
            bold = FALSE, italic = TRUE) %>% 
  pack_rows(" ", nrow(ape_mat)-1, nrow(ape_mat), indent = FALSE, latex_align = "l", latex_gap_space = "0em")  %>%
  column_spec(column = 1,
              width = "10em",
              latex_valign = "m") 


rm(ape_mat)




# Below are regressions of lagged deforestation on CPO price. 
# They do not produce a table for the paper but led to the robustness check above.
d <- rbind(d_30_suma, d_50_kali)

# make several lags of outcome var
for(ov in c("lucpfap_pixelcount_lag4", "ngb_ov_lag4")){
  for(LAG in 1:3){
    d <- 
      DataCombine::slide(d,
                         Var = ov,
                         TimeVar = "year",
                         GroupVar = "lonlat",
                         NewVar = paste0(ov, "_lag", LAG), 
                         slideBy = -LAG,
                         keepInvalid = TRUE)
  }}
for(LEAD in 1:3){
  d <- 
    DataCombine::slide(d,
                       Var = "lucpfap_pixelcount",
                       TimeVar = "year",
                       GroupVar = "lonlat",
                       NewVar = paste0("lucpfap_pixelcount_lead", LEAD), 
                       slideBy = LEAD,
                       keepInvalid = TRUE)
}
d <- 
  DataCombine::slide(d,
                     Var = "wa_cpo_price_imp1",
                     TimeVar = "year",
                     GroupVar = "lonlat",
                     NewVar = paste0("wa_cpo_price_imp1_lead8"), 
                     slideBy = 8,
                     keepInvalid = TRUE)
names(d)
d %>% filter(!is.na(lucpfap_pixelcount_lag4_lag3)) %>% pull(year) %>% summary()
d %>% filter(!is.na(wa_cpo_price_imp1_lead8)) %>% pull(year) %>% summary()

d <- 
  d %>% 
  mutate(lucpfap_pixelcount_lag4_4pya = rowMeans(across(.cols = c("lucpfap_pixelcount_lag4", 
                                                                  "lucpfap_pixelcount_lag4_lag1",
                                                                  "lucpfap_pixelcount_lag4_lag2",
                                                                  "lucpfap_pixelcount_lag4_lag3")), na.rm = FALSE),
         lucpfap_pixelcount_4fya = rowMeans(across(.cols = c("lucpfap_pixelcount", 
                                                             "lucpfap_pixelcount_lead1",
                                                             "lucpfap_pixelcount_lead2",
                                                             "lucpfap_pixelcount_lead3")), na.rm = FALSE),
         ngb_ov_lag4_4pya = rowMeans(across(.cols = c("ngb_ov_lag4", 
                                                      "ngb_ov_lag4_lag1",
                                                      "ngb_ov_lag4_lag2",
                                                      "ngb_ov_lag4_lag3")), na.rm = FALSE))

d %>% filter(!is.na(lucpfap_pixelcount_lag4_lag3)) %>% pull(year) %>% min()

# d$wa_cpo_price_imp1 <- d$lucpfap_pixelcount

# first end of period
fals_1 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1_lead8 ~ lucpfap_pixelcount + lucpfap_pixelcount_lead1 + lucpfap_pixelcount_lead2 + lucpfap_pixelcount_lead3 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))
fals_2 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1_lead8 ~ lucpfap_pixelcount_lead1 + lucpfap_pixelcount_lead2 + lucpfap_pixelcount_lead3 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))
fals_3 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1_lead8 ~ lucpfap_pixelcount_lead2 + lucpfap_pixelcount_lead3 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))
fals_4 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1_lead8 ~ lucpfap_pixelcount_lead3 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))
fals_5 <- 
  fixest::feols(fml = as.formula("wa_cpo_price_imp1_lead8 ~ lucpfap_pixelcount_4fya | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))


etable(
  fals_1,
  fals_2,
  fals_3,
  fals_4,
  fals_5,
  tex = F
)

# Second end of period
fals_1 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ lucpfap_pixelcount_lag4 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))

fals_2 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ lucpfap_pixelcount_lag4 + lucpfap_pixelcount_lag4_lag1 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))

fals_3 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ lucpfap_pixelcount_lag4 + lucpfap_pixelcount_lag4_lag1 + lucpfap_pixelcount_lag4_lag2 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))

fals_4 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ lucpfap_pixelcount_lag4 + lucpfap_pixelcount_lag4_lag1 + lucpfap_pixelcount_lag4_lag2 + lucpfap_pixelcount_lag4_lag3 | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))

fals_5 <- 
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ lucpfap_pixelcount_lag4_4pya | reachable + district^year"),
                data = d) # %>% filter(!is.na(lucpfap_pixelcount_lag4_4pya)))


etable(
  fals_1,
  fals_2,
  fals_3,
  fals_4,
  fals_5,
  tex = F
)

# with neighbors

fals_1 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ ngb_ov_lag4 | reachable + district^year"),
                data = d) # %>% filter(!is.na(ngb_ov_lag4_4pya)))

fals_2 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ ngb_ov_lag4 + ngb_ov_lag4_lag1 | reachable + district^year"),
                data = d) # %>% filter(!is.na(ngb_ov_lag4_4pya)))

fals_3 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ ngb_ov_lag4 + ngb_ov_lag4_lag1 + ngb_ov_lag4_lag2 | reachable + district^year"),
                data = d) # %>% filter(!is.na(ngb_ov_lag4_4pya)))

fals_4 <-
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ ngb_ov_lag4 + ngb_ov_lag4_lag1 + ngb_ov_lag4_lag2 + ngb_ov_lag4_lag3 | reachable + district^year"),
                data = d) # %>% filter(!is.na(ngb_ov_lag4_4pya)))

fals_5 <- 
  fixest::feols(fml = as.formula("wa_cpo_price_imp1 ~ ngb_ov_lag4_4pya | reachable + district^year"),
                data = d) # %>% filter(!is.na(ngb_ov_lag4_4pya)))


etable(
  fals_1,
  fals_2,
  fals_3,
  fals_4,
  fals_5,
  tex = F
)

rm(d)





# SPECIFICATION CHARTS -------------------------------------------------------------


# These charts compare coefficients from a single dependent variable across different models/specifications
# Therefore we produce one such chart for CPO TOTAL, FFB TOTAL (and possibly for LRs and SRs too)


# The function is written in base R and exactly copied from https://github.com/ArielOrtizBobea/spec_chart/blob/master/spec_chart_function.R
schart <- function(data, labels=NA, highlight=NA, n=1, index.est=1, index.se=2, index.ci=NA,
                   order="asis", ci=.95, ylim=NA, axes=T, heights=c(1,1), leftmargin=11, offset=c(0,0), ylab="Coefficient", lwd.border=1,
                   lwd.est=4, pch.est=21, lwd.symbol=2, ref=0, lwd.ref=1, lty.ref=2, col.ref="black", band.ref=NA, col.band.ref=NA,length=0,
                   col.est=c("grey60", "red3"), col.est2=c("grey80","lightcoral"), bg.est=c("white", "white"),
                   col.dot=c("grey60","grey95","grey95","red3"),
                   bg.dot=c("grey60","grey95","grey95","white"),
                   pch.dot=c(22,22,22,22), fonts=c(2,1), adj=c(1,1),cex=c(1,1)) {
  
  # Authors: Ariel Ortiz-Bobea (ao332@cornell.edu).
  # Version: March 10, 2020
  # If you like this function and use it, please send me a note. It might motivate
  # me to imporove it or write new ones to share.
  
  # Description of arguments
  
  # Data:
  # data: data.frame with data, ideally with columns 1-2 with coef and SE, then logical variables.
  # labels: list of labels by group. Can also be a character vector if no groups. Default is rownames of data.
  # index.est: numeric indicating position of the coefficient column.
  # index.se: numeric indicating position of the SE column.
  # index.ci: numeric vector indicating position of low-high bars for SE. Can take up to 2 CI, so vector can be up to length 4
  
  # Arrangement and basic setup:
  # highlight: numeric indicating position(s) of models (row) to highlight in original dataframe.
  # n: size of model grouping. n=1 removes groupings. A vector yields arbitrary groupings.
  # order: whether models should be sorted or not. Options: "asis", "increasing", "decreasing"
  # ci: numeric indicating level(s) of confidence. 2 values can be indicated.
  # ylim: if one wants to set an arbitrary range for Y-axis
  
  # Figure layout:
  # heights: Ratio of top/bottom panel. Default is c(1,1) for 1 50/50 split
  # leftmargin: amount of space on the left margin
  # offset: vector of numeric with offset for the group and specific labels
  # ylab: Label on the y-axis of top panel. Default is "Coefficient"
  # lwd.border: width of border and other lines
  
  # Line and symbol styles and colors:
  # lwd.est: numeric indicating the width of lines in the top panel
  # ref: numeric vector indicating horizontal reference lines(s) Default is 0.
  # lty.ref: Style of reference lines. Default is dash line (lty=2).
  # lwd.ref. Width of reference lines. Default is 1.
  # col.ref: vector of colors of reference lines. Default is black.
  # band.ref: vector of 2 numerics indicating upper abdn lower height for a band
  # col.band.ref: color of this band
  # col.est: vector of 2 colors indicating for "other" and "highlighted" models
  # col.est2: same for outer confidence interval if more than 1 confidence interval
  # col.dot: vector of 4 colors indicating colors for borders of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # bg.dot : vector of 4 colors indicating colors for background of symbol in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # pch.dot: style of symbols in bottom panel for "yes", "no", "NA", and "yes for highlighted model"
  # length: length of the upper notch on th vertical lines. default is 0.
  
  # Letter styles
  # fonts: numeric vector indicating font type for group (first) and other labels (second) (e.g. 1:normal, 2:bold, 3:italic)
  # adj: numeric vector indicating alignment adjustment for text label: 0 is left, .5 is center, 1 is right.
  # cex: numeric vector for size of fonts for top panel (first) and bottom panel (Second)
  
  # 1. Set up
  if (T) {
    # Arrange data
    d <- data
    rownames(d) <- 1:nrow(d)
    
    # Create ordering vector
    if (order=="asis")       o <- 1:length(d[,index.est])
    if (order=="increasing") o <- order(d[,index.est])
    if (order=="decreasing") o <- order(-d[,index.est])
    if (!is.numeric(d[,index.est])) {warning("index.est does not point to a numeric vector.") ; break}
    d <- d[o,]
    est <- d[,index.est] # Estimate
    if (length(index.ci)>1) {
      l1 <- d[,index.ci[1]]
      h1 <- d[,index.ci[2]]
      if (length(index.ci)>2) {
        l2 <- d[,index.ci[3]]
        h2 <- d[,index.ci[4]]
      }
    } else {
      if (!is.numeric(d[,index.se]))  {warning("index.se does not point to a numeric vector.") ; break}
      se  <- d[,index.se] # Std error
      ci <- sort(ci)
      a <- qnorm(1-(1-ci)/2)
      l1 <- est - a[1]*se
      h1 <- est + a[1]*se
      if (length(ci)>1) {
        l2 <- est - a[2]*se
        h2 <- est + a[2]*se
      }
    }
    
    # Table
    if (length(index.ci)>1) remove.index <- c(index.est,index.ci) else remove.index <- c(index.est,index.se)
    remove.index <- remove.index[!is.na(remove.index)]
    tab <- t(d[,-remove.index]) # get only the relevant info for bottom panel
    if (!is.list(labels) & !is.character(labels)) labels <- rownames(tab)
    
    # Double check we have enough labels
    if ( nrow(tab) != length(unlist(labels))) {
      print("Warning: number of labels don't match number of models.")
      labels <- rownames(tab)
    }
    
    # Plotting objects
    xs <- 1:nrow(d) # the Xs for bars and dots
    if (n[1]>1 & length(n)==1) xs <- xs + ceiling(seq_along(xs)/n) - 1 # group models by n
    if (length(n)>1) {
      if (sum(n) != nrow(d) ) {
        warning("Group sizes don't add up.")
      } else {
        idx <- unlist(lapply(1:length(n), function(i) rep(i,n[i])))
        xs <- xs + idx - 1
      }
    }
    h <- nrow(tab) + ifelse(is.list(labels),length(labels),0) # number of rows in table
    # Location of data and labels
    if (is.list(labels)) {
      index <- unlist(lapply(1:length(labels), function(i) rep(i, length(labels[[i]])) ))
      locs <- split(1:length(index),index)
      locs <- lapply(unique(index), function(i) {
        x <- locs[[i]]+i-1
        x <- c(x,max(x)+1)
      })
      yloc  <- unlist(lapply(locs, function(i) i[-1])) # rows where data points are located
      yloc2 <- sapply(locs, function(i) i[1]) # rows where group lables are located
    } else {
      yloc <- 1:length(labels)
    }
    
    # Range
    if (is.na(ylim[1]) | length(ylim)!=2) {
      if (length(index.ci)>2 | length(ci)>1) {
        ylim <- range(c(l2,h2,ref)) # range that includes reference lines
      } else {
        ylim <- range(c(l1,h1,ref))
      }
      ylim <- ylim + diff(ylim)/10*c(-1,1) # and a bit more
    }
    xlim <- range(xs) #+ c(1,-1)
  }
  
  # 2. Plot
  if (T) {
    #par(mfrow=c(2,1), mar=c(0,leftmargin,0,0), oma=oma, xpd=F, family=family)
    layout(t(t(2:1)), height=heights, widths=1)
    par(mar=c(0,leftmargin,0,0), xpd=F)
    
    # Bottom panel (plotted first)
    plot(est, xlab="", ylab="", axes=F, type="n", ylim=c(h,1), xlim=xlim)
    lapply(1:nrow(tab), function(i) {
      # Get colors and point type
      type <- ifelse(is.na(tab[i,]),3,ifelse(tab[i,]==TRUE,1, ifelse(tab[i,]==FALSE,2,NA)))
      type <- ifelse(names(type) %in% paste(highlight) & type==1,4,type) # replace colors for baseline model
      col <- col.dot[type]
      bg  <- bg.dot[type]
      pch <- as.numeric(pch.dot[type])
      sel <- is.na(pch)
      # Plot points
      points(xs, rep(yloc[i],length(xs)), col=col, bg=bg, pch=pch, lwd=lwd.symbol)
      points(xs[sel], rep(yloc[i],length(xs))[sel], col=col[sel], bg=bg[sel], pch=pch.dot[3]) # symbol for missing value
      
    })
    par(xpd=T)
    if (is.list(labels)) text(-offset[1], yloc2, labels=names(labels), adj=adj[1], font=fonts[1], cex=cex[2])
    # Does not accomodate subscripts
    text(-rev(offset)[1], yloc , labels=unlist(labels), adj=rev(adj)[1], font=fonts[2], cex=cex[2])
    # Accomodates subscripts at the end of each string
    if (F) {
      labels1 <- unlist(labels)
      lapply(1:length(labels1), function(i) {
        a  <- labels1[i]
        a1 <- strsplit(a,"\\[|\\]")[[1]][1]
        a2 <- rev(strsplit(a,"\\[|\\]")[[1]])[1]
        if (identical(a1,a2))  a2 <- NULL
        text(-rev(offset)[1], yloc[i], labels=bquote(.(a1)[.(a2)]), adj=adj[2], font=fonts[2], cex=cex[2])
      })
    }
    par(xpd=F)
    
    # Top panel (plotted second)
    colvec  <- ifelse(colnames(tab) %in% paste(highlight), col.est[2], col.est[1])
    bg.colvec  <- ifelse(colnames(tab) %in% paste(highlight), bg.est[2], bg.est[1])
    colvec2 <- ifelse(colnames(tab) %in% paste(highlight),col.est2[2], col.est2[1])
    plot(est, xlab="", ylab="", axes=F, type="n", ylim=ylim, xlim=xlim)
    # Band if present
    if (!is.na(band.ref[1])) {
      rect(min(xlim)-diff(xlim)/10, band.ref[1], max(xlim)+diff(xlim)/10, band.ref[2], col=col.band.ref, border=NA)
    }
    # Reference lines
    abline(h=ref, lty=lty.ref, lwd=lwd.ref, col=col.ref)
    # Vertical bars
    if (length(ci)>1 | length(index.ci)>2) arrows(x0=xs, y0=l2, x1=xs, y1=h2, length=length, code=3, lwd=rev(lwd.est)[1], col=colvec2, angle=90)
    arrows(x0=xs, y0=l1, x1=xs, y1=h1, length=length, code=3, lwd=lwd.est[1]     , col=colvec, angle=90)
    points(xs, est, pch=pch.est, lwd=lwd.symbol, col=colvec, bg=bg.colvec)
    # Axes
    if (axes) {
      axis(2, las=2, cex.axis=cex[1], lwd=lwd.border)
      axis(4, labels=NA, lwd=lwd.border)
    }
    mtext(ylab, side=2, line=3.5, cex=cex[1])
    box(lwd=lwd.border)
    
  }
  
} 


## We now produce the "data" argument for the function schart, i.e. we run one regression, and bind in a dataframe 
# the coefficients, the SE, and indicator variables for the correponding specifications. 
make_spec_chart_df <- function(island,
                               start_year = 2002, 
                               end_year = 2014,
                               outcome_variable = "lucpfip_pixelcount", # LHS. One of "lucfip_pixelcount", "lucfip_pixelcount_60th", "lucfip_pixelcount_90th", "lucpfip_pixelcount_intact", "lucpfip_pixelcount_degraded", "lucpfip_pixelcount"p
                               PS = 3000,
                               catchment = "CR",  
                               alt_cr = FALSE, # logical, if TRUE, Sumatra's catchment radius is 50000 meters, and Kalimantan's is 30000. 
                               nearest_mill = FALSE, # whether the ibs variables should be attributed to parcels as from the nearest mill or inverse distance weighted average.  
                               margin = "both",
                               restr_marg_def = TRUE, # should the restrictive definition of extensive margin be used, if margin is either "intensive" of "extensive"
                               commo = "cpo", # either "ffb", "cpo", or c("ffb", "cpo"), commodities the price signals of which should be included in the RHS
                               x_pya = 3, # either 2, 3, or 4. The number of past years to compute the average of rhs variables over. The total price signal is the average over these x_pya years and the current year. 
                               dynamics = FALSE, # Logical, should the total price signal(s) be split into current year and x_pya past year average. 
                               yoyg = FALSE, # logical, should the price variables be computed in year-on-year growth rate instead of level.
                               only_sr = FALSE,
                               log_prices = TRUE, # Logical, should the price variables be included as their logarithms instead of levels. No effect if yoyg is TRUE.    
                               short_run = "full", # either "full", or "dev". Used only if dynamics = TRUE and yoyg = FALSE. Should the short run (SR) measure of regressors be the price signal as such ("full") or be the deviation to past year average price signal ("dev"). 
                               imp = 1, # either 1 or 2. 1 selects the data cleaned with the stronger imputations. 
                               distribution = "quasipoisson", # either "poisson", "quasipoisson", or "negbin"
                               fe = "reachable + district_year", # fixed-effects, interactions should not be specified in {fixest} synthax with fe1^fe2
                               offset = FALSE, # Logical. Should the log of the remaining forest be added as an offset.  
                               lag_or_not = "_lag1", # either "_lag1", or  "", should the 
                               controls = c("n_reachable_uml"), # , "wa_prex_cpo_imp1"character vectors of names of control variables (don't specify lags in their names)
                               remaining_forest = FALSE, # Logical. If TRUE, the remaining forest is added as a control
                               interaction_terms = NULL, # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1"), # may be one or several of the controls specified above. 
                               interacted = "regressors",
                               interact_regressors = TRUE, # if there are two regressors (e.g. ffb and cpo), should their interaction be included in the model? 
                               pya_ov = FALSE, # logical, whether the pya (defined by x_pya) of the outcome_variable should be added in controls
                               illegal = "all", # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                               min_forest_2000 = 0, 
                               min_coverage = 0, # fraction, from 0 to 1. Minimum share of reachable IBS over all reachable (UML), for an obs.to be included in sample. 
                               weights = FALSE, # logical, should obs. be weighted by the share of our sample reachable mills in all (UML) reachable mills. 
                               # additional argument to this function: 
                               variable = "cpo",
                               cluster = "reachable" # "subdistrict"
){
  ### Get coefficient and SE for a given specification passed to make_base_reg.
  # call make_base_reg
  res_data_list <- make_base_reg(island = island,
                                       start_year = start_year, 
                                       end_year = end_year,
                                       outcome_variable = outcome_variable,
                                       PS = PS,
                                       catchment = catchment,
                                       alt_cr = alt_cr,  
                                       nearest_mill = nearest_mill, 
                                       margin = margin,
                                       restr_marg_def = restr_marg_def,
                                       commo = commo, 
                                       x_pya = x_pya, 
                                       dynamics = dynamics, 
                                       yoyg = yoyg, 
                                       only_sr = only_sr,
                                       log_prices = log_prices,
                                       short_run = short_run, 
                                       imp = imp, 
                                       distribution = distribution, 
                                       fe = fe, 
                                       offset = offset,
                                       lag_or_not = lag_or_not, 
                                       controls = controls,
                                       remaining_forest = remaining_forest, 
                                       interaction_terms = interaction_terms, # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp","n_reachable_uml", "wa_prex_cpo_imp1"), # may be one or several of the controls specified above. 
                                       interacted = interacted,
                                       interact_regressors = interact_regressors, # if there are two regressors (e.g. ffb and cpo), should their interaction be included in the model? 
                                       pya_ov = pya_ov, 
                                       illegal = illegal, # the default, "all" includes all data in the regression. "ill1" and "ill2" (resp. "no_ill1" and "no_ill2") include only illegal (resp legal) lucfp (two different definitions, see add_parcel_variables.R)
                                       min_forest_2000 = min_forest_2000, 
                                       min_coverage = min_coverage,
                                       weights = weights)
  
  # compute SE the way specified by cluster
  # get estimation results and data
  # For some reasons, this has to be done within the function, as before, unlike the main APE function which is called in the global env. directly... 
  # reg_res <- res_data_list[[1]]
  # d_clean <- res_data_list[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  # rm(d_clean, reg_res)
  ape_mat <- make_APEs_1regr(res_data = res_data_list, CLUSTER = cluster, stddev = FALSE) # and for the same reason, this cannot be wrapped in other functions
  rm(d_clean, reg_res, res_data_list)

  # ape_mat <- make_APEs_1regr(res_data_list, CLUSTER = cluster, stddev = FALSE) 
  # reg_summary <- summary(reg_res, se = "cluster", cluster = cluster)
  # get coeff and SE
  # coeff_and_se <- reg_summary$coeftable[match(variable, commo), 1:2]
  
  ### make indicator variables that will be used to label specifications. 
  # /!\ THE ORDER HERE MATTERS !
  ind_var <- data.frame(
    # outcome variable
    #"larger_forest_def" = FALSE,
    #"rapid" = FALSE,
    ## Data cleaning
    "imp1" = FALSE,
    "lag_or_not" = FALSE,
    ## Sampling
    "min_forest_2000" = FALSE,
    "min_coverage" = FALSE,
    ## Parcel size 
    "PS_1KM" = FALSE,
    "PS_3KM" = FALSE,
    "PS_5KM" = FALSE,
    ## Catchment model
    #"nearest_mill" = FALSE,
    #"intensive_only" = FALSE,
    "alt_cr" = FALSE,
    "CA" = FALSE,
    # Price signal past year average
    "only_sr" = FALSE,
    "ya_3" = FALSE,
    "ya_4" = FALSE,
    "ya_5" = FALSE,
    # distribution assumptions
    "poisson" = FALSE,
    "quasipoisson" = FALSE,
    "negbin" = FALSE,
    # controls
    # "two_commo" = FALSE,
    # "price_dynamics" = FALSE,
    "n_reachable_uml_control" = FALSE, 
    "control_own" = FALSE,
    "ov_lag1" = FALSE,
    # "ov_lag4" = FALSE,
    # "ngb_ov_lag4" = FALSE,
    "prex_cpo_control" = FALSE,
    "baseline_forest_trend" = FALSE,
    #"remaining_forest" = FALSE,
    # interactions
    # "own_interact" = FALSE,
    # "n_reachable_uml_interact" = FALSE,
    # "prex_cpo_interact" = FALSE,
    # "remaining_forest_interact"  = FALSE,
    # weights
    #"weights" = FALSE,
    # fixed effects
    "unit_fe" = FALSE,
    "tw_fe" = FALSE,
    "unit_provyear_fe" = FALSE,
    "unit_distryear_fe" = FALSE,
    "unit_subdistryear_fe" = FALSE,
    "plantation_distryear_fe" = FALSE,
    # "unit_villageyear_fe" = FALSE,
    # standard errors
    "unit_cluster" = FALSE, 
    "reachable_cluster" = FALSE,
    "village_cluster" = FALSE,
    "subdistrict_cluster" = FALSE,
    "district_cluster" = FALSE,
    "twoway_cluster" = FALSE
  )
  
  # change the OV indicator variable to TRUE 
  #if(grepl("lucf", outcome_variable)){ind_var[,"larger_forest_def"] <- TRUE}
  #if(grepl("rapid",outcome_variable)){ind_var[,"rapid"] <- TRUE}
  
  ## Data cleaning
  # set the indicator variable for the data imputation 
  if(imp == 1){ind_var[,"imp1"] <- TRUE}
  # set indicator variable for lagging or not
  if(lag_or_not == "_lag1"){ind_var[,"lag_or_not"] <- TRUE}
  
  ## sampling
  # minimums 
  if(min_forest_2000>0){ind_var[,"min_forest_2000"] <- TRUE}
  if(min_coverage>0){ind_var[,"min_coverage"] <- TRUE}
  
  if(PS==1000){ind_var[,"PS_1KM"] <- TRUE}
  if(PS==3000){ind_var[,"PS_3KM"] <- TRUE}
  if(PS==5000){ind_var[,"PS_5KM"] <- TRUE}
  
  ## Catchment model
  # nearest mill rather than inverse-distance weighted average
  if(nearest_mill){ind_var[,"nearest_mill"] <- TRUE}
  # margin
  #if(margin == "intensive"){ind_var[,"intensive_only"] <- TRUE}
  # Catchment 
  if(alt_cr){ind_var[,"alt_cr"] <- TRUE} 
  if(catchment == "CA"){ind_var[,"CA"] <- TRUE} 
  
  ## Past year average
  # only Short Run
  if(only_sr){ind_var[,"only_sr"] <- TRUE}
  # set # year average indicator variables (+1 because we look at the total effect)
  if(!only_sr){# otherwise, the default
    ind_var[,grepl(paste0("ya_", x_pya+1), colnames(ind_var))] <- TRUE
  }
  ## set the indicator variables for the distribution assumptions
  ind_var[,distribution] <- TRUE
  
  # set indicator variables for fixed effects
  if(fe == "reachable"){ind_var[,"unit_fe"] <- TRUE}
  if(fe == "reachable + year"){ind_var[,"tw_fe"] <- TRUE}
  if(fe == "reachable + province_year"){ind_var[,"unit_provyear_fe"] <- TRUE}
  if(fe == "reachable + district_year"){ind_var[,"unit_distryear_fe"] <- TRUE}
  if(fe == "reachable + subdistrict_year"){ind_var[,"unit_subdistryear_fe"] <- TRUE}
  if(fe == "reachable + village_year"){ind_var[,"unit_villageyear_fe"] <- TRUE}
  if(fe == "lonlat + district_year"){ind_var[,"plantation_distryear_fe"] <- TRUE}
  
  ## set indicator variables for controls
  # second commo 
  # if(length(commo) == 2){ind_var[,"two_commo"] <- TRUE}
  # if(dynamics){ind_var[,"price_dynamics"] <- TRUE}
  # ownership
  if(any(grepl("pct_own", controls))){ind_var[,"control_own"] <- TRUE}
  # n_reachable_uml
  if(any(grepl("n_reachable_uml", controls))){ind_var[,"n_reachable_uml_control"] <- TRUE}
  # 1-year lagged outcome variable
  if(any(grepl(paste0(outcome_variable, "_lag1"), controls))){ind_var[,"ov_lag1"] <- TRUE}
  # # 4-year lagged outcome variable
  # if(any(grepl(paste0(outcome_variable, "_lag4"), controls))){ind_var[,"ov_lag4"] <- TRUE}
  # # spatial lag deforestation
  # if(any(grepl("ngb_ov_lag4", controls))){ind_var[,"ngb_ov_lag4"] <- TRUE}
  # prex_cpo
  if(any(grepl("prex_cpo", controls))){ind_var[,"prex_cpo_control"] <- TRUE}
  # baseline forest trend
  if("baseline_forest_trend" %in% controls){ind_var[,"baseline_forest_trend"] <- TRUE}
  # remaining forest
  if("remain_pf_pixelcount" %in% controls){ind_var[,"remaining_forest"] <- TRUE}
  
  ## interactions
  # if(any(grepl("own_nat_priv",interaction_terms)) | any(grepl("own_for",interaction_terms))){ind_var[,"own_interact"] <- TRUE}
  # if(any(grepl("n_reachable_uml",interaction_terms))){ind_var[,"n_reachable_uml_interact"] <- TRUE}
  # if(any(grepl("prex_cpo",interaction_terms))){ind_var[,"prex_cpo_interact"] <- TRUE}
  # if(any(grepl("remain_pf",interaction_terms))){ind_var[,"remaining_forest_interact"] <- TRUE}
  
  
  # weights indicator variable
  #if(weights){ind_var[,"weights"] <- TRUE}
  #if(distribution == "negbin"){ind_var[,"weights"] <- FALSE} # turn it back to FALSE if negbin distribution
  
  # clustering
  if(length(cluster) == 1){
    if(cluster == "lonlat"){ind_var[,"unit_cluster"] <- TRUE}
    if(cluster == "reachable"){ind_var[,"reachable_cluster"] <- TRUE}
    if(cluster == "village"){ind_var[,"village_cluster"] <- TRUE}
    if(cluster == "district"){ind_var[,"district_cluster"] <- TRUE}
    if(cluster == "subdistrict"){ind_var[,"subdistrict_cluster"] <- TRUE}
  }else if(length(cluster) == 2){
    ind_var[,"twoway_cluster"] <- TRUE
  }
  
  ### Bind together coeff, SE, and specification labels
  spec_df <- cbind(ape_mat[1,],ape_mat[2,], ind_var)
  
  return(spec_df)
}


### COMPUTE THE DATASETS FOR EACH ISLAND, AND COMMODITY OF INTEREST

# There is no loop here, because it makes too much code to repeat. Just change values below and rerun everything. 
SIZE <- "sm" # "i" # "sm", "i" # 
ILL <- "all" # "ill2" # "ill2" # 
ISL <- "both"

# for(ISL in c("Sumatra", "Kalimantan", "both")){
#   for(SIZE in c("i","sm")){


reg_stats_indvar_list <- list()
i <- 1

### Add to the list the specifications that are to be compared in the chart. 
# These are particular departures from the preferred specification (which arguments are already set by default)
# Specifying only arguments for which it's changing + those defining the estimate we are interested in plotting the spec chart. 

## the main specification
reg_stats_indvar_list[["main"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"))
i <- i+1

## Sampling 
# Minimum forest coverage in 2000
reg_stats_indvar_list[["min_forest_2000"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                                 outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                 min_forest_2000 = 0.5)
i <- i+1

# Minimum IBS/UML ratio 
reg_stats_indvar_list[["min_coverage"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                              outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                              min_coverage = 0.5)
i <- i+1

# Alternative parcel size 
reg_stats_indvar_list[["PS_1km"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                        outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                        PS = 1000)
i <- i+1
reg_stats_indvar_list[["PS_5km"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                        outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                        PS = 5000)
i <- i+1

## Data cleaning
# Stronger imputations (imp2)
reg_stats_indvar_list[["imp"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                     outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                     imp = 2)
i <- i+1


# Not lagged controls
reg_stats_indvar_list[["lag_or_not"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                            outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                            lag_or_not = "")
i <- i+1

## For an alternative forest definition
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL,
#                                                  outcome_variable = paste0("lucf",SIZE,"p_pixelcount")) 
# i <- i+1

## Alternative past year average
reg_stats_indvar_list[["only_sr"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                         outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                         only_sr = TRUE)
i <- i+1

for(XPYA in c(2, 4)){
  reg_stats_indvar_list[[paste0(XPYA,"_pya")]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                                     outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                     x_pya = XPYA)
  i <- i+1
}

## Catchment model 

# Nearest mill
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  nearest_mill = TRUE)
# i <- i+1

# For intensive margin only
# reg_stats_indvar_list[["margin"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                         outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                         margin = "intensive") 
# i <- i+1

# For alternative catchment radius
reg_stats_indvar_list[["alt_cr"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                        outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                        alt_cr = TRUE) 
i <- i+1  

# For catchment area
reg_stats_indvar_list[["catchment"]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                           outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                           cluster = "reachable", # variable reachable is not computed in the CA stream. 
                                                           catchment = "CA") 
i <- i+1


i # 12

## RAPID LUCFP - this is not run anymore. 
# if(SIZE == "i"){
#   # loops over critical parameters (not repeating because the outcome is different)
#   for(IMP in c(1, 2)){
#     for(XPYA in c(2, 3, 4)){
#       for(LAG in c("_lag1", "")){
#         reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                          outcome_variable = paste0("lucpfip_rapid_pixelcount"),
#                                                          imp = IMP,
#                                                          x_pya = XPYA,
#                                                          lag_or_not = LAG)
#         i <- i+1
#       }
#     }
#   }    
#   i # 28
#   
#   # in CA
#   reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_rapid_pixelcount"),
#                                                    catchment = "CA") 
#   i <- i+1
#   
#   # With 2way clustering
#   reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_rapid_pixelcount"),
#                                                    SE = "twoway")
#   i <- i+1
#   
#   # negbin
#   # reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#   #                                                  outcome_variable = paste0("lucpfip_rapid_pixelcount"),
#   #                                                  distribution = "negbin")
#   # i <- i+1
#   
#   # unit FE
#   reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_rapid_pixelcount"),
#                                                    fe = "lonlat")
#   i <- i+1
# }


## For alternative distributional assumptions
for(DISTR in c("poisson")){#, "negbin"
  reg_stats_indvar_list[[paste0(DISTR,"_distribution")]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                                               outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                               distribution = DISTR)
  i <- i+1
}

## For alternative fixed effects
for(FE in c("reachable", "reachable + year", 
            "reachable + province_year", "reachable + subdistrict_year", "lonlat + district_year")){#, "lonlat + village_year"
  reg_stats_indvar_list[[paste0(FE," fe")]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                  fe = FE)
  i <- i+1
}
i  # 18
## For an alternative standard error computation (two-way clustering)
c <- 1
for(CLT in list("lonlat", "village", "subdistrict", "district", c("reachable","district_year"))){#"reachable" ,  
  reg_stats_indvar_list[[paste0("cluster_", c)]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                                       outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                       cluster = CLT)
  c <- c + 1
}


### For alternative control sets, (with their interaction each time)
# Without controls
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c(), 
#                                                  interaction_terms = NULL)

ov_lag1 <- paste0("lucpf",SIZE,"p_pixelcount_lag1")
ov_lag4 <- paste0("lucpf",SIZE,"p_pixelcount_lag4")

control_sets_list <- list(#c("n_reachable_uml"), # remoteness control only
  c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"), # Ownership control only
  c(ov_lag1), # lagged outcome variable
  c("wa_prex_cpo_imp1"), # prex cpo control only
  c("baseline_forest_trend"), # baseline forest trend only
  # c("remain_pf_pixelcount"), # remaining forest only
  
  c("n_reachable_uml",
    "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"), # remoteness and ownership only
  c("n_reachable_uml", 
    ov_lag1), # remoteness and lagged ov control                          
  c("n_reachable_uml", 
    "wa_prex_cpo_imp1"), # remoteness and wa_prex_cpo_imp1 control
  c("n_reachable_uml", 
    "baseline_forest_trend"), # remoteness and baseline_forest_trend control
  
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    ov_lag1), # Ownership, remoteness, and lagged ov control   
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    "wa_prex_cpo_imp1"), # Ownership, remoteness, and wa_prex_cpo_imp1 control
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    "baseline_forest_trend"), # Ownership, remoteness, and baseline_forest_trend control
  c("n_reachable_uml", ov_lag1, 
    "wa_prex_cpo_imp1"), # Remoteness, lagged ov, and export control
  c("n_reachable_uml", ov_lag1, 
    "baseline_forest_trend"), # Remoteness, lagged ov, and baseline_forest_trend control
  c("n_reachable_uml", "wa_prex_cpo_imp1", 
    "baseline_forest_trend"), # Remoteness, lagged ov, and baseline_forest_trend control
  # c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
  #   "remain_pf_pixelcount"), # Ownership, remoteness, and remain_pf_pixelcount control
  
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    "wa_prex_cpo_imp1", "baseline_forest_trend"), # all, controls, except lagged ov
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    ov_lag1, "wa_prex_cpo_imp1"), # all, controls, except baseline
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    ov_lag1, "baseline_forest_trend"), # all, controls, except export
  c("n_reachable_uml", ov_lag1,
    "wa_prex_cpo_imp1", "baseline_forest_trend"), # all, controls, except ownership 
  
  c("n_reachable_uml", "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", 
    "wa_prex_cpo_imp1", "baseline_forest_trend", ov_lag1) # all controls
)
length(control_sets_list) == length(unique(control_sets_list))

c <- 1
for(ctrl in control_sets_list){
  reg_stats_indvar_list[[paste0("ctrl_",c)]] <- make_spec_chart_df(island = ISL, illegal = ILL,
                                                                   outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
                                                                   controls = ctrl)
  c <- c + 1
}

# We don't feature these two specifications (lagged deforestation) anymore, as they are handled in a dedicated table. 
# ## Add the control set with the 4-year lagged deforestation 
# reg_stats_indvar_list[[paste0("ctrl_",c)]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                                  controls = c("n_reachable_uml",
#                                                                               ov_lag4))
# c <- c + 1
# 
# ## Add the control set with the 4-lagged deforestation in the neighbor grid cells
# reg_stats_indvar_list[[paste0("ctrl_",c)]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                                  controls = c("n_reachable_uml",
#                                                                               "ngb_ov_lag4"))
# c <- c + 1


# # two basic + wa_prex_cpo_imp1 and remain_pf_pixelcount
# ctrl <- c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
#           "wa_prex_cpo_imp1", "remain_pf_pixelcount")
# reg_stats_indvar_list[[paste0(ctrl, sep = "_")]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = ctrl)
# i <- i+1

# # two basic + baseline_forest_trend and remain_pf_pixelcount
# ctrl <- c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
#           "baseline_forest_trend", "remain_pf_pixelcount")
# reg_stats_indvar_list[[paste0(ctrl, sep = "_")]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = ctrl)
# i <- i+1

# ## All controls
# ctrl <- c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml",
#           "wa_prex_cpo_imp1", "baseline_forest_trend", "remain_pf_pixelcount")
# reg_stats_indvar_list[[paste0(ctrl, sep = "_")]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = ctrl)
# i <- i+1



## Now with less interactions
# # With all controls, but no interaction
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = NULL)
# i <- i+1
#   
# # With all controls, but only ownership interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp"))
# i <- i+1
# 
# # With all controls, but only remoteness interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("n_reachable_uml"))
# i <- i+1
# 
# # With all controls, but only wa_prex_cpo_imp1 interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_prex_cpo_imp1"))
# i <- i+1
# 
# # With all controls, but only ownership and remoteness interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml"))
# i <- i+1
# 
# # With all controls, but only ownership and wa_prex_cpo_imp1 interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "wa_prex_cpo_imp1"))
# i <- i+1
# 
# # With all controls, but only remoteness and wa_prex_cpo_imp1 interactions
# reg_stats_indvar_list[[i]] <- make_spec_chart_df(island = ISL, illegal = ILL,
#                                                  outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                  controls = c("wa_pct_own_nat_priv_imp","wa_pct_own_for_imp", "n_reachable_uml", "wa_prex_cpo_imp1"), 
#                                                  interaction_terms = c("n_reachable_uml", "wa_prex_cpo_imp1"))
# i <- i+1


# convert to dataframe to be able to chart
reg_stats_indvar <- bind_rows(reg_stats_indvar_list)
#reg_stats_indvar <- reg_stats_indvar[,-c(ncol(reg_stats_indvar))]

# save it 
if(sum(duplicated(reg_stats_indvar))==0 ){ # i.e. 50 currently & nrow(reg_stats_indvar)+1 == i
  saveRDS(reg_stats_indvar, file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",SIZE,"_",ILL,"_18072025")))
} else{print(paste0("SOMETHING WENT WRONG in spec_chart_df_",ISL,"_",SIZE,"_",ILL))}

#}else{reg_stats_indvar <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",ISL)))}
#   }
# }

# one_25 <- reg_stats_indvar[,c(1:25)] 
# nearest_mill <- one_25[,"nearest_mill"]
# one_25 <- one_25[,-4]
# one_25[,"nearest_mill"] <- nearest_mill
# new <- cbind(one_25, reg_stats_indvar[,c(26:34)])
# reg_stats_indvar <- new

# reg_stats_indvar <- dplyr::filter(reg_stats_indvar, nearest_mill == FALSE)
# reg_stats_indvar <- dplyr::select(reg_stats_indvar, -nearest_mill)

### PLOTTING 
### GIVE HERE THE ISLAND, THE OUTCOME AND THE DATE FOR WHICH YOU WANT THE SPEC CHART TO BE PLOTTED
scdf <- readRDS(file.path(paste0("temp_data/reg_results/spec_chart_df_",ISL,"_",SIZE,"_",ILL,"_18072025")))

# some modifications for now because scdf run with some "mistakes"
# scdf <- dplyr::select(scdf, -weights)
# scdf <- dplyr::mutate(scdf, own_interact = (own_for_interact | own_nat_priv_interact))
# scdf <- dplyr::select(scdf, -own_for_interact, -own_nat_priv_interact)
# scdf <- scdf[,c(1:23,29,24:28)]
# 
# scdf$ngb_ov_lag4 <- FALSE
# scdf <- cbind(scdf[,1:19], scdf[,ncol(scdf)], scdf[,20:31])
# names(scdf)[20] <- "ngb_ov_lag4"
# 
# colnames(scdf) %>% all.equal(colnames(reg_stats_indvar))
# scdf <- rbind(scdf, reg_stats_indvar)

# make labels for the chart
schart_labels <- list(#"Dependent variable:" = c("Larger forest definition"),
  "Data cleaning:" = c("Stronger imputations", 
                       "Lagged IBS variables"),
  
  "Sampling:" = c("Minimum forest cover in 2000",
                  "Minimum IBS to UML mill ratio"),
  
  "Plantation site size" = c("1x1 km",
                             "3x3 km", 
                             "5x5 km"),
  
  "Catchment model" = c(#"Nearest mill",
    #"Intensive margin expansion only",
    "Alternative catchment radius",
    "2-hour driving catchment area"), 
  
  "Price signal averaged over:" = c("1 year",
                                    "3 years", 
                                    "4 years", 
                                    "5 years"),        
  
  "Distribution assumption:" = c("Poisson",
                                 "Quasi-Poisson", 
                                 "Negative Binomial"),
  
  
  "Controls:" = c(#paste0(toupper(c("ffb", "cpo")[!grepl(VAR, c("ffb", "cpo"))]), " price signal"), 
    "# reachable mills", 
    "Ownership",
    "1-y lag deforestation",
    # "4-y lag deforestation",
    # "4-y neighbors' deforestation",
    "% CPO exported", 
    "Baseline forest trend"),#"Remaining forest""Neighbors' outcomes, 4-year lagged"
  # "Interaction with:" = c(#paste0(toupper(c("ffb", "cpo")[!grepl(VAR, c("ffb", "cpo"))]), " price signal"), 
  #   "Ownership",
  #   "# reachable UML mills", 
  #   "% CPO exported", 
  #   "Remaining forest"),
  #"Weights" = "", 
  
  "Fixed effects:" = c("Set of reachable mills", 
                       "Set of reachable mills and year", 
                       "Set of reachable mills and province-year",
                       "Set of reachable mills and district-year", 
                       "Set of reachable mills and subdistrict-year", 
                       "Plantation site and district-year"),
  
  "Level of SE clustering:" = c("Plantation", 
                                "Set of reachable mills",
                                "Village", 
                                "Subdistrict", 
                                "District",
                                "Set of reachable mills and district-year")
) 

# # remove rows where SEs could not be computed
# problems <- reg_stats_indvar[is.na(reg_stats_indvar$`Std. Error`),]
# reg_stats_indvar <- reg_stats_indvar[!is.na(reg_stats_indvar$`Std. Error`),]
# rownames(reg_stats_indvar) <- seq(1, nrow(reg_stats_indvar))


# find position of model to highlight in original data frame
a <- scdf
a <- a[#a$larger_forest_def==FALSE & 
  #a$rapid == FALSE &
  a$imp1 &
    a$lag_or_not &
    
    a$min_forest_2000 == FALSE & 
    a$min_coverage == FALSE & 
    
    a$PS_3km & 
    
    #a$nearest_mill == FALSE & 
    #a$intensive_only == FALSE & 
    a$alt_cr == FALSE &
    a$CA == FALSE &
    
    a$only_sr == FALSE &
    a$ya_4 &
    
    #a$CR_30km &  
    a$quasipoisson &
    
    a$unit_distryear_fe &
    #a$two_commo == FALSE & 
    
    a$n_reachable_uml_control  &
    a$control_own == FALSE & 
    a$ov_lag1 == FALSE &
    # a$ov_lag4 == FALSE &
    # a$ngb_ov_lag4 == FALSE & 
    a$prex_cpo_control == FALSE &
    a$baseline_forest_trend == FALSE &
    # a$remaining_forest == FALSE & 
    # a$ngb_ov_lag4 == FALSE & 
    # a$own_interact & 
    # a$n_reachable_uml_interact & 
    # a$prex_cpo_interact & 
    # a$remaining_forest_interact == FALSE & 
    #a$weights == FALSE & 
    a$reachable_cluster, ] 

model_idx <- row.names(a)
# this is our baseline model(s)
a[model_idx,]

# are there duplicated estimates (this is not expected) 
scdf[duplicated(scdf[,c(1,2)]),] # Poisson may be exactly the same as quasipoisson at the level of rounding CIs

schart(scdf, 
       labels = schart_labels,
       order="increasing",
       heights=c(.6,1.5),
       pch.dot=c(20,20,20,20),
       ci=c(.95),
       cex=c(0.6,0.7),
       highlight=model_idx, 
       #col.est=c(rgb(0,0,0,0.1),rgb(0,0.2,0.6, 0.1)),
       col.est=c("grey70", "red3"),
       #col.est2=c(rgb(0,0,0,0.08),"lightblue"),
       #col.dot=c(rgb(0,0,0,0.12),"grey95","grey95",rgb(0,0.4,0.6,0.3)),
       col.dot=c("grey70","grey95","grey95","red3"),
       #bg.dot=c(rgb(0,0,0,0.12),"grey95","grey95",rgb(0,0.4,0.6,0.3)),
       bg.dot=c("grey60","grey95","grey95","white"),
       adj=c(1,1), 
       #offset=c(10.4,10.1),
       leftmargin = 12,
       ylab = "Elasticity estimates",
       lwd.est = 5.8,
       lwd.symbol = 1
       #pch.est=20
)

# Should be saved in .pdf landscape A4 for correct vizualization under Latex.  





# COUNTERFACTUAL APPLICATION ---------------------------------------------------

# Produce the approximate population sample

# same sample selection but with some enlargements: 
# those within CR (intensive margin) of UML mills, 
# with positive deforestation at least one year, 
# with NAs in any covariate, 
# with remaining positive forest extent (but that should not change the cross-section size) 
# outside RSPO. 
# Rather annual than full period. This means that from the above selection, we take the number of grid cells. 
# Or the annual average number of grid cells. The former is better in line with the idea of effect on Indonesian oil palm sector as it was at the end of the study period. 


# Grid cells within 50km of a known (UML) mill. 
lucpfip <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",
                                    parcel_size/1000,"km_",
                                    "82km_UML_CR.rds")))
lucpfsmp <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfsmp_panel_",
                                     parcel_size/1000,"km_",
                                     "82km_UML_CR.rds")))
lucpfip_dyn <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_dynamics_",
                                        parcel_size/1000,"km_",
                                        "82km_UML_CR.rds")))

# keep only year before 2015 (after they mean nothing since we plantation data are from 2015)
# and after 2001 (because the d_clean final sample used for estimation runs from 2002 due 4-year average variables being only observed from 2002 on.)
lucpfip <- lucpfip[lucpfip$year < 2015 & lucpfip$year > 2001,] 
lucpfip_dyn <- lucpfip_dyn[lucpfip_dyn$year < 2015 & lucpfip_dyn$year > 2001,] 
lucpfsmp <- lucpfsmp[lucpfsmp$year < 2015 & lucpfsmp$year > 2001,] 

# remove coordinates, except in one data frame
lucpfsmp <- dplyr::select(lucpfsmp, -lat, -lon, -idncrs_lat, -idncrs_lon)
lucpfip_dyn <- dplyr::select(lucpfip_dyn, -lat, -lon, -idncrs_lat, -idncrs_lon)

# join them
d <- inner_join(lucpfip, lucpfip_dyn, by = c("lonlat", "year")) 
d <- inner_join(d, lucpfsmp, by = c("lonlat", "year")) 
rm(lucpfip, lucpfsmp, lucpfip_dyn)

## make a variable that counts rapid and slow lucpfp events (this is only computed for primary forest as of now)
# THIS IS GOING TO BE THE MAIN OUTCOME VARIABLE INSTEAD OF lucpfip_pixelcount_total
d$lucpfip_pixelcount <- d$lucpfip_rapid_pixelcount + d$lucpfip_slow_pixelcount

# make variable that counts lucpfp events on both small and medium sized plantations 
d$lucpfsmp_pixelcount <- d$lucpfsp_pixelcount_total + d$lucpfmp_pixelcount_total


## make a variable that counts lucpfp events on all types of plantations
d <- mutate(d, lucpfap_pixelcount = lucpfip_pixelcount + lucpfsmp_pixelcount) 


## MAKE THE lonlat 
uni_lonlat <- unique(d$lonlat)
d <- mutate(d, 
            lonlat = match(lonlat, uni_lonlat)) 
rm(uni_lonlat)

# some arrangements
d <- dplyr::arrange(d, lonlat, lonlat, year)
row.names(d) <- seq(1,nrow(d))


### RESTRICT THE GRID CELL SAMPLE 

## to islands of interest 
# d has no island variable so far. But it is the merger of prepare_ scripts run over Sumatra and Kalimantan, which are the islands of interest. 



## to within 50km a UML mill
uml <- read.dta13(file.path("temp_data/processed_UML/UML_valentin_imputed_est_year.dta"))
uml <- uml[!is.na(uml$lat),]
uml <- st_as_sf(uml, coords = c("lon", "lat"), remove = TRUE, crs = 4326)
uml <- st_transform(uml, crs = indonesian_crs)
uml <- st_buffer(uml, dist = 5e4)
uml <- st_union(uml)
# uml <- st_as_sf(uml)
# uml$CR <- TRUE

d <- st_as_sf(d, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs)

d_cs <- d[!duplicated(d$lonlat), ]
sgbp <- st_within(d_cs, uml)
d_cs$CR <- rep(FALSE, nrow(d_cs))
d_cs$CR[lengths(sgbp) > 0] <- TRUE
sum(d_cs$CR)

d_cs <- st_drop_geometry(d_cs)
d <- left_join(d, d_cs[, c("lonlat", "CR")], by = "lonlat")
d_save <- d
#d <- filter(d, d$CR)


## With some deforestation at least once 
# # If it's smallholders
# temp_est_removeObs <- feglm(fml = as.formula(paste0("lucpfsmp_pixelcount ~ lonlat")),
#                             data = d, 
#                             family = "poisson")
# 
# d$smallholders <- (d$lonlat %in% d[unlist(temp_est_removeObs$obs_selection),]$lonlat )
# rm(temp_est_removeObs)
# d$smallholders <- (d$lonlat %in% d[-obs2remove(fml = as.formula(paste0("lucpfsmp_pixelcount ~ lonlat")),
#                                                   data = d, 
#                                                   family = "poisson"),]$lonlat )


temp_est <- feglm(fml = as.formula(paste0("lucpfap_pixelcount ~ 1 | lonlat")),
                  data = d, 
                  family = "poisson")

d <- d[unlist(temp_est$obs_selection),] 
rm(temp_est)

# d <- d[-obs2remove(fml = as.formula(paste0("lucpfap_pixelcount ~ lonlat")),
#                    data = d, 
#                    family = "poisson"),]



## With some variables possibly NAs - ... need not do anything 



## That are not in RSPO
rspo <- st_read("input_data/RSPO_supply_bases/RSPO-certified_oil_palm_supply_bases_in_Indonesia.shp")
rspo <- st_transform(rspo, crs = indonesian_crs)
# set up a year variable $
# we assume that the year of probable behavior change is the year after the date of the icletdate variable
# BUT THIS MAY CHANGE 
rspo$year <- sub(pattern = ".*, ", 
                 replacement = "", 
                 x = rspo$icletdate)
rspo$year[rspo$year=="no letter"] <- "2015"
rspo$year <- rspo$year %>% as.integer()
rspo$year <- rspo$year +1

d$rspo_cert <- rep(FALSE, nrow(d))

n_with_rspo <- nrow(d)
gc_with_rspo <- length(unique(d$lonlat))

for(y in min(rspo$year):max(d$year)){
  # select d from the given year
  d_cs <- d[d$year == y, c("lonlat", "year", "rspo_cert")]
  # select supply bases already certified this year
  rspo_cs <- rspo[rspo$year <= y,] %>% st_geometry()
  
  sgbp <- st_within(x = d_cs, y = rspo_cs)
  
  d$rspo_cert[d$year == y][lengths(sgbp) == 1] <- TRUE
}

# - are not in an intended RSPO-certified supply base
d <- d[d$rspo_cert==FALSE,]

# how many in RSPO ? It removes observations but no grid cell, because all RSPO grid cells have obs. 
# before RSPO started
n_with_rspo - nrow(d)
gc_with_rspo - length(unique(d$lonlat))



## And finally identify illegal deforestation in this population, but do not restrict the sample to it yet. 
# OIL PALM CONCESSIONS
cns <- st_read(file.path("input_data/oil_palm_concessions"))
cns <- st_transform(cns, crs = indonesian_crs)

# LEGAL LAND USE 
llu <- st_read(file.path("input_data/kawasan_hutan/Greenorb_Blog/final/KH-INDON-Final.shp"))
llu <- st_transform(llu, crs = indonesian_crs)
unique(llu$Fungsi)
names(llu)[names(llu) == "Fungsi"] <- "llu"

# restrict llu to provinces of interest
llu <- llu[llu$Province == "Sumatra Utara" |
             llu$Province == "Riau" |
             llu$Province == "Sumatra Selantan" |
             llu$Province == "Papua Barat" |
             llu$Province == "Kalimantan Timur" |
             llu$Province == "Kalimantan Selatan" |
             llu$Province == "Kalimantan Tengah" |
             llu$Province == "Kalimantan Barat" |
             llu$Province == "Bengkulu" |
             llu$Province == "Lampung" |
             llu$Province == "Jambi" |
             llu$Province == "Bangka Belitung" |
             llu$Province == "Kepuluan Riau" |
             llu$Province == "Sumatra Barat" |
             llu$Province == "Aceh", ]


# OIL PALM CONCESSIONS
# We do not observe whether a grid cell is within a concession annually. 
# Therefore we only proceed with a cross section
d_cs <- d[!duplicated(d$lonlat),]
sgbp <- st_within(d_cs, cns)
d_cs$concession <- rep(FALSE, nrow(d_cs))
d_cs$concession[lengths(sgbp) > 0] <- TRUE

d_cs <- st_drop_geometry(d_cs)

d <- left_join(d, d_cs[,c("lonlat", "concession")], by = "lonlat")

# note that some d fall within more than one concession record. There may be several reasons for concession overlaps 
# like renewal of concession, with our withour aggrandisement. For our purpose, it only matters that there is at least one 
# concession record. 


# LEGAL LAND USE 
# this is quite long (~5min)
d_cs <- d[!duplicated(d$lonlat),] # the point of this is to give back spatial class to d_cs
d_cs <- st_join(x = d_cs, 
                y = st_make_valid(llu[,"llu"]), # st_make_valid bc thrown error otherwise 
                join = st_within, 
                left = TRUE)

# some grid cells seem to fall within overlapping llu shapes though. 
# It's really marginal (12 instances). Just remove the duplicates it produces.
d_cs <- d_cs[!duplicated(d_cs$lonlat),]

# merge back with panel 
d <- st_drop_geometry(d)
d_cs <- st_drop_geometry(d_cs)

d <- left_join(d, d_cs[,c("lonlat", "llu")], by = "lonlat")

unique(d$llu)

# ILLEGAL LUCFP 
# one possible link to shed light on accronyms http://documents1.worldbank.org/curated/pt/561471468197386518/pdf/103486-WP-PUBLIC-DOC-107.pdf

d <- dplyr::mutate(d,
                   illegal1 = (!concession & (llu != "HPK" )), # it's not in concession and not in a convertible forest zone. Don't add the following code, because these NAs are for all places outside the forest estate, and it changes exactly nothing to add this condition | llu == "<NA>"
                   illegal2 = (!concession & (llu == "HL" | # it's not in concession and it's in a permanent forest zone designation
                                                
                                                llu == "HP" | # production forest : " these areas may be selectively logged in a normal manner".
                                                llu == "HPT" | # limited production forest : "These areas be logged less intensively than is permitted in the Permanent Production Forest" 
                                                
                                                llu=="HK" | # below are all categories of HK
                                                llu=="KSA/KPA" |
                                                llu=="KSA" | 
                                                llu=="CA" | 
                                                llu=="SM" | 
                                                llu=="KPA" | 
                                                llu=="TN" | 
                                                llu=="TWA" | 
                                                llu=="Tahura" | 
                                                llu=="SML" | 
                                                llu=="CAL" | 
                                                llu=="TNL" | 
                                                llu=="TWAL" | 
                                                llu=="KSAL" | 
                                                llu=="TB" | 
                                                llu=="Hutan Cadangan")))

# yields many missing in illegal because many grid cells are within a mising land use legal classification
# parcels[!duplicated(parcels$lonlat) & !is.na(parcels$llu), c("lonlat", "concession", "llu", "illegal1", "illegal2")]

# restrict the sample to the pre-defined legal status 






### COUNT THE AVERAGE ANNUAL NUMBER OF GRID CELLS ###
aggr_factor_all <- ddply(d, "year", summarise,
                         annual_avg_n = length(unique(lonlat)))[,"annual_avg_n"] %>% mean() %>% round(0)

# aggr_factor_sm <- ddply(d[d$smallholders, ], "year", summarise,
#                         annual_avg_n = length(unique(lonlat)))[,"annual_avg_n"] %>% mean() %>% round(0)
# 
# aggr_factor_ill <- ddply(d[!is.na(d$illegal2) & d$illegal2 == TRUE, ], "year", summarise,
#                          annual_avg_n = length(unique(lonlat)))[,"annual_avg_n"] %>% mean() %>% round(0)
# because of the many NAs in illegal variables, we do not know the legal status for a large proportion of griod cells. 
# therefore we cannot say that "of the total aggregated effect, a given amount comes from illegal deforestation". 

### COMPUTE THE AVERAGE PARTIAL EFFECT ###




res_data_all <- make_base_reg(island = "both",
                              outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
                              output_full = FALSE)

# res_data_sm <- make_base_reg(island = "both",
#                              outcome_variable = paste0("lucpfsmp_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                              output_full = FALSE)
# 
# res_data_ill <- make_base_reg(island = "both",
#                               outcome_variable = paste0("lucpfap_pixelcount"), # or can be  lucpf",SIZE,"p_pixelcount"
#                               ill = "ill2",
#                               output_full = FALSE)


make_counterfactuals <- function(res_data, aggr_factor){
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it is necessary that d_clean is in the main environment so that vcov.fixest can retrieve cluster variable
  coeff <- reg_res$coefficients[1]
  names(coeff) <- NULL
  
  mean(d_clean$lucpfap_pixelcount)*pixel_area_ha # 5.172615 ha
  # this is exactly equal to the average fitted values
  mean(reg_res$fitted.values)*pixel_area_ha
  
  
  # but because reachable and district_year variables were not available to compute the aggregating factor, 
  # we used a different restriction on always null (only on lonlat), which is more restrictive. 
  # as a consequence, the scaling factor is too small compared to the average deforestation, 
  # so we adjust by taking the average outcome value in the analysis sample - additionally restricted. 
  temp_est <- feglm(fml = as.formula(paste0("lucpfap_pixelcount ~ 1 | lonlat")),
                    data = d_clean, 
                    family = "poisson")
  
  d_clean <- d_clean[unlist(temp_est$obs_selection),] 
  rm(temp_est)
  
  avg_ha <- mean(d_clean$lucpfap_pixelcount)*pixel_area_ha # 11.75246 ha
  
  print(paste0("Baseline deforestation is ", avg_ha*aggr_factor, " ha"))
  # avg_ha*aggr_factor
  
  # one std. dev. in price 
  # remove fixed effect variations from the regressor
  reg_sd <- fixest::feols(fml = as.formula(paste0(
    names(reg_res$coefficients)[1],
    " ~ 1 | ", 
    paste0(reg_res$fixef_vars, collapse = "+"))),
    data = d_clean)
  
  # and take the remaining standard deviation
  sd_price_change <- sd(reg_sd$residuals)
  print(paste0("std. dev. price change is: ",sd_price_change*100, "%"))
  
  # Price decrease necessary to halve deforestation (in pct.)
  delta_defo <- -0.29 # Indonesian unconditional target in NDC = 29% below BAU by 2030 incl. LULUCF
  tax <- (((1 + delta_defo)^(1/coeff)) - 1)
  
  # store application results in table
  price_changes <- c(sd_price_change, -0.01, tax)
  app_mat <- matrix(nrow = 3, ncol = length(price_changes))
  row.names(app_mat) <- c("Relative change (%)", "Total change (ha)", "Potential CO2 revenues (M$)")
  colnames(app_mat) <- as.character(price_changes)
  
  for(rel_price_change in price_changes){
    effect <- ((1+rel_price_change)^(coeff) - 1)
    app_mat[1, match(rel_price_change, price_changes)] <- effect*100
    
    defo <- (effect)*avg_ha*aggr_factor
    app_mat[2, match(rel_price_change, price_changes)] <- defo
    
    app_mat[3, match(rel_price_change, price_changes)] <- -defo*174*(44/12)*5*1e-6
    # 174 is the Carbon loss associated with deforestation, in ton (Guillaume et al. 2018)
    # 44/12 is the C to CO2 conversion factor 
    # 5 is the CO2 price paid by Norway in the first REDD+ payment in $/tCO2eq. https://www.regjeringen.no/en/aktuelt/noreg-betaler-530-millionar-for-redusert-avskoging-i-indonesia/id2722135/
    # We give the final value in million $
  }
  
  app_mat[1,] <- app_mat[1,] %>% round(2)
  app_mat[2,] <- app_mat[2,] %>% formatC(digits = 0, format = "f")
  app_mat[3,] <- as.numeric(app_mat[3,]) %>% formatC(digits = 1, format = "f")
  
  return(app_mat) 
}
cf_mat_all <- make_counterfactuals(res_data_all, aggr_factor = aggr_factor_all)
# cf_mat_sm <- make_counterfactuals(res_data_sm, aggr_factor = aggr_factor_sm)
# cf_mat_ill <- make_counterfactuals(res_data_ill, aggr_factor = aggr_factor_ill)

cf_mat <- cbind(cf_mat_all) # cf_mat_sm, cf_mat_ill, 
# 
# # opportunity cost for the oil palm sector associated with a 1% decrease in prices
# ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
# annual_avg_cpo_price <- mean(ibs$cpo_price_imp1, na.rm= T)
# annual_avg_cpo_production <- 19963.57143 * 1000 # between 2001 and 2014, according to USDA (which gives it in 1000 MT): https://www.indexmundi.com/agriculture/?country=id&commodity=palm-oil&graph=production
# cpo_oc <- 0.01 * annual_avg_cpo_price * annual_avg_cpo_production


# ACTUALLY IT IS NOT RELEVANT TO DISPLAY ILLEGAL AND SM, 
# because aggr_factor_ill is underestimated because many missings on illegal, so the participation of illegal to the total change
# is not meaningfull. 
cf_mat <- cf_mat_all
colnames(cf_mat) <- NULL
options(knitr.table.format = "latex")
kable(cf_mat, booktabs = T, align = "c",
      caption = "Counterfactual annual effects of different CPO price changes on deforestation in Indonesia") %>% #of 1 percentage change in medium-run price signal
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  add_header_above(c(" " = 1,
                     "+1 std. dev." = 1, # 
                     "-1%" = 1,
                     "-19%" = 1
                     # "+1 std. dev." = 1, # 
                     # "-19%" = 1,
                     # "+1 std. dev." = 1, # 
                     # "-19%" = 1
  ),
  bold = F,
  align = "c") %>%
  
  column_spec(column = 1,
              width = "9em",
              latex_valign = "m") %>% 
  column_spec(column = c(2:(ncol(cf_mat))),
              width = "5em",
              latex_valign = "m")

rm(cf_mat)
