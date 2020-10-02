### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
#             EXPLICATIVE VARIABLE DESCRIPTIVE STATISTICS 
# 
#   input   - panel data frame of parcels with the information on number of reachable uml mills
#           as outputted from make_n_reachable_uml.R
#           ---> temp_data/processed_parcels/wa_panel_parcels_reachable_uml_PS_CR.rds
# 
#           - panel data frame of IBS 
#           ---> temp_data/IBS_UML_panel_final.dta
# 
# 
# 
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

##### 0. PACKAGES, WD, OBJECTS #####


### WORKING DIRECTORY SHOULD BE CORRECT IF THIS SCRIPT IS RUN WITHIN R_project_for_individual_runs
### OR CALLED FROM LUCFP PROJECT master.do FILE.
### IN ANY CASE IT SHOULD BE (~/LUCFP/data_processing) 


### PACKAGES ###
# see this project's README for a better understanding of how packages are handled in this project. 

# These are the packages needed in this particular script. *** these are those that we now not install: "rlist","lwgeom","htmltools", "iterators", 
neededPackages = c("plyr", "dplyr", "data.table",  "stringr", "sjmisc", 
                   "foreign", "readxl", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr")
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
troublePackages <- c("kableExtra") 
# Attempt to load packages from user's default libraries.
lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...") 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 

### NEW FOLDERS USED IN THIS SCRIPT 

### INDONESIAN CRS ### 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator.
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### PARCEL SIZE ###
parcel_size <- 3000


##### IBS MILL PANEL STAT DES ##### 
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



# We want to produce, for a set of ibs variables, the mean, median, std.dev. min and max statistics, 
# for the geo_sample that we are going to use, and the is_mill sample. 

make_des_table_island <- function(ibs_isl){
  
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
  row.names(rhs_des) <- c("First year in IBS", "FFB muv (USD/ton)", "FFB input (ton)", 
                          "CPO muv (USD/ton", "CPO output (ton)", 
                          "PKO muv (USD/ton)", "PKO output (ton)", 
                          "CPO export share (%)", 
                          "Central government ownership (%)", 
                          "Local government ownership (%)", 
                          "National private ownership (%)", 
                          "Foreign ownership (%)")
  
  return(rhs_des)
}

#### Print the LateX table code ALL #### 
ibs_isl <- ibs
des_table <- make_des_table_island(ibs_isl)
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

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






#### Print the LateX table code SUMATRA #### 

ibs_isl <- ibs[ibs$island_name=="Sumatra",]
des_table <- make_des_table_island(ibs_isl)

# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, Sumatra") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "p-value" = 1,
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 464 mills" = 3,
                     "All IBS palm oil mills \n n = 677 mills" = 3, 
                     "t-test" = 1, 
                     "KS test" =1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "4em") %>% 
  column_spec(column = c(4,7),
              width = "10em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

#### Print the LateX table code KALIMANTAN #### 

ibs_isl <- ibs[ibs$island_name=="Kalimantan",]
des_table <- make_des_table_island(ibs_isl)
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, Kalimantan") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "p-value"=1,
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 98 mills" = 3,
                     "All IBS palm oil mills \n n = 200 mills" = 3, 
                     "t-test"=1, 
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


#### Print the LateX table code OTHERS #### 
ibs_isl <- ibs[!(ibs$island_name %in% c("Sumatra", "Kalimantan")),]
des_table <- make_des_table_island(ibs_isl)
# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  
length(unique(ibs_isl[ibs_isl$analysis_sample==TRUE, "firm_id"]))
length(unique(ibs_isl[ibs_isl$is_mill==TRUE, "firm_id"]))

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "IBS descriptive statistics, Maluku, Papua, Sulawesi") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "p-value" = 1,
                     "p-value" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Geo-localized IBS palm oil mills \n n = 25 mills" = 3,
                     "All IBS palm oil mills \n n = 53 mills" = 3, 
                     "t-test" = 1,
                     "KS test" = 1),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "4em") %>% 
  column_spec(column = c(4,7),
              width = "10em") %>% 
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

##### IBS-UML COMPARATIVE STATS ##### 
# all IBS oil palm related establishments. 
ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))


# UML 
##### PARCEL STAT DES ##### 


# We want to produce, for a set of variables, the mean, median, std.dev. min and max statistics, 
# for the three samples of 3km parcels within 10, 30 and 50km. 



make_des_table_parcels_island <- function(island){
  
  CR_list <- list()  
  catchment_radiuseS <- c(1e4, 3e4, 5e4)
  
  for(catchment_radius in catchment_radiuseS){
    
    parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                        parcel_size/1000,"km_",
                                        catchment_radius/1000,"CR.rds")))
    
    # filter parcels to island
    parcels_isl <- parcels[parcels$island %in% island,]
    rm(parcels)
    # filter parcels to primary forested
    parcels_isl <- parcels_isl[parcels_isl$any_fc2000_30th,]
  
    variables <- c("lucpfip_ha_intact", "lucpfip_ha_total",
                 "lucfip_ha_90th", "lucfip_ha_60th", "lucfip_ha_30th", 
                 "wa_est_year_imp", 
                 "wa_ffb_price_imp1",
                 "wa_cpo_price_imp1", 
                 "wa_pko_price_imp1",
                 "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                 #paste0("wa_concentration_",catchment_radius/1000), 
                 paste0("n_reachable_uml"),
                 paste0("sample_coverage"))
    
    statistics <- c("mean", "std.dev", "median", "min", "max")
    
    ## Matrix for 
    des_CR <- matrix(NA, nrow = length(variables), ncol = length(statistics))
    row.names(des_CR) <- variables
    colnames(des_CR) <- c(statistics)
    
    for(var in variables){
      des_CR[var,statistics] <- summarise(parcels_isl,
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
      des_CR[var, "median"] <- paste0(des_CR[var,"median"]," [",des_CR[var, "min"],"; ",des_CR[var,"max"],"]")
    }
    
    des_CR <- des_CR[,c("mean", "std.dev", "median")]
    
    CR_list[[match(catchment_radius, catchment_radiuseS)]] <- des_CR 
    length(unique(parcels_isl$parcel_id)) %>% print()
  }
  
  des_parcels <- cbind(CR_list[[1]], CR_list[[2]], CR_list[[3]])
  rm(parcels_isl, CR_list)
  
  # row names
  row.names(des_parcels) <- c("From intact primary forest (ha)",
                              "From total primary forest (ha)",
                              "From 90% tree cover forest (ha)",
                              "From 60% tree cover forest (ha)",
                              "From 30% tree cover forest (ha)",
                              "Establishment year", 
                              "FFB price signal (USD/ton)", 
                              "CPO price signal (USD/ton", 
                              "PKO price signal (USD/ton)", 
                              "Central government ownership (%)", 
                              "Local government ownership (%)", 
                              "National private ownership (%)", 
                              "Foreign ownership (%)", 
                              #"Competition", 
                              "# UML mills in catchment radius", 
                              "% of which are matched with IBS")
  
  return(des_parcels)
}


#### Print the LateX table code ALL ISLANDS ####

des_table <- make_des_table_parcels_island(c("Sumatra", "Kalimantan", "Papua"))

# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Grid cell descriptive statistics, Sumatra, Kalimantan, Papua") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1),
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Within 10km catchment radius \n of geo-localized IBS mills \n n = 11289 cells" = 3,
                     "Within 30km catchment radius \n of geo-localized IBS mills \n n = 41844 cells" = 3,
                     "Within 50km catchment radius \n of geo-localized IBS mills \n n = 66769 cells" = 3),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "3em") %>%
  column_spec(column = c(4,7,10),
              width = "9em") %>%
  pack_rows("LUC to industrial oil palm plantations", 1, 5, 
            italic = TRUE)  %>%
  pack_rows("Invert-distance weighted averages \n of mills in catchment radius", 6, 13, #14 
            italic = TRUE)  %>%
  pack_rows("Grid cell features", 14,15, #15, 16 with competition variable 
            italic = TRUE)  %>%
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 




#### Print the LateX table code SUMATRA####

des_table <- make_des_table_parcels_island("Sumatra")

# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Grid cell descriptive statistics, Sumatra") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Within 10km catchment radius \n of geo-localized IBS mills \n n = 8470 cells" = 3,
                     "Within 30km catchment radius \n of geo-localized IBS mills \n n = 27512 cells" = 3,
                     "Within 50km catchment radius \n of geo-localized IBS mills \n n = 39786 cells" = 3),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "3em") %>%
  column_spec(column = c(4,7,10),
              width = "9em") %>%
  pack_rows("LUC to industrial oil palm plantations", 1, 5, 
            italic = TRUE)  %>%
  pack_rows("Invert-distance weighted averages \n of mills in catchment radius", 6, 13, 
            italic = TRUE)  %>%
  pack_rows("Grid cell features", 14,15, 
            italic = TRUE)  %>%
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

#### Print the LateX table code KALIMANTAN ####

des_table <- make_des_table_parcels_island("Kalimantan")

# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Grid cell descriptive statistics, Kalimantan") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Within 10km catchment radius \n of geo-localized IBS mills \n n = 2609 cells" = 3,
                     "Within 30km catchment radius \n of geo-localized IBS mills \n n = 12706 cells" = 3,
                     "Within 50km catchment radius \n of geo-localized IBS mills \n n = 23257 cells" = 3),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "3em") %>%
  column_spec(column = c(4,7,10),
              width = "9em") %>%
  column_spec(column = c(5,6,7),
              strikeout = TRUE) %>% #background = "#BCD4E6"
  pack_rows("LUC to industrial oil palm plantations", 1, 5, 
            italic = TRUE)  %>%
  pack_rows("Invert-distance weighted averages \n of mills in catchment radius", 6, 13, 
            italic = TRUE)  %>%
  pack_rows("Grid cell features", 14,15, 
            italic = TRUE)  %>%
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 

#### Print the LateX table code PAPUA ####

des_table <- make_des_table_parcels_island("Papua")

# we cannot automate headers, with a paste0, it does not work, hence n mills manually specified...  

colnames(des_table) <- NULL

options(knitr.table.format = "latex") 
kable(des_table, booktabs = T, align = "c", 
      caption = "Grid cell descriptive statistics, Papua") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1, 
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1,
                     "mean" = 1, "std.dev." = 1, "median [min; max]" = 1), 
                   align = "c", 
                   strikeout = F) %>% 
  add_header_above(c(" " = 1,
                     "Within 10km catchment radius \n of geo-localized IBS mills \n n = 210 cells" = 3,
                     "Within 30km catchment radius \n of geo-localized IBS mills \n n = 1626 cells" = 3,
                     "Within 50km catchment radius \n of geo-localized IBS mills \n n = 3726 cells" = 3),
                   bold = T,
                   align = "c") %>%
  column_spec(column = c(2,3,5,6,8,9),
              width = "3em") %>%
  column_spec(column = c(4,7,10),
              width = "9em") %>%
  pack_rows("LUC to industrial oil palm plantations", 1, 5, 
            italic = TRUE)  %>%
  pack_rows("Invert-distance weighted averages \n of mills in catchment radius", 6, 13, 
            italic = TRUE)  %>%
  pack_rows("Grid cell features", 14,15, 
            italic = TRUE)  %>%
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 


#### Table to compare samples of grid cells within 10, 30 and 50 km catchment radii. #### 
# Interpretation
# In Sumatra and overall, the differences are significant. 
# In Kalimantan and Papua, the hypotheses that the samples are alike are more difficult to reject. 
# This may reflect the fact that in Sumatra we observe a larger share of all existing mills. 
# While in Kali and Papua, a larger share of the grid cells that we observe as within 50km 
# of a mill is actually within 10 or 30km of a mill we do not observe. 
# Therefore, in Kalimantan and Papua, the impression of similarity between areas close and away from mills 
# may be spuriously given by our non exhaustive observation set of mills. 

## TESTS BETWEEN 10-30km and between 30-50km catchment radii
make_tests_CR_samples <- function(island){
  parcels10 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                        parcel_size/1000,"km_10CR.rds")))
  parcels30 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                        parcel_size/1000,"km_30CR.rds")))
  parcels50 <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                        parcel_size/1000,"km_50CR.rds")))
  
  # filter parcels to island
  parcels10_isl <- parcels10[parcels10$island %in% island,]
  rm(parcels10)
  parcels30_isl <- parcels30[parcels30$island %in% island,]
  rm(parcels30)
  parcels50_isl <- parcels50[parcels50$island %in% island,]
  rm(parcels50)
  # filter parcels to primary forested
  parcels10_isl <- parcels10_isl[parcels10_isl$any_fc2000_30th,]
  parcels30_isl <- parcels30_isl[parcels30_isl$any_fc2000_30th,]
  parcels50_isl <- parcels50_isl[parcels50_isl$any_fc2000_30th,]
  
  
  variables <- c("lucpfip_ha_intact", "lucpfip_ha_total",
                 "lucfip_ha_90th", "lucfip_ha_60th", "lucfip_ha_30th", 
                 "wa_est_year_imp", 
                 "wa_ffb_price_imp1",
                 "wa_cpo_price_imp1", 
                 "wa_pko_price_imp1",
                 "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp", "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",
                 #paste0("wa_concentration_",catchment_radius/1000), 
                 paste0("n_reachable_uml"),
                 paste0("sample_coverage"))
  
  ## add the t-test column
  # t-test between 10CR and 30CR grid cell samples
  t_test10_30 <- matrix(NA, nrow = length(variables), ncol = 1) 
  row.names(t_test10_30) <- variables
  colnames(t_test10_30) <- "t test"
  
  for(var in variables){
    test <- t.test(x = parcels10_isl[,var],
                   y = parcels30_isl[,var], 
                   alternative = "two.sided", 
                   mu = 0, 
                   paired = F, 
                   var.equal = FALSE)
    
    t_test10_30[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  
  # t-test between 30CR and 50CR grid cell samples
  t_test30_50 <- matrix(NA, nrow = length(variables), ncol = 1) 
  row.names(t_test30_50) <- variables
  colnames(t_test30_50) <- "t test"
  
  for(var in variables){
    test <- t.test(x = parcels30_isl[,var],
                   y = parcels50_isl[,var], 
                   alternative = "two.sided", 
                   mu = 0, 
                   paired = F, 
                   var.equal = FALSE)
    
    t_test30_50[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  # interpretation: if the 95% CI includes 0, then the difference in means between two samples 
  # is not statistically different from 0. Hence, the two samples are "similar" in means 
  # with respect to the variable tested. 
  # (In other words, we cannot reject the null hypothesis that the difference is null
  # -i.e. the two samples are alike - with 95% confidence) 
  
  ## add the Smirnov test
  # KS test between 10CR and 30CR grid cell samples
  ks_test10_30 <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(ks_test10_30) <- variables
  colnames(ks_test10_30) <- "KS test"
  
  for(var in variables){
    test <- ks.test(x = parcels10_isl[,var],
                    y = parcels30_isl[,var], 
                    alternative = "two.sided", 
                    exact = FALSE)
    
    ks_test10_30[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  
  # KS test between 30CR and 50CR grid cell samples
  ks_test30_50 <- matrix(NA, nrow = length(variables), ncol = 1)
  row.names(ks_test30_50) <- variables
  colnames(ks_test30_50) <- "KS test"
  
  for(var in variables){
    test <- ks.test(x = parcels30_isl[,var],
                    y = parcels50_isl[,var], 
                    alternative = "two.sided", 
                    exact = FALSE)
    
    ks_test30_50[var,] <- test$p.value %>% formatC(digits = 3, format = "f")
  }
  
  # intepretation: if p-value < 0.05 we cannot reject with 95% confidence that 
  # two distributions are equal, implying that they are different. 
  
  tests_parcels <- cbind(t_test10_30, ks_test10_30, t_test30_50, ks_test30_50)
  
  rm(parcels10_isl, parcels30_isl, parcels50_isl)
  
  # row names
  row.names(tests_parcels) <- c("From intact primary forest (ha)",
                              "From total primary forest (ha)",
                              "From 90% tree cover forest (ha)",
                              "From 60% tree cover forest (ha)",
                              "From 30% tree cover forest (ha)",
                              "Establishment year", 
                              "FFB price signal (USD/ton)", 
                              "CPO price signal (USD/ton", 
                              "PKO price signal (USD/ton)", 
                              "Central government ownership (%)", 
                              "Local government ownership (%)", 
                              "National private ownership (%)", 
                              "Foreign ownership (%)", 
                              #"Competition", 
                              "# UML mills in catchment radius", 
                              "% of which are matched with IBS")
  

  return(tests_parcels)
}

tests_table_all <- make_tests_CR_samples(c("Sumatra", "Kalimantan", "Papua"))
tests_table_suma <- make_tests_CR_samples(c("Sumatra"))
tests_table_kali <- make_tests_CR_samples(c("Kalimantan"))
tests_table_papu <- make_tests_CR_samples(c("Papua"))

tests_table <- cbind(tests_table_all, tests_table_suma, tests_table_kali, tests_table_papu)

colnames(tests_table) <- NULL

options(knitr.table.format = "latex") 
kable(tests_table, booktabs = T, align = "c", 
      caption = "Comparisons of grid cell samples within different catchment radii") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  add_header_above(c(" " = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1,
                     "t-test" = 1,
                     "KS-test" = 1),
                   bold = F,
                   align = "c") %>%
  add_header_above(c(" " = 1,
                     "10km vs. 30km \n catchment radius" = 2,
                     "30km vs. 50km \n catchment radius" = 2,
                     "10km vs. 30km \n catchment radius" = 2,
                     "30km vs. 50km \n catchment radius" = 2,
                     "10km vs. 30km \n catchment radius" = 2,
                     "30km vs. 50km \n catchment radius" = 2,
                     "10km vs. 30km \n catchment radius" = 2,
                     "30km vs. 50km \n catchment radius" = 2),
                   bold = F, 
                   align = "c") %>% 
  add_header_above(c(" " = 1,
                     "All islands" = 4,
                     "Sumatra" = 4,
                     "Kalimantan" = 4,
                     "Papua" = 4), 
                   bold = T, 
                   align = "c") %>% 
  pack_rows("LUC to industrial oil palm plantations", 1, 5, 
            italic = TRUE)  %>%
  pack_rows("Invert-distance weighted averages \n of mills in catchment radius", 6, 13, #14 
            italic = TRUE)  %>%
  pack_rows("Grid cell features", 14,15, #15, 16 with competition variable 
            italic = TRUE)  %>%
  column_spec(column = c(1),
              width = "20em") %>%
  footnote(general = c("Note"),
           threeparttable = TRUE, 
           escape = TRUE) 











### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# OLD CODE to flag different categories of establishments. Now this is done in merging_geolocalization_works.do

# # make indicator of our sample of interest (geo-referenced mills) 
# ibs$geo_sample <- !is.na(ibs$lat)
# 
# # make indicator of being a mill 
# is_mill <- ddply(ibs, "firm_id", summarise, 
#                  is_mill = ((sum(out_ton_cpo, na.rm = TRUE) > 0 | 
#                             sum(out_val_cpo, na.rm = TRUE) > 0 | 
#                             sum(out_ton_pko, na.rm = TRUE) > 0 | 
#                             sum(out_val_pko, na.rm = TRUE) > 0 ) | 
#                            (sum(in_ton_ffb, na.rm = TRUE) > 0 | 
#                             sum(in_val_ffb, na.rm = TRUE) > 0 )))
# ibs <- merge(ibs, is_mill, by = "firm_id")
# rm(is_mill)
# # also, do not count mills that satisfy the above conditions but are in Java or Bali; 
# ibs[ibs$island_name == "Java" | ibs$island_name == "Bali Nusa Tenggara", "is_mill"] <- FALSE
# 
# ibs$analysis_sample <- (ibs$geo_sample & ibs$is_mill)
