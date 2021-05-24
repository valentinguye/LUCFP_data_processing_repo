
op <- raster("C:/Users/GUYE/Downloads/aopm_lower_version2/lower/2015_op.tif")

wp <- st_read("C:/Users/GUYE/Google Drive/indonesia_diverse/windprospecting/shape.shp")

wp <- st_read("C:/Users/GUYE/Google Drive/indonesia_diverse/windprospecting/shape.shp")

gwa <- raster("C:/Users/GUYE/Google Drive/indonesia_diverse/IDN_wind-speed_50m.tif")
)

CR <- 3e4
d_list <- list()
elm <- 1

parcel_size <- 3000 
catchment_radiuseS <- c(1e4, 3e4, 5e4) # (in meters)
for(CR in catchment_radiuseS){
 d <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                           parcel_size/1000,"km_",
                           CR/1000,"CR.rds")))
  
 # # THIS IS GOING TO BE THE MAIN OUTCOME VARIABLE INSTEAD OF lucpfip_pixelcount_total
 # d$lucpfip_pixelcount <- d$lucpfip_rapid_pixelcount + d$lucpfip_slow_pixelcount
 # d$lucfip_pixelcount <- d$lucfip_pixelcount_30th# pour l'instant on met lucfip_pixelcount_total dans "all producers" et pas rapid + slow, car on n'a 
 # # pas calculé rapid et slow pour ce type de forêt encore
 # 
 # # make variable that counts lucfp events on both small and medium sized plantations 
 # d$lucpfsmp_pixelcount <- d$lucpfsp_pixelcount_total + d$lucpfmp_pixelcount_total
 # d$lucfsmp_pixelcount <- d$lucfsp_pixelcount_30th + d$lucfmp_pixelcount_30th
 # 
 # ## make a variable that counts lucfp events on all types of plantations
 # d$lucpfap_pixelcount <- d$lucpfip_pixelcount + d$lucpfsmp_pixelcount
 # d$lucfap_pixelcount <- d$lucfip_pixelcount + d$lucfsmp_pixelcount 
 # 
 # 
 # d <- dplyr::select(d, -lucpfap_pixelcount_total, -lucfap_pixelcount_30th, -lucpfip_rapidslow_pixelcount, -lucpfsmp_pixelcount_total, -lucfsmp_pixelcount_30th)
 # 
 
 uni_lonlat <- unique(d$lonlat)
 d <- mutate(d, 
             parcel_id = match(lonlat, uni_lonlat)) 
  rm(uni_lonlat)
  
 saveRDS(d, file.path(paste0("temp_data/panel_parcels_ip_final_",
                             parcel_size/1000,"km_",
                             CR/1000,"CR.rds")))
 }

parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                    parcel_size/1000,"km_50CR.rds")))

variables <- c("wa_prex_cpo_imp1","wa_prex_cpo_imp2",       
               "wa_pct_own_cent_gov_imp", "wa_pct_own_loc_gov_imp",  "wa_pct_own_nat_priv_imp", "wa_pct_own_for_imp",     
               #"wa_concentration_10",     "wa_concentration_30", "wa_concentration_50",     
               "n_reachable_uml")#"n_reachable_ibs", "n_reachable_ibsuml",   "sample_coverage"

for(voi in variables){
  ## lags
  parcels <- dplyr::arrange(parcels, lonlat, year)
  parcels <- DataCombine::slide(parcels,
                                Var = voi, 
                                TimeVar = "year",
                                GroupVar = "lonlat",
                                NewVar = paste0(voi,"_lag1"),
                                slideBy = -1, 
                                keepInvalid = TRUE)
  parcels <- dplyr::arrange(parcels, lonlat, year)
}



saveRDS(parcels, file.path(paste0("temp_data/panel_parcels_ip_final_",
                               parcel_size/1000,"km_50CR.rds")))




  CR <- 3e4
  
  parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                      parcel_size/1000,"km_",CR/1000,"CR.rds")))
  names(parcels)
  
  # ONCE  wa_at_parcels_distances.R IS RERUN ONCE, THE PREPARATION CODE WILL RATHER BE 
  # nothing actually. 
  
  # make spatial
  parcels_cs <- parcels[!duplicated(parcels$lonlat),c("lonlat", "year", "idncrs_lat", "idncrs_lon", "lon", "lat")]
  parcels_cs <- st_as_sf(parcels_cs, coords = c("idncrs_lon", "idncrs_lat"), remove = FALSE, crs = indonesian_crs)
  row.names(parcels_cs) <- parcels_cs$lonlat
  
  # define the set of mills within a distance. 
  sgbp <- st_is_within_distance(parcels_cs, ibs_cs, dist = catchment_radius, sparse = TRUE)
  
  usgbp <- unique(sgbp)
  
  # takes ~20 minutes
  parcels_cs$reachable <- rep(NA,nrow(parcels_cs))
  for(i in 1:length(usgbp)){
    parcels_cs$reachable[sgbp %in% usgbp[i]] <- i
  }
  
  
  # Define the nearest mill 
  nearest_mill_idx <- st_nearest_feature(parcels_cs, ibs_cs)
  
  parcels_cs$nearest_firm_id <- ibs_cs$firm_id[nearest_mill_idx]
  
  parcels_cs <- st_drop_geometry(parcels_cs)
  
  parcels <- left_join(parcels, parcels_cs[,c("lonlat", "reachable", "nearest_firm_id")], by = "lonlat")

  saveRDS(parcels, file.path(paste0("temp_data/panel_parcels_ip_final_",
                                    parcel_size/1000,"km_",CR/1000,"CR.rds")))






