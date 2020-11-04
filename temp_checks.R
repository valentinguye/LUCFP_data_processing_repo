


island <- "Sumatra"


###   2000FC EXTENT
in_fc_type_file_name <- paste0("parcel_fc2000_",island,"_30th_",parcel_size/1000,"km_IBS_masked")
m.parcels <- raster(file.path(paste0("temp_data/processed_lu/",in_fc_type_file_name,".tif")))
m.df <- raster::as.data.frame(m.parcels, na.rm = TRUE, xy = TRUE, centroids = TRUE)
m.df <- m.df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
m.df <- st_as_sf(m.df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
m.df <- st_transform(m.df, crs = 4326)
m.df$lon <- st_coordinates(m.df)[,"X"]
m.df$lat <- st_coordinates(m.df)[,"Y"]
m.df <- st_drop_geometry(m.df)



### LUCPFIP
parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_total")
parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
mask(x = parcels_brick, mask = total_ibs_ca_sp,
     filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")),
     datatype = "INT4U",
     overwrite = TRUE)

rm(parcels_brick)
# Turn the masked raster to a sf data frame
parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))
# note the na.rm = TRUE. Due to the mask, this keeps in the df object only the parcel within 80km of a mill at least one year. 
ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
ibs_msk_df <- ibs_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
ibs_msk_df <- st_as_sf(ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
ibs_msk_df <- st_transform(ibs_msk_df, crs = 4326)
ibs_msk_df$lon <- st_coordinates(ibs_msk_df)[,"X"]
ibs_msk_df$lat <- st_coordinates(ibs_msk_df)[,"Y"]
ibs_msk_df <- st_drop_geometry(ibs_msk_df)
# do not need to make it spatial here. But make IDs
island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
ibs_msk_df$parcel_id <- paste0(island_id, c(1:nrow(ibs_msk_df))) %>% as.numeric()
row.names(ibs_msk_df) <- ibs_msk_df$parcel_id

### DURATION MATRIX BASED ON LUCPFIP
dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_IBS")))
dur_mat_log <- (dur_mat$durations/60) < travel_time
dur_mat_log <- as.data.frame(dur_mat_log)
dur_mat_log <- cbind(dur_mat_log, dur_mat$sources)
dur_mat_log$parcel_id <- row.names(dur_mat_log)
  
nrow(m.df)
nrow(ibs_msk_df)
nrow(dur_mat_log)

head(m.df[,c("lon","lat")])
head(ibs_msk_df[,c("lon","lat")])
head(dur_mat_log[,])
head(dur_mat$sources)
head(dur_mat$durations[,1])
head(ibs_msk_df_lonlat[,c("lon","lat")])

head(m.df[,c("idncrs_lon","idncrs_lat")])
head(ibs_msk_df[,c("idncrs_lon","idncrs_lat")])


m.df <- merge(m.df, dur_mat_log, by = c("lon", "lat"), all = FALSE)

