island <- "Kalimantan"
travel_time <- 2
catchment_radius <- 10000
parcel_size <- 3000

Island <- island
CR <- catchment_radius
sample <- "IBS"
TT <- travel_time
dyna <- "rapid"
pf_type <- "total"
forest <- pf_type
size <- "i"

if(forest %in% c("intact", "degraded", "total")){
  lucf <- "lucpf"
}
if(forest %in% c("30th", "60th", "90th")){
  lucf <- "lucf"
}

#### WITHOUT DYNAMICS #### 

#### CATCHMENT RADIUS STREAM ####

# just total forest
parcels_brick_name <- paste0("parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type)
parcels_brick_cr <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif"))) # 236174 grid cells

m.df_wide_cr <- raster::as.data.frame(parcels_brick_cr, na.rm = TRUE, xy = TRUE, centroids = TRUE)
nrow(m.df_wide_cr) # 86259

# the file that is save at the end of the selection on CA
m.df_cr <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",island,"_",parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR_", "total.rds")))
nrow(m.df_cr) # 173412

# the code at the end of prepare_lucpfip_dyn, that cbind together intact, degraded and total
df_intact   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_intact.rds")))
df_degraded <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_degraded.rds")))
df_total    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_total.rds")))

df_degraded <- dplyr::select(df_degraded, -lon, -lat)
df <- inner_join(df_intact, df_degraded, by = c("parcel_id", "year"))

nrow(df_intact) == nrow(df_degraded)# 173412
nrow(df_degraded) == nrow(df)# 173412






#### CATCHMENT AREA STREAM ####
# the parcels_brick (masked) from the **CA** stream of dynamics (so in make_osrm_CA)

parcels_brick_name <- paste0("parcel_",lucf,size,"p_",island,"_",parcel_size/1000,"km_",forest)
parcels_brick_ca <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif"))) # 236174

ibs_msk_df_ca <- raster::as.data.frame(parcels_brick_ca, na.rm = TRUE, xy = TRUE, centroids = TRUE)
nrow(ibs_msk_df_ca) # 86259

# the file that is save at the end of the selection on CA
m.df_ca <- readRDS(file.path(paste0("temp_data/processed_parcels/",lucf,size,"p_panel_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_CA_", forest,".rds")))
nrow(m.df_ca) # 1181088 THIS IS THE ONE THAT IS DIFFERENT FROM WITH DYNAMICS CA

df_indus <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_total.rds")))
nrow(df_indus) #1181088
# so cross section of parcels that one year at least are within 2 hours from a mill is 65616 grid cells. 

# that's 11106 records more than with dynamics (CA), so 617 more cells by cross section. 
# although m.df_wide data sets have the same number of parcels. 








#### WITH DYNAMICS #### 
#### CATCHMENT RADIUS STREAM ####
# the parcels_brick (masked) from the **CR** stream of dynamics (so in 3. of prepare_lucpfip_dyn)

# RAPID
parcels_brick_name <- paste0("parcel_lucpfip_rapid_",island,"_",parcel_size/1000,"km_total")
parcels_brick_cr_rapid <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif"))) # 236174 grid cells

m.df_wide_cr_rapid <- raster::as.data.frame(parcels_brick_cr_rapid, na.rm = TRUE, xy = TRUE, centroids = TRUE)
nrow(m.df_wide_cr_rapid) # 86259

# the file that is save at the end of the selection on CA
m.df_cr_rapid <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_rapid_",island,"_",parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR_", "total.rds")))
nrow(m.df_cr_rapid) # 173412



# REPLACE
parcels_brick_name <- paste0("parcel_lucpfip_replace_",island,"_",parcel_size/1000,"km_total")
parcels_brick_cr_replace <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif"))) # 236174 grid cells

m.df_wide_cr_replace <- raster::as.data.frame(parcels_brick_cr_replace, na.rm = TRUE, xy = TRUE, centroids = TRUE)
nrow(m.df_wide_cr_replace) # 86259 

# the file that is save at the end of the selection on CA
m.df_cr_replace <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_replace_",island,"_",parcel_size/1000,"km_",catchment_radius/1000,"km_IBS_CR_", "total.rds")))
nrow(m.df_cr_replace) # 173412


# the code at the end of prepare_lucpfip_dyn, that cbind together replace, rapid and slow
df_replace   <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_replace_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_total.rds")))
df_rapid <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_rapid_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_total.rds")))
df_slow    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_slow_",Island,"_",PS/1000,"km_",CR/1000,"km_",sample,"_CR_total.rds")))

df_rapid <- dplyr::select(df_rapid, -lon, -lat)
df <- inner_join(df_replace, df_rapid, by = c("parcel_id", "year"))

nrow(df_replace) == nrow(df_rapid) == nrow(df)# 173412
nrow(df_rapid) == nrow(df)# 173412




#### CATCHMENT AREA STREAM ####
# the parcels_brick (masked) from the **CA** stream of dynamics (so in make_osrm_CA)

# RAPID
parcels_brick_name <- paste0("parcel_lucpfip_rapid_",island,"_",parcel_size/1000,"km_total")
parcels_brick_ca_rapid <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif"))) # 236174

ibs_msk_df_ca_rapid <- raster::as.data.frame(parcels_brick_ca_rapid, na.rm = TRUE, xy = TRUE, centroids = TRUE)
nrow(ibs_msk_df_ca_rapid) # 86259

# the file that is save at the end of the selection on CA
m.df_ca_rapid <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_rapid_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_CA_", "total.rds")))
nrow(m.df_ca_rapid) # 1169982



# REPLACE
parcels_brick_name <- paste0("parcel_lucpfip_replace_",island,"_",parcel_size/1000,"km_total")
parcels_brick_ca_replace <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif"))) # 236174

ibs_msk_df_ca_replace <- raster::as.data.frame(parcels_brick_ca_replace, na.rm = TRUE, xy = TRUE, centroids = TRUE)
nrow(ibs_msk_df_ca_replace) # 86259

# the file that is save at the end of the selection on CA
m.df_ca_replace <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_replace_",island,"_",parcel_size/1000,"km_",travel_time,"h_IBS_CA_", "total.rds")))
nrow(m.df_ca_replace) # 1169982


# the code at the end of prepare_lucpfip_dyn, that cbind together replace, rapid and slow
df_replace    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_replace_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_total.rds")))
df_rapid    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_rapid_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_total.rds")))
df_slow    <- readRDS(file.path(paste0("temp_data/processed_parcels/lucpfip_panel_slow_",Island,"_",PS/1000,"km_",TT,"h_",sample,"_CA_total.rds")))

df_replace <- dplyr::select(df_replace, -idncrs_lon, -idncrs_lat, -lon, -lat)
df <- inner_join(df_slow, df_replace, by = c("parcel_id", "year"))
nrow(df_replace) == nrow(df_rapid) # 1169982
nrow(df_rapid) == nrow(df)# 1169982


1169982/18











## comparing the two blocks of code within make_osrm_CA, without and with dynamicz
dur_mat <- readRDS(file.path(paste0("input_data/local_osrm_outputs/osrm_driving_durations_",island,"_",parcel_size/1000,"km_IBS")))
dur_mat <- dur_mat$durations

#### Without dynamics, whole code ####
  if(forest %in% c("intact", "degraded", "total")){
    lucf <- "lucpf"
  }
  if(forest %in% c("30th", "60th", "90th")){
    lucf <- "lucf"
  }
  
  
  # Now mask rasters for all outcomes
  parcels_brick_name <- paste0("parcel_",lucf,size,"p_",island,"_",parcel_size/1000,"km_",forest)
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
  
  mask(x = parcels_brick, mask = total_ibs_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(parcels_brick)
  
  # Turn the masked raster to a sf data frame
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))
  
  # note the na.rm = TRUE 
  nodyn_ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  # the whole point of the below is to merge more safely the duration matrix and the parcel data frame. 
  nodyn_ibs_msk_df <- nodyn_ibs_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
  nodyn_ibs_msk_df <- st_as_sf(nodyn_ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  nodyn_ibs_msk_df <- st_transform(nodyn_ibs_msk_df, crs = 4326)
  nodyn_ibs_msk_df$lon <- st_coordinates(nodyn_ibs_msk_df)[,"X"]%>% round(6)
  nodyn_ibs_msk_df$lat <- st_coordinates(nodyn_ibs_msk_df)[,"Y"]%>% round(6)
  nodyn_ibs_msk_df <- mutate(nodyn_ibs_msk_df, lonlat = paste0(lon, lat))
#  nodyn_ibs_msk_df <- mutate(nodyn_ibs_msk_df, lonlat = paste0(round(lon,6), round(lat,6)))
  length(unique(nodyn_ibs_msk_df$lonlat)) == nrow(nodyn_ibs_msk_df)
  #nodyn_ibs_msk_df <- st_drop_geometry(nodyn_ibs_msk_df)
  
  # Besides, make IDs 
  island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
  nodyn_ibs_msk_df$parcel_id <- paste0(island_id, c(1:nrow(nodyn_ibs_msk_df))) %>% as.numeric()
  
  
  
  # keep only the parcels that satisfy the travel time condition for at least one mill in the whole period.
  dur_mat_log <- dur_mat/(60) < travel_time
  dur_mat_log <- as.data.frame(dur_mat_log)
  
  # The row.names have been computed the same way the nodyn_ibs_msk_df$lonlat column was is computed above. 
  colnames(dur_mat_log) <- paste0("firm_id_",colnames(dur_mat_log))
  dur_mat_log$lonlat <- row.names(dur_mat_log)
  
  # THIS IS THE LINES THAT SELECT PARCELS ON THEIR BEING WITHIN THE CATCHMENT AREA
  nodyn_ibs_msk_TT_df <- merge(nodyn_ibs_msk_df, dur_mat_log, by = "lonlat", all = FALSE) # 86259
  nodyn_ibs_msk_TT_df <- nodyn_ibs_msk_TT_df[base::rowSums(nodyn_ibs_msk_TT_df[,grepl("firm_id",colnames(nodyn_ibs_msk_TT_df))], na.rm = TRUE)>0,]
  # 65616




#### With dynamics, whole code ####
  parcels_brick_name <- paste0("parcel_lucpfip_",dyna,"_",island,"_",parcel_size/1000,"km_total")
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, ".tif")))
  
  mask(x = parcels_brick, mask = total_ibs_ca_sp,
       filename = file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")),
       datatype = "INT4U",
       overwrite = TRUE)
  
  rm(parcels_brick)
  
  # Turn the masked raster to a sf data frame
  parcels_brick <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name, "_IBS_masked.tif")))
  
  # note the na.rm = TRUE 
  dyn_ibs_msk_df <- raster::as.data.frame(parcels_brick, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  # the whole point of the below is to merge more safely the duration matrix and the parcel data frame. 
  dyn_ibs_msk_df <- dyn_ibs_msk_df %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
  dyn_ibs_msk_df <- st_as_sf(dyn_ibs_msk_df, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  dyn_ibs_msk_df <- st_transform(dyn_ibs_msk_df, crs = 4326)
  dyn_ibs_msk_df$lon <- st_coordinates(dyn_ibs_msk_df)[,"X"]%>% round(6)
  dyn_ibs_msk_df$lat <- st_coordinates(dyn_ibs_msk_df)[,"Y"]%>% round(6)
  dyn_ibs_msk_df <- mutate(dyn_ibs_msk_df, lonlat = paste0(lon, lat))
  #dyn_ibs_msk_df <- mutate(dyn_ibs_msk_df, lonlat = paste0(round(lon,6), round(lat,6))) 
  length(unique(nodyn_ibs_msk_df$lonlat)) == nrow(dyn_ibs_msk_df)
  #dyn_ibs_msk_df <- st_drop_geometry(dyn_ibs_msk_df)
  
  # Besides, make IDs 
  island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
  dyn_ibs_msk_df$parcel_id <- paste0(island_id, c(1:nrow(dyn_ibs_msk_df))) %>% as.numeric()
  
  
  # keep only the parcels that satisfy the travel time condition for at least one mill in the whole period.
  dur_mat_log <- dur_mat/(60) < travel_time
  dur_mat_log <- as.data.frame(dur_mat_log)
  
  # The row.names have been computed the same way the dyn_ibs_msk_df$lonlat column was is computed above. 
  colnames(dur_mat_log) <- paste0("firm_id_",colnames(dur_mat_log))
  dur_mat_log$lonlat <- row.names(dur_mat_log)
  
  dyn_ibs_msk_TT_df <- merge(dyn_ibs_msk_df, dur_mat_log, by = "lonlat", all = FALSE) # 85221
  dyn_ibs_msk_TT_df <- dyn_ibs_msk_TT_df[base::rowSums(dyn_ibs_msk_TT_df[,grepl("firm_id",colnames(dyn_ibs_msk_TT_df))], na.rm = TRUE)>0,]
  # 64999
  
  # find those 1038 parcels that do not share the lonlat
  # those in nodyn that are not in dyn
  nodyn_ibs_msk_df$lonlat[!(nodyn_ibs_msk_df$lonlat %in% dyn_ibs_msk_df$lonlat)] %>% head()
 
  # those in dyn that are not in nodyn
  dyn_ibs_msk_df$lonlat[!(dyn_ibs_msk_df$lonlat %in% nodyn_ibs_msk_df$lonlat)] %>% head()
  
  
  
  
  plot(st_geometry(dyn_ibs_msk_df[!(dyn_ibs_msk_df$lonlat %in% nodyn_ibs_msk_df$lonlat),]), add = T) 
  plot(st_geometry(island_sf))
  plot(st_geometry(nodyn_ibs_msk_df))
  

  all.equal(nodyn_ibs_msk_df$parcel_id, dyn_ibs_msk_df$parcel_id)
  all.equal(nodyn_ibs_msk_df$lonlat, dyn_ibs_msk_df$lonlat)
  length(unique(nodyn_ibs_msk_df$lonlat))
  nodyn_ibs_msk_df$lonlat[!(dyn_ibs_msk_df$lonlat %in% nodyn_ibs_msk_df$lonlat)]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# Now between replace and the two others, why are the parcel_id not the same (and why in Kalimantan only?)
  
  ### REPLACE
  
  # Now mask rasters for all outcomes
  parcels_brick_name_replace <- paste0("parcel_lucpfip_replace_",island,"_",parcel_size/1000,"km_total")

  
  # Turn the masked raster to a sf data frame
  parcels_brick_replace <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name_replace, "_IBS_masked.tif")))
  
  # note the na.rm = TRUE 
  ibs_msk_df_replace <- raster::as.data.frame(parcels_brick_replace, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  # the whole point of the below is to merge more safely the duration matrix and the parcel data frame. 
  ibs_msk_df_replace <- ibs_msk_df_replace %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
  ibs_msk_df_replace <- st_as_sf(ibs_msk_df_replace, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  ibs_msk_df_replace <- st_transform(ibs_msk_df_replace, crs = 4326)
  ibs_msk_df_replace$lon <- st_coordinates(ibs_msk_df_replace)[,"X"]%>% round(6) # the rounding is bc otherwise there are very little differences in the decimals of the coordinates... 
  ibs_msk_df_replace$lat <- st_coordinates(ibs_msk_df_replace)[,"Y"]%>% round(6) 
  ibs_msk_df_replace <- mutate(ibs_msk_df_replace, lonlat = paste0(lon, lat))
  ibs_msk_df_replace <- st_drop_geometry(ibs_msk_df_replace)
  
  # Besides, make IDs 
  island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
  ibs_msk_df_replace$parcel_id <- paste0(island_id, c(1:nrow(ibs_msk_df_replace))) %>% as.numeric()
  
  
  ### RAPID
  # Now mask rasters for all outcomes
  parcels_brick_name_rapid <- paste0("parcel_lucpfip_rapid_",island,"_",parcel_size/1000,"km_total")
  
  
  # Turn the masked raster to a sf data frame
  parcels_brick_rapid <- brick(file.path("temp_data/processed_lu", paste0(parcels_brick_name_rapid, "_IBS_masked.tif")))
  
  # note the na.rm = TRUE 
  ibs_msk_df_rapid <- raster::as.data.frame(parcels_brick_rapid, na.rm = TRUE, xy = TRUE, centroids = TRUE)
  
  # the whole point of the below is to merge more safely the duration matrix and the parcel data frame. 
  ibs_msk_df_rapid <- ibs_msk_df_rapid %>% dplyr::rename(idncrs_lon = x, idncrs_lat = y)
  ibs_msk_df_rapid <- st_as_sf(ibs_msk_df_rapid, coords = c("idncrs_lon", "idncrs_lat"), crs = indonesian_crs, remove = FALSE)
  ibs_msk_df_rapid <- st_transform(ibs_msk_df_rapid, crs = 4326)
  ibs_msk_df_rapid$lon <- st_coordinates(ibs_msk_df_rapid)[,"X"]%>% round(6) # the rounding is bc otherwise there are very little differences in the decimals of the coordinates... 
  ibs_msk_df_rapid$lat <- st_coordinates(ibs_msk_df_rapid)[,"Y"]%>% round(6) 
  ibs_msk_df_rapid <- mutate(ibs_msk_df_rapid, lonlat = paste0(lon, lat))
  ibs_msk_df_rapid <- st_drop_geometry(ibs_msk_df_rapid)
  
  # Besides, make IDs 
  island_id <- if(island == "Sumatra"){1} else if(island == "Kalimantan"){2} else if (island == "Papua"){3}
  ibs_msk_df_rapid$parcel_id <- paste0(island_id, c(1:nrow(ibs_msk_df_rapid))) %>% as.numeric()
  
  
  
  
  ibs_msk_df_replace$parcel_id %>% unique() %>%  length()
  ibs_msk_df_rapid$parcel_id %>% unique() %>%  length()
  
  ibs_msk_df_replace$lonlat %>% unique() %>%  length()
  ibs_msk_df_rapid$lonlat %>% unique() %>%  length()
  
  
  
  
  
  
  

