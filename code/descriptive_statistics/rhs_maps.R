### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
#             EXPLICATIVE VARIABLE DESCRIPTIVE MAP FIGURES 
# 
#   input   - panel data frame of parcels with the information on number of reachable uml mills
#           as outputted from add_parcel_variables.R
#           ---> temp_data/processed_parcels/wa_panel_parcels_reachable_uml_PS_CR.rds
# 
#           - panel data frame of IBS 
#           ---> temp_data/IBS_UML_panel_final.dta
# 
#   output  none automatically, but figures generated with this code were saved. 
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
# PB: on n'arrive pas installed leaflet, kableExtra (qui n'est peut être pas nécessaire) et velox (qui n'est pas nécessaire)
neededPackages = c("plyr", "dplyr", "data.table",  "stringr", "sjmisc", 
                   "foreign", "readxl", "readstata13", 
                   "raster", "rgdal",  "sp", "sf",
                   "knitr", 
                   "parallel", "foreach","doParallel")
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
troublePackages <- c("leaflet", "leaflet.providers", "png")
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


# ISLAND SHAPES
island_sf <- st_read(file.path("temp_data/processed_indonesia_spatial/island_sf"))
names(island_sf)[names(island_sf)=="island"] <- "shape_des"

island_sf_prj <- st_transform(island_sf, crs = indonesian_crs)

# Accumulated LUCPFIP raster (template to rasterize RHS)

prepare_accu_lucpfip <- function(island){
  parcel_size <- 3000
  pf_type <- "total"
  brick_lucpfip <- brick(file.path(paste0("temp_data/processed_lu/parcel_lucpfip_",island,"_",parcel_size/1000,"km_",pf_type,".tif")))
  
  # remove layers of years aftr 2015
  brick_lucpfip <- raster::subset(brick_lucpfip, c(1:15))
  # Add up annual aggregated lucpfip
  accu_lucpfip <- calc(brick_lucpfip, fun = sum, na.rm = TRUE)
  # convert the sum value from number of pixels of resolution 27.7m (i.e; 767.29 square meters) to million hectare 
  pixel_area <- (27.8*27.6)/(1e4) # 1e4 is the convertion factor between a square meter and ONE hectare.  
  accu_lucpfip <- calc(accu_lucpfip, fun = function(x){x*pixel_area})
  # turn 0 (which are most parcels in the raster) to NA for transparence
  accu_lucpfip <- reclassify(accu_lucpfip, rcl = cbind(0,NA))
  
  return(accu_lucpfip)
}

accu_lucpfip_suma <- prepare_accu_lucpfip("Sumatra")
accu_lucpfip_kali <- prepare_accu_lucpfip("Kalimantan")
accu_lucpfip_papu <- prepare_accu_lucpfip("Papua")


parcel_size <- 3000
catchment_radius <- 30000


  # Prepare catchment radii polygons 
  ibs <- read.dta13(file.path("temp_data/IBS_UML_panel_final.dta"))
  # keep only a cross section of those that are geolocalized mills
  ibs <- ibs[ibs$analysis_sample==1,]
  ibs <- ibs[!duplicated(ibs$firm_id),]

  ibs <- st_as_sf(ibs,	coords	=	c("lon",	"lat"), crs=4326)
  ibs <- st_geometry(ibs)
  ibs_prj <- st_transform(ibs, crs = indonesian_crs)
  ibs_cr <- st_buffer(ibs_prj, dist = catchment_radius)
  ibs_cr <- st_union(ibs_cr)
  
  # we don't do that here because we want the catchment radius to encompass parcels within
  # and we did not remove parcels in the sea / on the coast. 
  # ibs_cr <- st_intersection(x = ibs_cr, y = island_sf_prj)
  
  # keep only the part of this total catchment area that is on our island of interest
  #ibs_cr <- st_intersection(x = ibs_cr, y = island_sf_prj[island_sf_prj$shape_des == island,])
  
  # un-project to leaflet lon lat 
  ibs_cr_lonlat <- st_transform(ibs_cr, crs = 4326)
  
  
  parcels <- readRDS(file.path(paste0("temp_data/panel_parcels_ip_final_",
                                      parcel_size/1000,"km_",
                                      catchment_radius/1000,"CR.rds")))  
  
  time_mean <- ddply(parcels, 
                     "parcel_id", 
                     summarise, 
                     tm_est_year_imp = mean(wa_est_year_imp, na.rm = TRUE),
                     tm_pffb_imp1 = mean(wa_ffb_price_imp1, na.rm = TRUE),
                     tm_pffb_imp2 = mean(wa_ffb_price_imp2, na.rm = TRUE),
                     tm_pcpo_imp1 = mean(wa_cpo_price_imp1, na.rm = TRUE),
                     tm_pcpo_imp2 = mean(wa_cpo_price_imp2, na.rm = TRUE),
                     tm_prex_cpo_imp1 = mean(wa_prex_cpo_imp1, na.rm = TRUE),
                     tm_prex_cpo_imp2 = mean(wa_prex_cpo_imp2, na.rm = TRUE),
                     tm_pct_own_for = mean(wa_pct_own_for_imp, na.rm = TRUE),
                     tm_pct_own_nat_priv = mean(wa_pct_own_nat_priv_imp, na.rm = TRUE),
                     tm_pct_own_cent_gov = mean(wa_pct_own_cent_gov_imp, na.rm = TRUE),
                     # tm_concentration_10 = mean(wa_concentration_10, na.rm = TRUE),
                     # tm_concentration_30 = mean(wa_concentration_30, na.rm = TRUE),
                     # tm_concentration_50 = mean(wa_concentration_50, na.rm = TRUE),
                     tm_n_reachable_ibs = mean(n_reachable_ibs, na.rm = TRUE),
                     tm_n_reachable_uml = mean(n_reachable_uml, na.rm = TRUE),
                     tm_n_reachable_uml = mean(n_reachable_uml, na.rm = TRUE),
                     tm_n_reachable_uml = mean(n_reachable_uml, na.rm = TRUE),
                     tm_sample_coverage = mean(sample_coverage, na.rm = TRUE),
                     tm_sample_coverage = mean(sample_coverage, na.rm = TRUE),
                     tm_sample_coverage = mean(sample_coverage, na.rm = TRUE)
                     )
  
  # round time mean variables for display purpose
  # all to 2 digits
  time_mean <- round(time_mean, digits = 2)
  # parcel_id and est_year to 0. 
  time_mean$parcel_id <- round(time_mean$parcel_id, digits = 0)
  time_mean$tm_est_year_imp <- round(time_mean$tm_est_year_imp, digits = 0)
  
  
  parcels <- merge(parcels, time_mean, by = "parcel_id")
  
  
  # keep one record per parcel
  parcels_cs <- parcels[!duplicated(parcels$parcel_id),]
  
  # those with NA all years produce a NaN despite na.rm = TRUE (it is the mean of nothing...)
  # and this makes rasterize not work. 
  for(nc in 1:ncol(parcels_cs)){
    parcels_cs[,nc][is.nan(parcels_cs[,nc])] <- -0.01 
    # this negative number enables us to identify these cases because we have only variables defined on positive segment.
  }
  
  # make spatial and split into islands
  parcels_cs <- st_as_sf(parcels_cs, coords = c("lon", "lat"), crs = indonesian_crs)
  
  parcels_suma <- parcels_cs[startsWith(as.character(parcels_cs$parcel_id), "1"),]
  parcels_kali <- parcels_cs[startsWith(as.character(parcels_cs$parcel_id), "2"),]
  parcels_papu <- parcels_cs[startsWith(as.character(parcels_cs$parcel_id), "3"),]
  
  # drop geometry of parcels_cs because need it as a mere df in the leaflet afterwards. 
  parcels_cs <- st_drop_geometry(parcels_cs)
  
  parcels_suma_sp <- as(parcels_suma, "Spatial")  
  parcels_kali_sp <- as(parcels_kali, "Spatial")  
  parcels_papu_sp <- as(parcels_papu, "Spatial")  
  
  template_suma <- projectExtent(accu_lucpfip_suma, crs = indonesian_crs) 
  template_kali <- projectExtent(accu_lucpfip_kali, crs = indonesian_crs) 
  template_papu <- projectExtent(accu_lucpfip_papu, crs = indonesian_crs) 
  
  raster_suma <- rasterize(parcels_suma_sp, template_suma, 
                      field = c("tm_est_year_imp", 
                                "tm_pffb_imp1", "tm_pffb_imp2",             
                                 "tm_pcpo_imp1", "tm_pcpo_imp2",
                                 "tm_prex_cpo_imp1", "tm_prex_cpo_imp2",
                                 "tm_pct_own_for", "tm_pct_own_nat_priv", "tm_pct_own_cent_gov",
                                 #"tm_concentration_10", "tm_concentration_30","tm_concentration_50",
                                 "tm_n_reachable_ibs", 
                                 "tm_n_reachable_uml",
                                 "tm_sample_coverage" ))
  
  raster_kali <- rasterize(parcels_kali_sp, template_kali, 
                      field = c("tm_est_year_imp", 
                                "tm_pffb_imp1", "tm_pffb_imp2",             
                                "tm_pcpo_imp1", "tm_pcpo_imp2",
                                "tm_prex_cpo_imp1", "tm_prex_cpo_imp2",
                                "tm_pct_own_for", "tm_pct_own_nat_priv", "tm_pct_own_cent_gov",
                                #"tm_concentration_10", "tm_concentration_30","tm_concentration_50",
                                "tm_n_reachable_ibs", 
                                "tm_n_reachable_uml",
                                "tm_sample_coverage"))
  
  raster_papu <- rasterize(parcels_papu_sp, template_papu, 
                      field = c("tm_est_year_imp", 
                                "tm_pffb_imp1", "tm_pffb_imp2",             
                                "tm_pcpo_imp1", "tm_pcpo_imp2",
                                "tm_prex_cpo_imp1", "tm_prex_cpo_imp2",
                                "tm_pct_own_for", "tm_pct_own_nat_priv", "tm_pct_own_cent_gov",
                                #"tm_concentration_10", "tm_concentration_30","tm_concentration_50",
                                "tm_n_reachable_ibs", 
                                "tm_n_reachable_uml",
                                "tm_sample_coverage"))

  # reclassify -0.01 to NAs (all these NAs before prevent rasterize to work)
  raster_suma <- reclassify(raster_suma, rcl = cbind(-0.01, NA))
  raster_kali <- reclassify(raster_kali, rcl = cbind(-0.01, NA))
  raster_papu <- reclassify(raster_papu, rcl = cbind(-0.01, NA))
  
  
  #### MAP FFB PRICE ####
  displayed_var <- "tm_pffb_imp1"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("BuPu", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                 domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                 #bins = 4, 
                 na.color = "transparent", 
                 reverse = FALSE)
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    #addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  
              values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
              bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = " USD/ton FFB"),
              #labels = bin_labels, # this lign does not do anything
              title = "Local FFB prices, <br/> 1998-2015 average",
              position = "bottomright") 

  
  #### MAP CPO PRICE ####
  displayed_var <- "tm_pcpo_imp1"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("YlOrRd", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    #addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var],
              bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = " USD/ton CPO"),
              #labels = bin_labels, # this lign does not do anything
              title = "Local CPO prices, <br/> 1998-2015 average",
              position = "bottomright") 
  
  
  
  
  #### MAP PCT OWN PUBLIC ####
  displayed_var <- "tm_pct_own_cent_gov"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = "%"),
              #labels = bin_labels, # this lign does not do anything
              title = "Central gvt. ownership, <br/> 1998-2015 average",
              position = "bottomright") 
    
 
  
  #### MAP PCT OWN PRIVATE ####
  displayed_var <- "tm_pct_own_for"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = "%"),
              #labels = bin_labels, # this lign does not do anything
              title = "Foreign ownership, <br/> 1998-2015 average",
              position = "bottomright") 
  
  
  #### MAP  CONCENTRATION ####
  displayed_var <- "tm_concentration_30"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("plasma", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP 
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = ""),
              #labels = bin_labels, # this lign does not do anything
              title = "Competing mills, <br/> 1998-2015 average",
              position = "bottomright") %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") 
  

  
  
  #### MAP  N REACHABLE UML ####
  displayed_var <- "tm_n_reachable_uml_30"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("plasma", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP 
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], bins = 5, opacity = 0.7,
              labFormat = labelFormat(suffix = ""),
              #labels = bin_labels, # this lign does not do anything
              title = "Reachable UML mills, <br/> 1998-2015 average",
              position = "bottomright")

  
  #### MAP  SAMPLE COVERAGE ####
  displayed_var <- "tm_sample_coverage_30km"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("viridis", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP 
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = "%"),
              #labels = bin_labels, # this lign does not do anything
              title = "Share of sample in <br/> all reachable mills <br/> 1998-2015 average",
              position = "bottomright") 

  
  
  
  
  #### MAP  EXPORT SHARE ####
  displayed_var <- "tm_prex_cpo_imp2"
  parcels_cs[parcels_cs[,displayed_var] == -0.01, displayed_var] <- NA
  
  # settings for lucfp legend 
  summary(parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var])
  # bins <- seq(from = 40, to = 180, by = 20)
  cb <- colorNumeric("BuPu", # "viridis" (green-purple), "magma" (yellow-purple), "inferno" (like magma), or "plasma", "BuPu", "Greens"
                     domain = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], 
                     #bins = 4, 
                     na.color = "transparent")
  # "viridis", "magma", "inferno", or "plasma".
  # cb <- colorNumeric(c("#0C2C84", "#41B6C4", "#FFFFCC"), values(accu_lucfp_suma),
  #                     na.color = "transparent")
  
  # this does not do anything: 
  # bin_labels <- c("0-300ha", "300-600ha", "600-900ha")
  # bin_labels <- paste0("<div style='display: inline-block;height: ", 
  #                      size, "px;margin-top: 4px;line-height: ", 
  #                      size, "px;'>", bin_labels, "</div>")
  
  # settings for catchment radius legend
  color <- "transparent"
  label <- paste0(catchment_radius/1000,"km catchment radius<br/>of IBS sample mills")
  shape <- "circle"
  border <- "black"
  size <- 10
  
  shape <- gsub("circle", "50%", shape)
  legend_colors <- paste0(color, "; width:", size, "px; height:", size, "px; border:3px solid ", border, "; border-radius:", shape)
  
  # MAP 
  ibs_cr_lonlat%>% 
    leaflet() %>% 
    addTiles()%>%
    addProviderTiles(providers$Esri.WorldImagery, group ="ESRI") %>%
    addPolygons(opacity = 0.5, color = border, weight = 2, fill = FALSE) %>%
    addRasterImage(raster::subset(raster_suma, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_kali, displayed_var), project = TRUE, colors = cb) %>% 
    addRasterImage(raster::subset(raster_papu, displayed_var), project = TRUE, colors = cb) %>% 
    addLegend(colors = legend_colors, labels = label, position = "bottomright") %>% 
    addLegend(pal = cb,  values = parcels_cs[parcels_cs$island=="Sumatra"|parcels_cs$island=="Kalimantan",displayed_var], bins = 4, opacity = 0.7,
              labFormat = labelFormat(suffix = "%"),
              #labels = bin_labels, # this lign does not do anything
              title = "CPO export share, <br/> 1998-2015 average",
              position = "bottomright") 

  
  
  
  
  


# Figures 


