d_clean %>% names()
head(d_clean)
d_clean_sf <- st_as_sf(d_clean, coords = c("lon","lat"), crs = 4326)
d_clean_sf <- st_transform(d_clean_sf, crs = indonesian_crs)

d_clean_sp <- as(d_clean_sf, "Spatial")

library(rgdal)

library(spdep)

moran_t <- list()
for(t in sort(unique(d_clean$year))){
  
  d_clean_cs <- d_clean_sp[d_clean_sp$year == t,]
  
  #work out the neighbours
  d_clean_knn <- knearneigh(d_clean_cs, k = 8)
  
  d_clean_nb <- knn2nb(d_clean_knn)
  
  #reformat the list of neighbours to a weight matrix
  d_clean_lw <- nb2listw(d_clean_nb)
  
  #run the test
  moran_t[[match(t,sort(unique(d_clean$year)))]] <- moran.test(d_clean_cs@data$lucpfap_pixelcount, d_clean_lw, alternative = "greater")

}

