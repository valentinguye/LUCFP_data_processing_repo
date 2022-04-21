
# These are the packages needed in this particular script. 
neededPackages = c("tibble", "plyr", "dplyr", "data.table",
                   "foreign", "readstata13", "readxl",
                   "raster", "rgdal",  "sp", "spdep", "sf",
                   "DataCombine",
                   "knitr", "kableExtra",
                   "car",  "fixest", "sandwich", "boot", "multcomp", "urca",# 
                   "ggplot2","leaflet")
#install.packages("sf", source = TRUE)
# library(sf)
# 
# neededPackages = c("tidyverse","data.table", "readxl","foreign", "data.table", "readstata13", "here",
#                    "rgdal", "raster", "velox","sp", "lwgeom", "rnaturalearth", 
#                    "rlist", "parallel", "foreach", "iterators", "doParallel" )
# 

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
# troublePackages <- c("")
# # Attempt to load packages from user's default libraries. 
# lapply(troublePackages, library, lib.loc = default_libraries, character.only = TRUE)

# 3. If the troubling packages could not be loaded ("there is no package called ...) 
#   you should try to install them, preferably in their versions stated in the renv.lock file. 
#   see in particular https://rstudio.github.io/renv/articles/renv.html 



### NEW FOLDERS USED IN THIS SCRIPT 

### INDONESIAN CRS 
#   Following http://www.geo.hunter.cuny.edu/~jochen/gtech201/lectures/lec6concepts/map%20coordinate%20systems/how%20to%20choose%20a%20projection.htm
#   the Cylindrical Equal Area projection seems appropriate for Indonesia extending east-west along equator. 
#   According to https://spatialreference.org/ref/sr-org/8287/ the Proj4 is 
#   +proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs
#   which we center at Indonesian longitude with lat_ts = 0 and lon_0 = 115.0 
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

### LEGAL LAND USE 
llu <- st_read(file.path("input_data/kawasan_hutan/Greenorb_Blog/final/KH-INDON-Final.shp"))
# llu <- st_transform(llu, crs = indonesian_crs)

names(llu)[names(llu) == "Fungsi"] <- "llu"
unique(llu$llu)
llu <- llu[,"llu"]

length(unique(llu$geometry))
llu <- llu[!duplicated(llu$geometry),]

# According to MoF 2020: https://kemlu.go.id/download/L1NoYXJlZCUyMERvY3VtZW50cy9Eb2t1bWVuX0luZm9ybWFzaS9UaGUlMjBTdGF0ZSUyMG9mJTIwSW5kb25lc2lhJTIwRm9yZXN0JTIwMjAyMCUyMChSZXByaW50ZWQpLnBkZg==

# 3 main classes: HP, HL, HK

# Sub classes: 

# HP: production forest
#   - HP: Permanent Production Forest (Hutan Produksi Tetap)
#   - HPK: Convertible Production Forest (Hutan Produksiyang Dapat Dikonversi)
#   - HPT: Limited production forest (Hutan Produksi Terbatas)

# HK: conservation forest 
#   - KSA: Sanctuary Reserve Areas (Kawasan Suaka Alam)
#           - CA: Strict Nature Reserves (Cagar Alam)
#           - SM: Wildlife Sanctuaries (Suaka Margasatwa)
#   - KPA: Nature Conservation Areas (Kawasan Pelestarian Alam)
#           - TN: National Parks (Taman Nasional)
#           - TWA: Nature Recreation Parks (Taman Wisata Alam)
#           - Tahura: Grand Forest Parks (Taman Hutan Raya, Tahura)

# HL: protected forest
# no sub class

# Apparently stand-alone classes:
# TB: Hunting Park (Taman Buru)

# Some are unknown
# Among which some seem to have a "L" suffix - they would all belong to HK. 
# It seems they are the marine part of conservation forest ("marine KSA/KPA have islands with forests that are classified as marine KSA/KPA, because the majority of the area in these marine KSA/KPA are sea waters.", in MoF 2019)

# SML
# CAL 
# TNL 
# TWAL
# KSAL

# Hutan Cadangan: is not in MOEF but means "reserve forest". 
# Hutan Pangonan: ??? don't know what it means
# HSA: ???
# SA ?? 


HKL <- llu[( llu$llu=="SML" | 
               llu$llu=="CAL" | 
               llu$llu=="TNL" | 
               llu$llu=="TWAL" | 
               llu$llu=="KSAL" ), ]
plot(llu[llu$llu == "Hutan Cadangan",])

# this is in a rather small region (bbox: xmin: 106.4455 ymin: -7.786833 xmax: 108.7177 ymax: -6.373507) 
# and not many features (237)
plot(llu[llu$llu == "Hutan Pangonan", "llu"])
plot(countries, add = T)


llu <- mutate(llu, llu3 = if_else( (llu=="HK" | 
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
                                    # llu=="TNL" | # remove juste this one for visibility 
                                    llu=="TWAL" | 
                                    llu=="KSAL" | 
                                    llu=="Hutan Cadangan"), true = "HK", false = ""))
llu[llu$llu=="HP" | 
    llu$llu=="HPK" | 
    llu$llu=="HPT", "llu3"] <- "HP"

llu[llu$llu=="HL", "llu3"] <- "HL"

llu[llu$llu3=="",]$llu%>%unique()

unique(llu$llu3)

plot(llu[,"llu3"])

#dissolve
# llu_dissolve <- llu %>% group_by("llu3") %>% summarize() 

# this takes ~1 hour, especially HP
HK <- st_union(llu[llu$llu3=="HK",]) %>% st_sf() %>% mutate(llu3="HK")
HL <- st_union(llu[llu$llu3=="HL",]) %>% st_sf() %>% mutate(llu3="HL")
HP <- st_union(llu[llu$llu3=="HP",]) %>% st_sf() %>% mutate(llu3="HP")
#OTHER <- st_union(llu[llu$llu3=="",]) %>% st_sf() %>% mutate(llu3="")

llu_final <- rbind(HK, HL, HP) # ,OTHER

st_write(llu_final, "temp_data/legal_lu_4categories.shp")

plot(llu_final)


# prepare backgroud layers with other countries
countries <- st_read(file.path("input_data/Global_LSIB_Polygons_Detailed"))
countries <- countries[countries$COUNTRY_NA == "Indonesia" | 
                         countries$COUNTRY_NA == "Malaysia" | 
                         countries$COUNTRY_NA == "Brunei", "COUNTRY_NA"]
# these two lines to speed up mapping
countries <- st_transform(countries, crs = indonesian_crs) %>% st_simplify(dTolerance = 1000)
countries <- st_transform(countries, crs = 4326)

## Prepare llu_final for the plot

# First simplify
llu_final <-  st_transform(llu_final, crs = indonesian_crs) %>% st_simplify(dTolerance = 100)
llu_final <- st_transform(llu_final, crs = 4326)


# Change the values for the legend
llu_final[llu_final$llu3=="HL", "llu3"] <- "Protection forest"
llu_final[llu_final$llu3=="HK", "llu3"] <- "Conservation forest"
llu_final[llu_final$llu3=="HP", "llu3"] <- "Production forest"
# llu_final[llu_final$llu3=="", "llu3"] <- "No-status forest"


ggplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_brewer(type = "div", palette="Dark2") +
  #geom_sf(data = countries[countries$COUNTRY_NA!="Indonesia",], fill = NA, col = "black") +
  # geom_sf(data = countries[countries$COUNTRY_NA=="Indonesia",], col = alpha("grey",0.2)) +

  geom_sf(data = llu_final, aes(col = NA, fill = llu_final$llu3)) + 
  #geom_sf(data=st_geometry(d_geo), color=alpha("grey",0.2))+
  coord_sf(xlim = c(90, 145), ylim = c(-12, 7), expand = FALSE) # 145


plot(st_geometry(llu))












