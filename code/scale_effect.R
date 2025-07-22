
### scale effect ###

### LIBRARY ----
Sys.setenv(LANG = "en")
setwd("C:/Users/User/OneDrive - UGent/Desktop/Rstuff/master_paper")

library(dplyr)
library(ggplot2)
library(sf)            # spatial data manipulation (shp)
library(terra)         # spatial data manipulation (raster), successor of raster
library(ncdf4)         # netcdf manipulation
library(lubridate)     # working with date-time objects
library(rnaturalearth) # land shape
library(rlang)         # to dinamicly name colums in summarize

## haul positions
# import community data
load("C:/Users/User/OneDrive - UGent/Desktop/Rstuff/master_thesis/data/intermediate_data/community/NOR-BTS_filtered.RData") 

# retain only one raw per haul, transform into 
# a point vector and projects to equal area
nor_point <-  nor %>%
  group_by(haul_id, latitude, longitude, year) %>%
  summarize(depth = mean(depth)) %>%  # this allows to have one raw per haul 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326) %>%            # transform to point vector (sf object)
  st_transform(crs = "EPSG:3035")     # Lambert Azimuthal Equal-Area

rm(nor) # cleaning

## fishing data
load("C:/Users/User/OneDrive - UGent/Desktop/Rstuff/master_thesis/data/intermediate_data/fishing/trawling.RData")

# extraction from raster ----

## raster extraction
# objects for the extraction
raster_list <- list()       # empty list to store rasters
ext <- ext(6, 48, 69, 83)   # extent of the output rasters
y <- unique(trawling$year)  # index for loops

# loop to have a raster layer per year
for (i in y) {
  yr <- trawling %>%
    filter(year == i) %>%
    select(Lon, Lat, fishing_hours)  # reorder lat and long appropriately
  r <- rast(x = yr, type = "xyz", crs = "EPSG:4326")  # convert to raster
  r[is.na(r)] <- 0 # NA values are actually 0 fishing
  # the values need to be aligned and in the same extent. resample can do both
  r.1 <- resample(r, rast(ext = ext, crs = "EPSG:4326", res = 0.01),
                  method = "near") # "near" preserves original (which is a sum)
  raster_list[[as.character(i)]] <- r.1
}

# stack all layers into a single SpatRaster
fish_stack <- rast(raster_list) |> project("EPSG:3035")


# looped extraction

# scales
radius <- c(2e3, 4e3, 8e3, 16e3, 32e3)

extracted <- list()

# outer loop to extract values at different radia
for(k in seq_along(radius)){
  
# inner loop to extract values at each year
for(i in seq_along(fish_stack)) {
  
  # the circle vector as SpatVector 
  circle <- nor_point %>%
    filter(year == as.integer(names(fish_stack)[i])) %>%
    st_buffer(dist = radius[k]) 
  
  # the corresponding raster layer
  fishing <- fish_stack[[i]]
  
  # extraction
  e <- extract(fishing, vect(circle), fun = mean, na.rm = TRUE, ID = TRUE)
  
  e$haul_id <- circle$haul_id[e$ID]
  
  # placement in the list
  extracted[[length(extracted) + 1]] <- data.frame(haul_id = e$haul_id,
                                                   mean_hours = e$,
                                                   radius = radius[k])
  rm(e, fishing, circle)
}
}

fishing_scales <- bind_rows(extracted)
