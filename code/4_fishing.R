## File: fishing.R
## Purpose: compute fishing effort
## Author: Filippomaria Cassarino
## Date: 

## Notes ----


# fish raster resolution is too high to be true
#
## Library ---- 

library(dplyr)
library(tidyr) # replace_na function
library(ggplot2)
library(GGally) # collinearity plot
library(terra)
library(sf)
library(lubridate) # to work with dates
library(rnaturalearth) # land shape, better and faster than mapdata
library(modi) # to compute cwv
library(janitor) # clean_names 
library(usethis) # edit Rprofile

#
## Import Global Fishing Watch data directly (Requires external actions) ----

# Global Fishing Watch data can be imported as data.frame from raster files
# directly into RStudio. The function to retrieve the data is in the gfwr package
install.packages("gfwr",
                 repos = c("https://globalfishingwatch.r-universe.dev",
                           "https://cran.r-project.org")) 

# To use the function, an API token is needed, which can be requested.
# The API token must then be set as a new variable in the Renviron or Rprofile

# Edit the r profile, setting the token as a variable 
usethis::edit_r_profile() 

# Check that the token is there - it should return a very long code
Sys.getenv("GFW_TOKEN") 

# Create the object key with the token
key <- gfw_auth()

# Data importation: recover data by gear from 2012 (start of gfw data series) 
# to 2024 

# Create an sf object to define the spatial extent of the data
polygon <- sf::st_bbox(c(xmin = 6, xmax = 48, ymin = 69, ymax = 83),
                       crs = 4326) |>
  sf::st_as_sfc() |>
  sf::st_as_sf()

# Create the start and end date vectors for a loop using the required format
start_date <- paste0(2012:2024, "-01-01") # 13 years
end_date <- paste0(2012:2024, "-12-31")

# Effort monthly values at 0.01°
# I get the first year (2021) outside the loop to then rbind the rest to it
effort <- get_raster(
  spatial_resolution = "HIGH", #HIGH = 0.01°; LOW = 0.1°
  temporal_resolution = "MONTHLY", #DAILY, MONTHLY, YEARLY
  start_date = start_date[1], # this must be 366 days or less`
  end_date = end_date[1],
  region = polygon, #sf shapefile or GFW code
  region_source = "USER_SHAPEFILE", #EEZ, RFMO, MPA, or USER_SHAPEFILE
  group_by = "GEARTYPE" #VESSEL_ID, FLAG, GEARTYPE, or FLAGANDGEARTYPE
)

# Now a loop for the other years
for(i in 2:13){
  more_effort <- get_raster(
    spatial_resolution = "HIGH", #HIGH = 0.01°; LOW = 0.1°
    temporal_resolution = "MONTHLY", #DAILY, MONTHLY, YEARLY
    start_date = start_date[i], # this must be 366 days or less`
    end_date = end_date[i],
    region = polygon, #sf shapefile or GFW code
    region_source = "USER_SHAPEFILE", #EEZ, RFMO, MPA, or USER_SHAPEFILE
    group_by = "GEARTYPE" #VESSEL_ID, FLAG, GEARTYPE, or FLAGANDGEARTYPE
  )
  
  effort <- rbind(effort, more_effort)
}

rm(more_effort, polygon, end_date, start_date, key, i) 

# Save
write.csv(effort, "data/original/effort_gfw_0.01_monthly.csv", row.names=F)



## Compute total yearly fishing effort ----

# Load data
effort <- read.csv("data/original/effort_gfw_0.01_monthly.csv")

# Filter to keep only what is need now and summarize
trawling <- effort %>%
  filter(geartype == "trawlers") %>%
  mutate(year = year(ym(Time.Range))) %>%
  filter(year < 2023) %>%  
  group_by(year, Lat, Lon) %>%
  summarize(fishing_hours = sum(Apparent.Fishing.Hours)) %>%
  as.data.frame()

rm(effort) 

# Save for plotting
save(trawling, file = "data/intermediate/trawling.RData")

## Extract raster values ----

# Load data
load("data/intermediate/trawling.RData")
load("data/intermediate/haul_point.RData")

# Objects for the extraction
raster_list <- list()       # empty list to store rasters
ext <- ext(6, 48, 69, 83)   # extent of the output rasters
y <- unique(trawling$year)  # index for loops

# Loop to have a rater layer per year
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

# Stack all layers into a single SpatRaster
fish_stack <- rast(raster_list)

# Assign proper layer names (years)
names(fish_stack) <- unique(trawling$year)

# Project (new resolution is 414.4856 x 414.4856 m)        
equal <- project(fish_stack, 
                 "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs",
                 res = 2) 

# Filter the vector and raster by year
fishing_effort_terra <- dplyr::filter(haul_point, year == y[1])  # this will be the product
raster <- equal[[1]] 

# Extract the raster values in the points and add them to the vector
var <- terra::extract(x = raster, y = vect(fishing_effort_terra), # extract raster values
                      method = "bilinear") # interpolates from closest 4 cells

# Add the values as a  new column to the vector
fishing_effort_terra$fishing_effort <- var[, -1] 

# Loop to paste all subsequent years
for (i in 2:length(y)){ # the first was done outside the loop
  point <- dplyr::filter(haul_point, year == y[i]) # this will be the final df
  raster <- equal[[i]] 
  var <- terra::extract(x = raster, y = vect(point), # extract raster values
                        method = "bilinear") 
  point$fishing_effort <- var[, -1]
  fishing_effort_terra <- rbind(fishing_effort_terra, point)
}

fishing_effort <- fishing_effort_terra %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  select(haul_id, fishing_effort)

# Save
save(fishing_effort, file = "data/intermediate/fishing_effort.RData")

## End