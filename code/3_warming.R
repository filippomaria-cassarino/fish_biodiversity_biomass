## File: warming.R
## Purpose: compute warming metrics
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

# decide where to put haul position (haul_point)

#
## Library ---- 

library(dplyr)
library(tidyr) # replace_na function
library(ggplot2)
library(GGally) # collinearity plot
library(terra)
library(sf)
library(rnaturalearth) # land shape, better and faster than mapdata
library(modi) # to compute cwv
library(janitor) # clean_names 

# Haul positions
load("data/intermediate/community_filtered.RData")

haul_point <-  community_filtered %>%
  group_by(haul_id, latitude, longitude, year) %>%
  summarize(depth = mean(depth)) %>%  # this allows to have one raw per haul 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326) %>%            # transform to point vector (sf object)
  st_transform(crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs")   

save(haul_point, file = "data/intermediate/haul_point.RData")

#
## Load and filter Copernicus files ----

# sst, sbt and sic (2003 - Jun 2021)  
phy_stack_1 <- rast("data/original/cmems_mod_glo_phy_my_0.083deg_P1M-m_1738920541876.nc")

# sst, sbt and sic (Jul 2021 - 2024) 
phy_stack_2 <- rast("data/original/cmems_mod_glo_phy_myint_0.083deg_P1M-m_1741103917033.nc")

# Stack the stacks
phy_stack <- c(phy_stack_1, phy_stack_2) # using terra::c 

# Extract year from time stamp
phy_times <- time(phy_stack)
years <- format(phy_times, "%Y")

# Retain layers from Jan 2004 to Dec 2022
phy_stack <- phy_stack[[!years %in% c("2003", "2023", "2024")]] # 2003 was downloaded to use "12 months before sampling"

# Subset to compute each variable on its own
sst_stack <- phy_stack[[grep("thetao", names(phy_stack))]]  # sst
sbt_stack <- phy_stack[[grep("bottomT", names(phy_stack))]] # sbt
sic_stack <- phy_stack[[grep("siconc", names(phy_stack))]]  # sic

rm(phy_stack, phy_stack_1, phy_stack_2, phy_times, years) # cleaning

# Visual check for raster time continuity 
#plot(sic_stack[[210:221]])

# Chla (2003-2024) 
chla_stack <- rast("data/original/cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_1738579225409.nc")

# Extract year from time stamp
chla_times <- time(chla_stack)
years <- format(chla_times, "%Y")

# Retain layers from Jan 2004 to DEc 2022
chla_stack <- chla_stack[[!years %in% c("2003", "2023", "2024")]] # 2003 was downloaded to use "12 months before sampling"

# Vector for loops
stacks <- c("sst_stack", "sbt_stack", "sic_stack", "chla_stack")
vars <- c("sst", "sbt", "sic", "chla")

# Check the rasters
for (i in 1:length(stacks)){
  print(paste0("-------", stacks[i], "-------"))
  print(get(stacks[i]))
} 

#
## Compute yearly mean and sd values ----

# Vectors for the loops
stacks <- c("sst_stack", "sbt_stack", "sic_stack", "chla_stack")
vars <- c("sst", "sbt", "sic", "chla")

# Compute mean values by year 
for (i in 1:length(stacks)){
  t <- time(get(stacks[i]))  # extract the time dimension format
  y <- format(t, "%Y")    # convert dates to years (index for mean calculation)
  m <- tapp(get(stacks[i]), index = y, # compute the mean per cell per year
            fun = mean, na.rm = TRUE)
  unique_y <- as.numeric(gsub("X", "", names(m))) # years of the aggregation 
  time(m) <- unique_y # assign the new time values to the raster stack
  writeRaster(m, file = paste0(  
    "data/intermediate/", vars[i], "_mean.tif"), overwrite = TRUE)  # save raster
  rm(y, t, m, unique_y)
}

# Compute sd values by year 
for (i in 1:length(stacks)){
  t <- time(get(stacks[i]))  # extract the time dimension format
  y <- format(t, "%Y")    # convert dates to years (index for mean calculation)
  m <- tapp(get(stacks[i]), index = y, # compute the sd per cell per year
            fun = sd, na.rm = TRUE)
  unique_y <- as.numeric(gsub("X", "", names(m))) # years of the aggregation 
  time(m) <- unique_y # assign the new time values to the raster stack
  writeRaster(m, file = paste0(  
    "data/intermediate/", vars[i], "_sd.tif"), overwrite = TRUE)  # save raster
  rm(y, t, m, unique_y)
}

rm(sst_stack, sbt_stack, sic_stack, chla_stack, i, stacks, vars) 

## Extract raster values at haul level ----

# Load haul positions
load("data/intermediate/haul_point.RData")

# Output data frame
warming <- data.frame(haul_id = haul_point$haul_id,
                              geometry = haul_point$geometry,
                              year = haul_point$year,
                              depth = haul_point$depth)

# Vectors for the loops
vars <- c("sst", "sbt", "sic", "chla") # outer loop
vars_mean <- c("sst_mean", "sbt_mean", "sic_mean", "chla_mean") # outer loop
vars_sd <- c("sst_sd", "sbt_sd", "sic_sd", "chla_sd") # outer loop
y <- 2004:2022 # inner loop 

# Loop to extract mean values
for (k in 1:length(vars_mean)){
  # Import the raster 
  var_mean <- rast(paste0(   # import mean rasters
    "data/intermediate/", vars[k], "_mean.tif"))
  
  # Project as haul_point (new resolution is 3464.996 x 3464.996 m for phy and 1727.024 x 1727.024 m for chla)
  equal <- project(var_mean, 
                   "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs",
                   res = 2) 
  
  # Filter the vector and raster by year
  point_vars <- dplyr::filter(haul_point, year == y[1]) # this will be the product
  raster <- equal[[1]] 
  
  # Extract the raster values at the points 
  var <- terra::extract(x = raster, y = vect(point_vars), # extract raster values
                        method = "bilinear") # interpolates from closest 4 cells
  
  # Add the values as a  new column to the vector
  point_vars[[vars_mean[k]]] <- var[, -1] 
  
  # Loop to paste all subsequent years
  for (i in 2:length(y)){ # the first was done outside the loop
    point <- dplyr::filter(haul_point, year == y[i]) # this will be the final df
    raster <- equal[[i]] 
    var <- terra::extract(x = raster, y = vect(point), # extract raster values
                          method = "bilinear") # interpolates from closest 4 cells
    point[[vars_mean[k]]] <- var[, -1]
    point_vars <- rbind(point_vars, point)
  }
  
  # Remove unwanted columns
  point_vars <- select(point_vars, - c(depth, year))
  
  # Join to the final dataset
  warming <- left_join(warming, point_vars,
                               by = c("haul_id" = "haul_id",
                                      "geometry" = "geometry"))
  
  rm(raster, var, point, equal, var_mean, point_vars) # cleaning
}

# Loop to extract sd values
for (k in 1:length(vars_sd)){
  # Import the raster 
  var_sd <- rast(paste0(   # import mean rasters
    "data/intermediate/", vars[k], "_sd.tif"))
  
  # Project as haul_point 
  equal <- project(var_sd, 
                   "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs",
                   res = 2) 
  
  # Filter the vector and raster by year
  point_vars <- dplyr::filter(haul_point, year == y[1]) 
  raster <- equal[[1]] 
  
  # Extract the raster values at the points 
  var <- terra::extract(x = raster, y = vect(point_vars), # extract raster values
                        method = "bilinear") # interpolates from closest 4 cells
  
  # Add the values as a  new column to the vector
  point_vars[[vars_sd[k]]] <- var[, -1] 
  
  # Loop to paste all subsequent years
  for (i in 2:length(y)){ # the first was done outside the loop
    point <- dplyr::filter(haul_point, year == y[i]) # this will be the final df
    raster <- equal[[i]] 
    var <- terra::extract(x = raster, y = vect(point), # extract raster values
                          method = "bilinear") # interpolates from closest 4 cells
    point[[vars_sd[k]]] <- var[, -1]
    point_vars <- rbind(point_vars, point)
  }
  
  # Remove unwanted columns
  point_vars <- select(point_vars, - c(depth, year))
  
  # Join to the final dataset
  warming <- left_join(warming, point_vars,
                               by = c("haul_id" = "haul_id",
                                      "geometry" = "geometry"))
  
  rm(raster, var, point, equal, var_sd, point_vars) # cleaning
}

rm(vars, vars_mean, vars_sd, y, i, k)

# Save 
save(warming, file = "data/intermediate/warming.RData")

rm(list = ls())

#

### MONTHS ------------------
# extract months
chla_times <- time(chla_stack)
months <- format(chla_times, "%m")  # Extract months

# retain layers from april (04) to july (07)
chla_stack <- chla_stack[[months %in% c("04", "05", "06", "07")]] 
rm(months, chla_times)  # cleaning the environment