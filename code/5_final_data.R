## File: final_data.R
## Purpose: joining the variables in a single dataset and prepare for analysis
## Author: Filippomaria Cassarino
## Date: 

## Notes ----
#
## Library ---- 

library(dplyr)
library(ggplot2)
library(tidyr) # pivot to wide format
library(terra) # ratser 
library(sf)    # vector
library(rnaturalearth) # land shape
library(modi) # compute cwv
library(janitor) # clean column names 

#
## Load datasets ----

# Biomass
load("data/intermediate/biomass.Rdata")

# Taxonomic diversity
load("data/intermediate/taxonomic_diversity.Rdata")

# Functional diversity
load("data/intermediate/functional_diversity.Rdata")

# Haul position, depth and warming variables
load("data/intermediate/warming.RData")

# Fishing effort
load("data/intermediate/fishing_effort.RData")

## Join data, correct NAs and duplicated coordinates ----

# Join, correct NAs and project to equidistant
plotting_data <- biomass %>% 
  left_join(functional_diversity, by = "haul_id") %>%
  left_join(taxonomic_diversity, by = "haul_id") %>%
  left_join(warming, by = "haul_id") %>%
  left_join(fishing_effort, by = "haul_id") %>%
  mutate(
    across(c(sic_mean, sic_sd), ~ replace_na(., 0)), # sic NAs are 0s
    fishing_effort = if_else(                        # fishing NAs from 2012 are 0s
      year >= 2012,
      replace_na(fishing_effort, 0),  
      fishing_effort))%>%
  mutate(log_fishing_effort = log(fishing_effort + 0.00001),
         log_chla_mean = log(chla_mean)) %>%
  st_as_sf(crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") %>% 
  st_transform(crs = "+proj=aeqd +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") 

# Save as RData for plotting 
save(plotting_data, file = "data/final/plotting_data.RData")

# Fix duplicated coordinates

# Remove geometry
r <- plotting_data %>%
  mutate(latitude = st_coordinates(.)[, 2],      # extract latitude
         longitude = st_coordinates(.)[, 1]) %>% # extract longitude
  st_drop_geometry() 

# Check identical coordinates
dup_coord <- r[c(which(duplicated(r$longitude, fromLast = TRUE)),
                 which(duplicated(r$longitude, fromLast = FALSE))),] |>
  group_by(haul_id, year, latitude, longitude) |> summarize(depth = mean(depth))

# Are the duplicates from the original data? - Yes
load("data/intermediate/community_filtered.RData")
check_coord <- community_filtered[community_filtered$haul_id %in% dup_coord$haul_id,
                   c("haul_id", "year", "latitude", "longitude")] 

# This may create problems in the analysis, so longitude is modified 
# by an insignificant amount
sample(r$longitude, 5) # 5 decimals
r[which(duplicated(r$longitude, fromLast = TRUE)),]$longitude <- 
  r[which(duplicated(r$longitude, fromLast = TRUE)),]$longitude + 0.00001

# Check
r[which(duplicated(r$longitude, fromLast = TRUE)),]$longitude

# Save csv
write.table(r, file = "data/final/final_data.csv", sep = ",") 

rm(list = ls())

## End