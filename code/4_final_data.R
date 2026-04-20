## File: 4_final_data.R
## Purpose: join the variables into a single data set and explore for analysis
## Author: Filippomaria Cassarino
## Date: 30 Mar 2026

## Notes ----

## Library ---- 
load("tools/install_load_packages_function.RData")

# Install/load required packages
install_load_packages(
  c(
    "dplyr",
    "sf"
  ))

## Load data ----

# Total and top 10 species biomass by haul
load("data/intermediate/biomass.RData")

# Taxonomic diversity by haul
load("data/intermediate/taxonomic_diversity.RData")

# Functional diversity by haul
load("data/intermediate/functional_diversity.RData")

# Temperature, fishing, depth, year, latitude and longitude by haul
load("data/intermediate/environmental_data.RData")

## Join data ----
full_data <- biomass %>% 
  left_join(functional_diversity, by = "haul_id") %>%
  left_join(taxonomic_diversity, by = "haul_id") %>%
  left_join(environmental_data, by = "haul_id") 

## Check distributions ----

# Variables distribution 
vars <- setdiff(colnames(full_data), c(
  "year",
  "latitude",
  "longitude",
  "haul_id"
))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))

for(i in seq_along(vars)) {     # biomass, chla, and fishing are skewed
  boxplot(full_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i])
}

par(mfrow = c(1, 1))

# Transform where necessary
full_data_2 <- full_data %>% 
  mutate(
    log_chla_mean = log(chla_mean),
    log_sum_fishing = log(sum_fishing + .01)
  )

vars <- c("log_chla_mean", "log_sum_fishing")

for(i in seq_along(vars)) {     # biomass, chla, and fishing are skewed
  boxplot(full_data_2[[vars[i]]],
          horizontal = TRUE,
          main = vars[i])
}

## Check duplicate coordinates ----

# Transform coordinates to aequidistant
y_center <- round(mean(range(full_data_2$latitude, na.rm = TRUE)))
x_center <- round(mean(range(full_data_2$longitude, na.rm = TRUE)))

full_data_3 <- full_data_2 %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = paste0("+proj=aeqd +lat_0=", # project to aequidistant
                                y_center,
                                " +lon_0=",
                                x_center,
                                " +datum=WGS84 +units=km +no_defs")) %>%
  mutate(latitude = st_coordinates(.)[, 2],      # extract latitude
         longitude = st_coordinates(.)[, 1]) %>% # extract longitude
  st_drop_geometry()

# Check identical coordinates
dup_coord <- full_data_3[c(
  which(
    duplicated(full_data_3$longitude, fromLast = TRUE)
  ),
  which(
    duplicated(full_data_3$longitude, fromLast = FALSE))
), ] |>
  group_by(
    haul_id,
    year,
    latitude,
    longitude) |> 
  summarize(depth = mean(depth)) # three couples of hauls 

# Are the duplicates from the original data? - Yes
load("data/intermediate/community_filtered.RData")

check_coord <- community_filtered[
  community_filtered$haul_id %in% dup_coord$haul_id,
  c("haul_id", "year", "latitude", "longitude")
] |>
  distinct(haul_id, latitude, longitude)

# Modify identical coordinates by an insignificant amount
full_data_3[
  which(
    duplicated(full_data_3$longitude, fromLast = TRUE)),]$longitude <- 
  full_data_3[
    which(
      duplicated(full_data_3$longitude, fromLast = TRUE)),]$longitude + 1e-5

## Save final data ----
write.table(full_data_3, file = "data/final/final_data.csv", sep = ",") 

rm(list = ls())

## End