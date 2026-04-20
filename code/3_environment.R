## Name: 3_environment.R
## Purpose: Add environmental variables (temperature and fishing)
## Author: Filippomaria Cassarino
## Date: 30/03/26

## Notes ----

## Library ----
load("tools/install_load_packages_function.RData")

# Install/load required packages
install_load_packages(
  
  # list of required packages
  c("dplyr",
    "tidyr",
    "ggplot2",
    "reticulate", # Python interface
    "copernicusR",
    "usethis", # edit R profile or environment
    "lubridate", # work with raster dates
    "sf",
    "terra"
  ))

## Download Copernicus data ----

# This section downloads temperature and sea ice data from Copernicus, 
# requiring to register at the link below. 

# Install miniconda (Python) -
# https://www.anaconda.com/docs/getting-started/miniconda/main
install_miniconda()

# Check that miniconda is installed
conda_list()

# Then install the copernicus marine toolbox
reticulate::py_install("copernicusmarine")

# Tell R which ython environment to use
use_condaenv("r-reticulate", required = TRUE)
py_config() # check path and version

# Access the data using your credentials (Copernicus account)
# Register at: https://data.marine.copernicus.eu/register?pk_vid=11ceae9d3160528d17743415560fb363
setup_copernicus(username = "filippomariacassa@gmail.com",
                 password = "v7w!vT_AS8W_S52",
                 store_credentials = TRUE)

# Load original data for extent
load("data/intermediate/community_filtered.RData")

ye <- range(community_filtered$year)
la <- range(community_filtered$latitude)
lo <- range(community_filtered$longitude)

# Download temperature data from Global Ocean Physics Reanalysis
# https://doi.org/10.48670/moi-00021
f <- copernicus_download(
  dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1M-m", # monthly
  variables = c("thetao",  # SST
                "bottomT", # SBT
                "siconc"), # Sea Ice Concentration
  bbox = c(lo[1] - 1, # add a bit of margin (1 degree lat = 111 km)
           lo[2] + 1,
           la[1] - 1,
           la[2] + 1),
  start_date = paste0(ye[1] - 1, "-01-01"),
  end_date = paste0(ye[2], "-12-31"),
  depth = c(0.49402499198913574, 0.49402499198913574), # first depth layer
  dataset_version = "202311", # until end of 2023
)

# Move to data folder (may fail across drives)
file.rename(f, "data/original/temperature_data.nc")

# Read data as raster 
raster_t_data <- terra::rast("data/original/temperature_data.nc")

# Check layers
raster_t_data 

# Visual check
terra::plot(raster_t_data[[1:12]])

# Download temperature data from Global Ocean Colour
# https://doi.org/10.48670/moi-00281
f_2 <- copernicus_download(
  dataset_id = "cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M", # monthly
  variables = c("CHL"), # Chla [mg/m3]
  bbox = c(lo[1] - 1, # add a bit of margin (1 degree lat = 111 km)
           lo[2] + 1,
           la[1] - 1,
           la[2] + 1),
  start_date = paste0(ye[1] - 1, "-01-01"),
  end_date = paste0(ye[2], "-12-31"),
  dataset_version = "202603" 
)

# Move to data folder (may fail across drives)
file.rename(f_2, "data/original/chlorophyll_data.nc")

# Read data as raster
chla_stack <- rast("data/original/chlorophyll_data.nc")

# Check layers
chla_stack

# Visual check
terra::plot(chla_stack[[1:12]])

rm(list = ls())

## Download Global Fishing Watch data ----

# This section downloads fishing effort data from Global Fishing Watch, 
# requiring to register and request an api token (see below)

# Install gfwr package
remotes::install_github("GlobalFishingWatch/gfwr",
                        dependencies = TRUE)

# To access the API, a token is needed. It can be requested here:
# https://gateway.api.globalfishingwatch.org/v3/auth?client=gfw&callback=https%3A%2F%2Fglobalfishingwatch.org%2Four-apis%2Ftokens&locale=en

# Then, the token must be saved to the .Renviron file using 
usethis::edit_r_environ() # write GFW_TOKEN="YOUR_TOKEN" and save .Renviron

# Restart R for changes to take effect

# Check that the token variable has been set
Sys.getenv("GFW_TOKEN") # should return the token

library(gfwr)

# Create the object key with the token
key <- gfw_auth()

# Create an sf object to define the spatial extent of the data
load("data/intermediate/community_filtered.RData")

ye <- range(community_filtered$year)
la <- range(community_filtered$latitude)
lo <- range(community_filtered$longitude)

polygon <- sf::st_bbox(c(xmin = lo[1] - 1,
                         xmax = lo[2] + 1,
                         ymin = la[1] - 1,
                         ymax = la[2] + 1),
                       crs = 4326) |>
  sf::st_as_sfc() |>
  sf::st_as_sf()

# Only one year of data can be imported at a time, so a loop is necessary
years <- (ye[1] - 1):ye[2]

fishing_list <- vector("list", length(years))

for (i in seq_along(years)) {
  
  yr <- years[i]
  
  fishing_list[[i]] <- gfw_ais_fishing_hours(
    spatial_resolution = "HIGH",    # 0.01 degree or ~ 1km
    temporal_resolution = "MONTHLY",
    group_by = "GEARTYPE",
    start_date = paste0(yr, "-01-01"),
    end_date = paste0(yr, "-12-31"),
    region_source = "USER_SHAPEFILE",
    region = polygon
  )
}

# Bind all data frames
fishing_data <- do.call(rbind, fishing_list)

# Save
save(fishing_data, file = "data/original/fishing_data.RData")

# What gears are used in the area?
prop.table(table(fishing_data$geartype))

# Keep only trawling and prepare for raster format
trawling_data <- fishing_data %>%
  filter(geartype == "trawlers") %>% # focus on trawling (60% of observations)
  mutate(time = as.Date(paste0(`Time Range`, "-01"))) %>% # prepare time column
  group_by(time, Lat, Lon) %>%
  summarize(fishing_hours = sum(`Apparent Fishing Hours`, na.rm = TRUE),
            .groups = "drop") 

# Prepare a raster grid to extract values
r_grid <- rast(
  xmin = min(trawling_data$Lon),
  xmax = max(trawling_data$Lon),
  ymin = min(trawling_data$Lat),
  ymax = max(trawling_data$Lat),
  resolution = 0.01  # fishing_data resolution
)

# Split data by time (month and year)
time_list <- split(trawling_data, trawling_data$time) # these will be the layers

# Rasterize each time layer
r_list <- lapply(time_list, function(df) {
  
  pts <- vect(df, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  
  rasterize(
    pts,
    r_grid,
    field = "fishing_hours",
    fun = "sum"   # sum if multiple points per cell
  )
})

# Stack the layers and name them
fishing_stack <- rast(r_list)

time(fishing_stack) <- as.Date(names(time_list))

# Check
fishing_stack

terra::plot(log(fishing_stack[[100:106]]))

# Save intermediate
writeRaster(fishing_stack, file = "data/intermediate/fish_stack.tiff")

rm(list = setdiff(ls(), "fishing_stack"))

## Compute mean, sd, and sum values ----

# This section computes statistics of environmental variables and extracts 
# their raster values at the haul positions

# Load rasters 
temperature_stack <- terra::rast("data/original/temperature_data.nc")
chla_stack <- terra::rast("data/original/chlorophyll_data.nc")
fishing_stack <- terra::rast("data/intermediate/fish_stack.tiff")

# Subset temperature data to compute each variable on its own
sst_stack <- temperature_stack[[grep("thetao", names(temperature_stack))]] # sst
sbt_stack <- temperature_stack[[grep("bottomT", names(temperature_stack))]]# sbt
sic_stack <- temperature_stack[[grep("siconc", names(temperature_stack))]] # sic

# Load community data for haul locations in space and time
load("data/intermediate/community_filtered.RData")

# Output data frame
environmental_data <- community_filtered %>%
  group_by(haul_id, year, latitude, longitude) %>%
  summarize(depth = median(depth), .groups = "drop")

# CRS for extraction (equal area)
y_center <- round(mean(range(environmental_data$latitude, na.rm = TRUE)))
x_center <- round(mean(range(environmental_data$longitude, na.rm = TRUE)))

equal_crs <- paste0("+proj=laea +lat_0=",
                    y_center,
                    " +lon_0=",
                    x_center,
                    " +datum=WGS84 +units=km +no_defs")

# Haul point positions for extraction
haul_points <- environmental_data %>% 
  vect(geom = c("longitude", "latitude"),
       crs = "EPSG:4326") %>% # should be the correct crs
  project(equal_crs)

# Month to divide years
earliest_sampling_month <- min(as.numeric(unique(community_filtered$month)))

# Loop to compute yearly mean and sd, project, and extract
stacks <- list(
  sst  = "sst_stack",
  sbt  = "sbt_stack",
  sic  = "sic_stack",
  chla = "chla_stack"
)

for (v in names(stacks)) {
  
  r <- get(stacks[[v]])
  
  # Shift year to identify the 12 months before sampling ----
  dates <- time(r)
  yrs   <- as.integer(format(dates, "%Y"))
  mons  <- as.integer(format(dates, "%m"))
  
  yrs_new <- ifelse(mons >= earliest_sampling_month, # > or >= ??
                    yrs + 1, yrs) 
  
  # Compute yearly mean and sd ----
  r_mean <- terra::tapp(r, index = yrs_new, fun = mean, na.rm = TRUE)
  r_sd   <- terra::tapp(r, index = yrs_new, fun = sd,   na.rm = TRUE)
  
  # Assign time
  terra::time(r_mean) <- as.numeric(gsub("X", "", names(r_mean)))
  terra::time(r_sd)   <- as.numeric(gsub("X", "", names(r_sd)))
  
  # Project to an equal area projection centered on the data ----
  r_mean <- project(r_mean, equal_crs, res = 4) # 4km resolution
  
  r_sd <- project(r_sd, equal_crs, res = 4) # 4km resolution
  
  mean_vals <- terra::extract(r_mean, haul_points, 
                              method = "bilinear") # 8km resolution
  sd_vals   <- terra::extract(r_sd, haul_points,
                              method = "bilinear") # 8km resolution
  
  # Remove ID column from extract
  mean_vals <- mean_vals[, -1]
  sd_vals   <- sd_vals[, -1]
  
  # Match raster layers to haul year ----
  layer_years <- as.numeric(gsub("X", "", names(r_mean)))
  idx <- match(environmental_data$year, layer_years)
  
  # Add new variable to data ----
  environmental_data[[paste0(v, "_mean")]] <- mean_vals[cbind(1:nrow(mean_vals),
                                                              idx)]
  environmental_data[[paste0(v, "_sd")]] <- sd_vals[cbind(1:nrow(sd_vals),
                                                          idx)]
}

# Add fishing (same as the loop, but just sum)
r <- fishing_stack

dates <- terra::time(r) # shift years
yrs   <- as.integer(format(dates, "%Y"))
mons  <- as.integer(format(dates, "%m"))
yrs_new <- ifelse(mons >= earliest_sampling_month, yrs + 1, yrs) 

r_sum <- terra::tapp(r, index = yrs_new, fun = sum, na.rm = TRUE) # compute sum

terra::time(r_sum) <- as.numeric(gsub("X", "", names(r_sum))) # assign time

r_sum <- project(r_sum, equal_crs, res = 8) # project - 8km resolution

sum_vals <- terra::extract(r_sum, haul_points, # extract
                           method = "simple") # 8km resolution

sum_vals <- sum_vals[, -1] # remove id

layer_years <- as.numeric(gsub("X", "", names(r_sum))) # match raster to haul
idx <- match(environmental_data$year, layer_years)

environmental_data[["sum_fishing"]] <- sum_vals[cbind(1:nrow(sum_vals),
                                                      idx)]
rm(list = setdiff(ls(), c("environmental_data",
                          "sst_stack",
                          "sic_stack",
                          "sbt_stack",
                          "chla_stack",
                          "fishing_stack"))) 

## Clean and save ----

# Final cleanup 
environmental_data <- environmental_data %>%
  select(!chla_sd) %>% # chla_sd makes little sense
  mutate(
    across(c(sic_mean, sic_sd), # NA sic are actually 0
           ~ replace_na(., 0)),
    sum_fishing = if_else(   # fishing NAs from 2012 are 0s
      year >= 2012,
      replace_na(sum_fishing, 0),  
      sum_fishing)
    ) 

# Save
save(environmental_data, file = "data/intermediate/environmental_data.RData")

rm(list = ls())

## End