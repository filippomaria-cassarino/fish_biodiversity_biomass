## File: biomass.R
## Purpose: clean and summarize standardized data into biomass and diversity metrics;
#  plot their value over the study area
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

# the mapping objects should be created in a single, seperate script

## Library ---- 

library(dplyr)
library(ggplot2)
library(tidyr) # pivot to wide format
library(terra) # ratser 
library(sf)    # vector
library(rnaturalearth) # land shape
library(modi) # compute cwv
library(janitor) # clean column names 

# Selecting land contours for the plots
land <- st_crop(ne_countries(scale = "medium", returnclass = "sf"),
                xmin = -10, xmax = 65, 
                ymin = 55, ymax = 90) |>
  st_transform(
    crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") 

# Creating a custom theme
theme_custom <- function() {
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray", color = "gray45", linewidth = 1),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 8)
  )
}

# Function to plot mean values
grid_plot <- function(data, variable, land) {
  
  # plotting devices
  load("tools/mapping_objects.RData")
  
  # create grid
  g <- data %>%
    st_make_grid(cellsize = 64.82,  # 35 nautical miles in km
                 what = "polygons",
                 square = TRUE, 
                 offset = c(-483, -573)) %>% # center the cells
    st_sf()
  
  # Spatial join
  grid_r <- st_join(g, data, left = FALSE)
  
  # Plot limits
  xlim <- st_bbox(data)[c(1, 3)] + c(-65, 65)
  ylim <- st_bbox(data)[c(2, 4)] + c(-65, 65)
  
  # Mean per cell
  grid_mean <- grid_r %>%
    group_by(geometry) %>%
    summarize(x = mean(.data[[variable]], na.rm = TRUE)) %>%
    st_as_sf()
  
  # Range
  max_val <- max(grid_mean$x, na.rm = TRUE)
  min_val <- min(grid_mean$x, na.rm = TRUE)
  
  # Plot
  ggplot() +
    geom_sf(data = grid_mean, aes(fill = x), color = NA) +
    geom_sf(data = land, fill = "gray45", color = "gray45") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                         values = scales::rescale(c(min_val, max_val)),
                         limits = c(min_val, max_val),
                         name = variable,
                         breaks = seq(min_val, max_val, length.out = 3)) +
    labs(x = "\nLongitude", y = "Latitude\n", title = "Map") +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 10, ticks.linewidth = .1)) +
    theme_custom()
}

# Save 
save(theme_custom, land, grid_plot, file = "tools/mapping_objects.RData")

#
## Community data cleaning ----
 
# Load data for 2004-2021 (FISHGLOB_data) AND 2022 (Laurene Pecouchet)
load("data/original/NOR-BTS_clean.RData") 
load("data/original/BESS_data2022.RData")

# Bind and filter to retain only the Barents Sea Ecosystem Survey data (https://github.com/AquaAuma/FishGlob_data/tree/main/metadata_docs#norway-nor-bts)
nor_full <- data %>% 
  rbind(data2022) %>%                 
  filter(
    gear %in% c("3270","3271"), 
    latitude > 70,              
    month %in% c("08", "09"))  

length(unique(nor_full$haul_id)) # 3932 hauls

#  Check haul duration and depth distribution
hist(nor_full$haul_dur, breaks = 100) # haul_dur should be ~ 15 minutes
hist(nor_full$depth, breaks = 100)    # continental shelf ends at 500m

# Exlcude deviations from 15 minutes and non-shelf areas
community_filtered <- nor_full %>% filter(  
             between(haul_dur, 10, 20), 
             depth < 500)                

length(unique(community_filtered$haul_id)) # 3452 hauls, 480 were lost from filtering

# Save filtered community data
save(community_filtered, file = "data/intermediate/community_filtered.RData")

rm(list = setdiff(ls(), "community_filtered")) # cleaning


## Compute biomass ----

# Unique taxa 
length(unique(community_filtered$accepted_name)) # n = 106

# Total biomass per haul 
total_biomass <- community_filtered %>%    
  group_by(haul_id) %>%
  summarize(total_biomass = sum(wgt_cpua, na.rm = TRUE)) 

# Species biomass per haul
species_biomass <- community_filtered %>%  
  group_by(haul_id, accepted_name) %>%
  summarize(biomass = sum(wgt_cpua, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = accepted_name,
    values_from = biomass,
    values_fill = 0  
  )

# Most abundant species and proportion of total biomass
top_10 <- community_filtered %>%
  group_by(accepted_name) %>%
  summarize(biomass = sum(wgt_cpua, na.rm = TRUE)) %>%
  arrange(desc(biomass)) %>%
  slice_head(n = 10) %>%
  mutate(fraction = biomass/sum(community_filtered$wgt_cpua, na.rm = TRUE))

sum(top_10$fraction) # together, these 10 account for ~ 90% of the sampled biomass

# Merge and keep only top 10 species
biomass <- total_biomass %>%
  left_join(species_biomass, by = "haul_id") %>%
  select(c("haul_id", "total_biomass", top_10$accepted_name)) %>%
  janitor::clean_names()

# Save clean biomass data
save(biomass, file = "data/intermediate/biomass.RData")

## Plotting biomass distribution ----

# Load theme and land contours
load("tools/mapping_objects.RData")

# Data for the plot
biomass_plot_data <- community_filtered %>%
  group_by(haul_id, year, latitude, longitude) %>%
  summarize(depth = median(depth)) %>%
  left_join(biomass, by = "haul_id") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") 

# Continuous grid to bin the observations (65 km cells)
g <- biomass_plot_data %>%
  st_make_grid(cellsize = 64.82,  # 35 nautical miles in km
               what = "polygons",
               square = TRUE, 
               offset = c(-483, -573)) %>% # center the cells
  st_sf() 

# Spatial join (inner join) to match observations to cells
grid_r <- st_join(g, biomass_plot_data, left = FALSE) 

# Check observation number per cell (ideally = n of years)
cell_counts <- grid_r%>%
  group_by(geometry) %>%
  summarize(n_obs = n(), .groups = "drop")

prop.table(table((cell_counts$n_obs > 5))) # ~ 80% has more than 5 observations

# responses to plot
res <- names(biomass[2:12])
titles <- c("Total biomass", top_10$accepted_name)

# Plot limits
xlim <- st_bbox(biomass_plot_data)[c(1, 3)] + c(-65, 65)
ylim <- st_bbox(biomass_plot_data)[c(2, 4)] + c(-65, 65)

# Compute mean values over the study area, plot and save
for (i in 1:length(res)) { 
  
  grid_mean <- grid_r %>%
    group_by(geometry) %>%
    summarize(x = mean(get(res[i]), na.rm = TRUE)) %>%
    st_as_sf()
  
  max <- round(max(log(grid_mean$x + 1)) * 10)/10
  min <- round(min(log(grid_mean$x + 1)) * 10)/10
  
  p <- ggplot() +
    geom_sf(data = grid_mean, aes(fill = log(x + 1)), color = NA) +
    geom_sf(data = land, fill = "gray45", color = "gray45") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                         values = scales::rescale(c(min, max)),
                         limits = c(min, max), 
                         name = "log(kg/km2)",
                         breaks = seq(min, max, length = 3)) + 
    labs(x = "\nLongitude", y = "Latitude\n",
         title = titles[i]) + 
    guides(fill = guide_colorbar(barwidth = .5,
                                 barheight = 10,
                                 ticks.linewidth = .1)) +
    theme_custom() 
  
  ggsave(p, file = paste0("figures/biomass_plots/", i - 1, "_mean_", res[i], ".png"),
         width = 15, height = 10, units = "cm")
}

# Plot the number of hauls included in each cell and save
p <- ggplot() +
  geom_sf(data = cell_counts, fill = "white", color = "black", size = 0.3) +  
  geom_sf(data = land, fill = "gray45", color = "gray45", alpha = 0.3) +  
  geom_sf_text(data = cell_counts, aes(label = n_obs), size = 2, color = "black") +  
  scale_x_continuous(limits = xlim) +
  scale_y_continuous(limits = ylim) +
  labs(
    x = "\nLongitude", y = "Latitude\n",
    title = "Hauls per Grid Cell"
  ) +
  theme_custom() 

ggsave(p, file = "figures/biomass_plots/haul_number.png",
         width = 15, height = 10, units = "cm")

## End