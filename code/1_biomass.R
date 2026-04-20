## File: biomass.R
## Purpose: clean community data, compute and visualize biomass 
#  plot their value over the study area
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

## Library ---- 
load("tools/install_load_packages_function.RData")
load("tools/model_selection_functions.RData")

# Install/load required packages
install_load_packages(
  c(
    "dplyr",
    "tidyr",
    "ggplot2",
    "ggbreak", # breaks in x axis scale
    "tibble",
    "GGally",
    "mFD",
    "BAT",
    "janitor",
    "sf",
    "sdmTMB"
  ))

## Community data cleaning ----

# This section downloads and cleans the community data

# Download data from github for 2004-2021
download.file(
  url = "https://github.com/fishglob/FishGlob_data/raw/refs/heads/main/outputs/Cleaned_data/NOR-BTS_clean.RData",
  destfile = "data/original/NOR-BTS_clean.RData",
  mode = "wb")

# Load into environment
load("data/original/NOR-BTS_clean.RData")

# Load data for 2022 if not on github yet (shared directly by Laurene Pecuchet)
load("data/original/BESS_data2022.RData")

# Bind and filter to retain only the Barents Sea Ecosystem Survey data
# https://github.com/AquaAuma/FishGlob_data/tree/main/metadata_docs#norway-nor-bts
nor_full <- data %>% 
  rbind(data2022) %>%                 
  filter(
    gear %in% c("3270","3271"), 
    latitude > 70,              
    month %in% c("08", "09")) %>%
  dplyr::rename("taxon" = "accepted_name") # rename for clarity

length(unique(nor_full$haul_id)) # 3932 hauls
length(unique(nor_full$taxon))   # 108 species

# Check haul_dur and depth distribution 
haul_unique <- nor_full %>%
  distinct(haul_id, haul_dur, depth)

hist(haul_unique$haul_dur, breaks = 100) # haul_dur should be ~ 15 minutes
hist(haul_unique$depth, breaks = 100) #  focus on depth below 500

# Exclude large deviations from 15 minutes and non-shelf areas
community_filtered <- nor_full %>% 
  filter(  
    between(haul_dur, 10, 20), 
    depth < 500) 

length(unique(community_filtered$haul_id)) # 3452 hauls
length(unique(community_filtered$taxon))   # 106 species

## Check the effect of haul duration on species richness ----

# Haul duration may affect the number of species encountered, so here
# tric is models as response to haul duration, including spatiotemporal
# random fields for possible autocorrelations

# Extract midpoint of the spatial distribution
y_center <- round(mean(range(community_filtered$latitude, na.rm = TRUE)))
x_center <- round(mean(range(community_filtered$longitude, na.rm = TRUE)))

# Prepare data
duration_data <- community_filtered %>% 
  group_by(
    haul_id,
    haul_dur,
    year,  
    latitude,
    longitude
  ) %>%
  summarize(tric = length(num_cpua), .groups = "drop") %>% # richness
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = paste0("+proj=aeqd +lat_0=", # project to aequidistant
                                y_center,
                                " +lon_0=",
                                x_center,
                                " +datum=WGS84 +units=km +no_defs")) %>%
  mutate(latitude = st_coordinates(.)[, 2],      # extract latitude
         longitude = st_coordinates(.)[, 1]) %>% # extract longitude
  st_drop_geometry() %>% 
  distinct(latitude, longitude, .keep_all = TRUE)

# Check distributions
hist(duration_data$haul_dur, breaks = 100) 
hist(duration_data$tric, breaks = 100) 

# Fit model
duration_model <- model_selection(
  directory = "models/haul_duration",
  name = "haul_duration",
  data = duration_data,  
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = tric ~ haul_dur
)

# Wald test - H0: coefficient = 0
w <- ((tidy(duration_model)[2, 2] / tidy(duration_model)[2, 3])[1, 1])^2
1 - pchisq(w, df = 1) # retain H0

# Validate
load("tools/model_validation_function.RData")

model_validation(
  model = duration_model,
  name = "duration",
  directory = "models/haul_duration"
)

# Predict
load("tools/model_prediction_function.RData")

model_prediction(
  model = duration_model,
  vars = "haul_dur",
  name = "duration",
  directory = "models/haul_duration"
)

# Model validation shows that assumptions are largely met. 
# haul_dur has a minimal effect

## Save filtered community data ----
save(community_filtered, file = "data/intermediate/community_filtered.RData")

rm(list = setdiff(ls(), "community_filtered")) # cleaning

## Compute and visualize biomass ----

# Compute total fish biomass per unit of area for each haul, then identify the
# most abundan species and dominance structure

# Unique taxa 
length(unique(community_filtered$taxon)) # n = 106

# Individual encountered by taxon
rare_species <- community_filtered %>%    
  group_by(taxon) %>%
  summarize(number = sum(num, na.rm = TRUE)) %>%
  arrange(number)

# Total biomass per haul 
total_biomass <- community_filtered %>%    
  group_by(haul_id) %>%
  summarize(total_biomass = sum(wgt_cpua, na.rm = TRUE), .groups = "drop") 

# Total biomass by species
spe_biomass <- community_filtered %>%    
  group_by(taxon) %>%
  summarize(biomass = sum(wgt_cpua, na.rm = TRUE)) %>%
  arrange(desc(biomass)) 

# Progressively remove top species
no_top_10 <- tail(spe_biomass, nrow(spe_biomass) - 10)
no_top_30 <- tail(spe_biomass, nrow(spe_biomass) - 30)
no_top_50 <- tail(spe_biomass, nrow(spe_biomass) - 50)

new_spe <- bind_rows(spe_biomass,
                     no_top_10, 
                     no_top_30, 
                     no_top_50)

new_spe$exclude <- c(rep(paste0("All species (n = ", nrow(spe_biomass), ")"),
                                nrow(spe_biomass)),
                     rep("Top 10 species excluded", nrow(spe_biomass) - 10),
                     rep("Top 30 species excluded", nrow(spe_biomass) - 30),
                     rep("Top 50 species excluded", nrow(spe_biomass) - 50))

# Visualize dominance structure
p <- ggplot(data = new_spe, aes(x = biomass)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                     fill = "plum", color = "black", bins = 50) +
  facet_wrap(vars(exclude), nrow = 1, scales = "free") +
  labs(x = "\nBiomass (kg/km2)", y = "Proportion of species\n") + # check unit
  theme_bw()+
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid = element_blank())

# Save plot
ggsave(p, file = "figures/biomass_plots/dominance_figure.png", unit = "cm", 
       height = 7, width = 27)

# Species biomass per haul
species_biomass <- community_filtered %>%  
  group_by(
    haul_id,
    year,
    taxon
    ) %>%
  summarize(biomass = sum(wgt_cpua, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = taxon,
    values_from = biomass,
    values_fill = 0  
  ) %>%
  left_join(total_biomass, by = "haul_id")

# Most abundant species 
top_10 <- spe_biomass$taxon %>% head(10)

top_10

# Top 10 proportion of total
sum(head(spe_biomass$biomass, 10)) / sum(spe_biomass$biomass) # 90% 

# Cod proportion of total
sum(head(spe_biomass$biomass, 1)) / sum(spe_biomass$biomass) # 31% 

# Visualize changes in species biomass
plot_data <- species_biomass %>%
  group_by(year) %>%
  summarise(
    across(all_of(top_10), \(x) sum(x, na.rm = TRUE)),
    total_biomass = sum(total_biomass, na.rm = TRUE)
  ) %>%
  mutate(
    Others = total_biomass - rowSums(across(all_of(top_10)), na.rm = TRUE)
  ) %>%
  select(year, all_of(top_10), Others) %>%
  pivot_longer(
    -year,
    names_to = "species",
    values_to = "biomass"
  )

plot_data$species <- factor(plot_data$species, levels = c(top_10, "Others"))

levels(plot_data$species) <- c(
  "G. morhua",
  "M. aeglefinus",
  "S. mentella",
  "M. poutassou",
  "H. platessoides",
  "M. villosus",
  "P. virens",
  "B. saida",
  "T. esmarkii",
  "R. hippoglossoides",
  "Others")

cols <- c(colorRampPalette(c("#2C3E50","#03A9F4", "#E2F4C7"))(10),
          "gray50")

p <- ggplot(
  plot_data,
  aes(
    x = year,
    y = biomass,
    fill = species
  )
) +
  geom_area(position = "fill") +
  geom_line(position = "fill", color = "gray30", linewidth = 0.1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Year",
    y = "Proportion of total biomass",
    fill = "Species"
  ) +
  theme_minimal(base_size = 14)

# Save plot
ggsave(p, file = "figures/biomass_plots/top_species_biomass.png",
       unit = "cm", height = 16, width = 25)

# Merge to have total biomass and top 10 species by haul
biomass <- species_biomass %>%
  select(c("haul_id", "total_biomass", all_of(top_10))) %>%
  janitor::clean_names()

## Save biomass data ----
save(biomass, file = "data/intermediate/biomass.RData")

rm(list = ls())

## Plotting biomass distribution (WORK IN PROGRESS) ----

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
titles <- c("Total biomass", top_10$taxon)

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