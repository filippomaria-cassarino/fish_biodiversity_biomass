## File: biomass.R
## Purpose: clean community data, compute and visualize biomass 
#  plot their value over the study area
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

## Library ---- 
load("tools/install_load_packages_function.RData")

# Install/load required packages
install_load_packages(
  c(
    "dplyr",
    "tidyr",
    "ggplot2",
    "ggbreak", # breaks in x axis scale
    "tibble",
    "janitor",
    "sf"
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
length(unique(nor_full$taxon))   # 108 taxa

# Check haul_dur and depth distribution 
haul_unique <- nor_full %>%
  distinct(haul_id, haul_dur, depth)

hist(haul_unique$haul_dur, breaks = 100) # haul_dur should be ~ 15 minutes
range(haul_unique$haul_dur) # 3.8 to 58 minutes
sd(haul_unique$haul_dur) # 3.9 minutes
hist(haul_unique$depth, breaks = 100) #  focus on depth below 500

# Exclude large deviations from 15 minutes and non-shelf areas
community_filtered <- nor_full %>% 
  filter(  
    between(haul_dur, 10, 20), 
    depth < 500) 

length(unique(community_filtered$haul_id)) # 3452 hauls
length(unique(community_filtered$taxon))   # 106 taxa

# Study spatial and temporal range
range(community_filtered$year)
range(community_filtered$latitude)
range(community_filtered$longitude)

y_center <- round(mean(range(community_filtered$latitude, na.rm = TRUE)))
x_center <- round(mean(range(community_filtered$longitude, na.rm = TRUE)))

s <- community_filtered  %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = paste0("+proj=laea +lat_0=", # equal area projection
                                y_center,
                                " +lon_0=",
                                x_center,
                                " +datum=WGS84 +units=km +no_defs")) %>%
  st_union() %>%
  st_convex_hull() %>%
  st_area # 962217.3 [km^2]

## Save filtered community data ----
save(community_filtered, file = "data/intermediate/community_filtered.RData")

rm(list = setdiff(ls(), "community_filtered")) # cleaning

## Check if Sebastes species may inflate richness  ----

# Hauls with Sebastes genus
sebastes <- community_filtered %>%
  group_by(haul_id) %>%
  summarise(
    has_sebastes_genus = "Sebastes" %in% taxon) %>%
  filter(has_sebastes_genus) %>%
  summarise(n_hauls = n())

# Number of hauls containing any sebastes
has_any_sebastes <- community_filtered %>%
  group_by(haul_id) %>%
  summarise(
    has_sebastes = any(
      taxon %in% c(
        "Sebastes",
        "Sebastes viviparus",
        "Sebastes norvegicus",
        "Sebastes mantella"
      ))) %>%
  filter(has_sebastes) %>%
  summarise(n_hauls = n())

# Number of hauls containing both Sebastes genus and species 
overlap_sebastes <- community_filtered %>%
  group_by(haul_id) %>%
  summarise(
    has_sebastes_genus = "Sebastes" %in% taxon,
    has_species_species = any(
      taxon %in% c(
        "Sebastes viviparus",
        "Sebastes norvegicus",
        "Sebastes mantella"
      ))) %>%
  filter(has_sebastes_genus & has_species_species) %>%
  summarise(n_hauls = n())

# % of hauls where richness may be inflated by sebastes
overlap_sebastes/length(unique(community_filtered$haul_id)) * 100 # 9.7%

# % of hauls where any sebastes occurs where there may be richness inflation
overlap_sebastes/has_any_sebastes * 100 # 21.1% 

## Compute and visualize biomass ----

# Compute total fish biomass per unit of area for each haul, then identify the
# most abundant species and dominance structure

# Unique taxa 
length(unique(community_filtered$taxon)) # n = 106

# Individuals encountered by taxon
rare_species <- community_filtered %>%    
  group_by(taxon, rank) %>%
  summarize(number = sum(num, na.rm = TRUE)) %>%
  arrange(number)

table(rare_species$rank)

# Number of individuals resolved to each taxonomic level
number_at_species <- rare_species %>%
  group_by(rank)%>%
  summarise(sum = sum(number)/sum(rare_species$number) * 100)

# Number of greenland sharks
rare_species[rare_species$taxon == "Somniosus microcephalus", ] # 11

# Total biomass per haul 
total_biomass <- community_filtered %>%    
  group_by(haul_id) %>%
  summarize(total_biomass = sum(wgt_cpua, na.rm = TRUE), .groups = "drop") 

# Inspect 0 biomass haul
no_biomass <- filter(community_filtered,
                     haul_id == "2006702/4/2006/1173/2-25-2025")

# Just 1 Macrourus berglax and no weight data. Remove haul. 

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
ggsave(p, file = "figures/dominance_figure.png", unit = "cm", 
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

save(top_10, file = "data/intermediate/top_10_species.RData")

# Top 10 proportion of total
sum(head(spe_biomass$biomass, 10)) / sum(spe_biomass$biomass) # 90% 

# Cod proportion of total
sum(head(spe_biomass$biomass, 1)) / sum(spe_biomass$biomass) # 31% 

# Cod and haddock proportion of total
sum(head(spe_biomass$biomass, 2)) / sum(spe_biomass$biomass) # 44% 

# Merge to have total biomass and top 10 species by haul
b <- species_biomass %>%
  select(c("haul_id", "total_biomass", all_of(top_10))) %>%
  filter(!haul_id == "2006702/4/2006/1173/2-25-2025") 

# Compute total biomass when progressively excluding top 10 species
new_res <- paste0("remove_", 1:10, "_biomass")

for (i in 1:length(new_res)) {
  b[[new_res[i]]] <- b[["total_biomass"]] - rowSums(b[top_10[1:i]])
}

biomass <- b %>%
  janitor::clean_names()

## Save biomass data ----
save(biomass, file = "data/intermediate/biomass.RData")

rm(list = ls())

## End