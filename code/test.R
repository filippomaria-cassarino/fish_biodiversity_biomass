
### Find out if there is a mistake in Teve

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

# Exlcude deviations from 15 minutes and non-shelf areas
community_filtered <- nor_full %>% filter(  
  between(haul_dur, 10, 20), 
  depth < 500)                

# Save filtered community data
save(community_filtered, file = "data/test/community_filtered.RData")

rm(list = setdiff(ls(), "community_filtered")) # cleaning

## Compute biomass 

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

top_10 <- community_filtered %>%
  group_by(accepted_name) %>%
  summarize(biomass = sum(wgt_cpua, na.rm = TRUE)) %>%
  arrange(desc(biomass)) %>%
  slice_head(n = 10) %>%
  mutate(fraction = biomass/sum(community_filtered$wgt_cpua, na.rm = TRUE))

# Merge and keep only top 10 species
biomass <- total_biomass %>%
  left_join(species_biomass, by = "haul_id") %>%
  select(c("haul_id", "total_biomass", top_10$accepted_name)) %>%
  janitor::clean_names()

# Save clean biomass data
save(biomass, file = "data/test/biomass.RData")

# Load filtered community data and trait data
load("data/test/community_filtered.RData")

# Taxonomic richness, evenness (Pielou's), and diversity (Shannon-Wiener)
taxonomic_diversity <- community_filtered %>%
  group_by(haul_id) %>%
  summarize(Tric = vegan::specnumber(accepted_name),
            Teve = vegan::diversity(num_cpua)/log(specnumber(accepted_name)),
            Tdiv = vegan::diversity(num_cpua))

save(taxonomic_diversity, file = "data/test/taxonomic_diversity.RData")

# Haul positions
load("data/test/community_filtered.RData")

haul_point <-  community_filtered %>%
  group_by(haul_id, latitude, longitude, year) %>%
  summarize(depth = mean(depth)) %>%  # this allows to have one raw per haul 
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326) %>%            # transform to point vector (sf object)
  st_transform(crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs")   

save(haul_point, file = "data/test/haul_point.RData")

#
## Copernicus files ----

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

#
## Compute yearly mean and sd values 

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
    "data/test/", vars[i], "_mean.tif"), overwrite = TRUE)  # save raster
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
    "data/test/", vars[i], "_sd.tif"), overwrite = TRUE)  # save raster
  rm(y, t, m, unique_y)
}

rm(sst_stack, sbt_stack, sic_stack, chla_stack, i, stacks, vars) 

## Extract raster values at haul level 

# Load haul positions
load("data/test/haul_point.RData")

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
    "data/test/", vars[k], "_mean.tif"))
  
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
    "data/test/", vars[k], "_sd.tif"))
  
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
save(warming, file = "data/test/warming.RData")

rm(list = ls())

## Merge data ----
# Biomass
load("data/test/biomass.Rdata")

# Taxonomic diversity
load("data/test/taxonomic_diversity.Rdata")

# Haul position, depth and warming variables
load("data/test/warming.RData")

# Functional diversity
load("data/test/functional_diversity.Rdata")

## Join data, correct NAs and duplicated coordinates ----

# Join, correct NAs and project to equidistant
plotting_data <- biomass %>% 
  left_join(taxonomic_diversity, by = "haul_id") %>% 
  left_join(functional_diversity, by = "haul_id") %>%
  left_join(warming, by = "haul_id") %>%
  mutate(
    log_biomass = log(total_biomass),
    log_cod_biomass = log(gadus_morhua + 0.00001),
    across(c(sic_mean, sic_sd), ~ replace_na(., 0)))%>%
  mutate(
         log_chla_mean = log(chla_mean)) %>%
  st_as_sf(crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") %>% 
  st_transform(crs = "+proj=aeqd +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") 


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


# This may create problems in the analysis, so longitude is modified 
# by an insignificant amount
sample(r$longitude, 5) # 5 decimals
r[which(duplicated(r$longitude, fromLast = TRUE)),]$longitude <- 
  r[which(duplicated(r$longitude, fromLast = TRUE)),]$longitude + 0.00001

# Check
r[which(duplicated(r$longitude, fromLast = TRUE)),]$longitude

# Save csv
write.table(r, file = "data/test/final_data.csv", sep = ",") 

rm(list = ls())


## Compare datasets ----

old_data <- read.csv("data/test/final_data.csv", sep = ",")
new_data <- read.csv("data/final/final_data.csv", sep = ",")

evenness_data <- new_data %>% rename(
  "New_Fric" = "kde_fric",
  "New_Feve" = "kde_feve",
  "Teve" = "teve",
  "Tric" = "tric",
  "P_Teve" = "p_teve"
) %>%
  left_join(select(old_data, "haul_id", "Old_Fric" = "Fric", "Old_Feve" = "Feve"))


plot(evenness_data$Old_Fric, evenness_data$P_Teve)

plot(evenness_data$New_Fric, evenness_data$P_Teve)

## Test through models ----


# Prepare data
b_data <- evenness_data

sp <- c("gadus_morhua",
        "melanogrammus_aeglefinus",
        "sebastes_mentella",
        "micromesistius_poutassou",
        "hippoglossoides_platessoides",
        "mallotus_villosus",
        "pollachius_virens",
        "boreogadus_saida",
        "trisopterus_esmarkii",
        "reinhardtius_hippoglossoides")

new_res <- paste0("remove_", 1:10)

for (i in 1:length(new_res)) {
  b_data[[new_res[i]]] <- b_data[["total_biomass"]] - rowSums(b_data[sp[1:i]])
}

biomass_data <- b_data %>%
  select(c(total_biomass, all_of(new_res), # response
           Tric, Teve, P_Teve,     
           Old_Fric, Old_Feve,
           New_Fric, New_Feve,
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean,  log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(starts_with("remove_"), ~ .x + 0.0001), # a few observations are slightly negative due to rounding
         across(-c(haul_id, total_biomass, all_of(new_res), year, latitude, longitude),
                base::scale))

# Models
load("tools/model_selection_functions.RData")

# Test one model
model <- manual_random_model_selection(
  directory = "data/test",
  name = "test_model",
  data = biomass_data,
  random_selection = "spatiotemporal_random_fields_ar1",
  mesh_cutoff = 100, 
  family = lognormal(link = "log"),
  fixed_formula = remove_8 ~  
    #P_Teve +
    Teve +
    Tric +
    New_Feve +
    New_Fric +
    #Old_Fric +
    #Old_Feve +
    sbt_mean +
    sbt_sd +
    sst_sd +
    sic_mean +        
    log_chla_mean +
    depth
)

summary(model)


for (i in 1:length(new_res)) {
  
  name <- paste0(new_res[i], "_biomass")
  
  # Fitting
  model <- random_model_selection(
    directory = "data/test",
    name = name,
    data = biomass_data,
    mesh_cutoff = 100, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i],
      "  ~  
      Tric +
      P_Teve + 
      Old_Fric +
      Old_Feve +
      sbt_mean +
      sbt_sd +
      sst_sd +
      sic_mean +        
      log_chla_mean +
      depth")))


}

## Plot ----

# Import csv files of model coefficients
csv_name <- paste0("remove_", 2:10, "_biomass_model_coefficients.csv")

species <- read.csv(
  file = "data/test/remove_1_biomass_model_coefficients.csv", sep = ",")

for (i in 1:length(csv_name)) { # bind all the others
  data <- read.csv(file = paste0("data/test/",
                                 csv_name[i]), sep = ",")
  species <- bind_rows(species, data)
  
  rm(data)
}

# Data to plot
plot_data <- species %>%
  filter(term != "(Intercept)") %>%
  drop_na() %>%
  mutate(across(-c(model, term), as.numeric))

plot_data$model <- factor(plot_data$model, 
                          levels = c(paste0("remove_", 1:10, "_biomass")))

levels(plot_data$model) <- c("Pielou model", 
                             paste0("Top ", 1:10, " removed"))

plot_data$term <- factor(plot_data$term, levels = c(
  "depth", 
  "log_chla_mean",
  "sst_sd",
  "sbt_sd",
  "sic_mean",
  "sbt_mean",
  "Fric",
  "Feve",
  "Tric",
  "Teve"
))

# Plot
cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(10)

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = cols) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .9),
                     breaks = seq(-.5, .75, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 1.75,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_pielou_species_effect.png",
       unit = "cm", height = 20, width = 25) 

rm(list = setdiff(ls(), keep))
