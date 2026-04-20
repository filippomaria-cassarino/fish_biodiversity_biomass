## File: 6_diversity_models.R
## Purpose: model the response of diversity to warming and fishing
## Author: Filippomaria Cassarino
## Date: 08 Apr 2026
## Notes ----

## Library ---- 

# Import data, objects and functions
final_data <- read.csv("data/final/final_data.csv")

# Load functions
load("tools/model_selection_functions.RData")
load("tools/vif_function.RData")
load("tools/install_load_packages_function.RData")

# Load required packages
required_pakages <- c("dplyr",
                      "ggplot2",
                      "tidyr",   
                      "sdmTMB", 
                      "DHARMa") 

install_load_packages(required_pakages)

# Objects to keep when cleaning
keep <- c(
  "keep",
  "final_data",
  "general_data",
  "fishing_data",
  "vif",
  "model_selection",
  "manual_model_selection"
)

rm(list = setdiff(ls(), keep))

## General and fishing data ----

# Prepare data for all responses excluding fishing
general_data <- final_data %>%
  select(
    kde_fric,
    kde_feve,       
    tric,
    teve,           
    sst_mean,       
    sbt_mean,       
    sic_mean, 
    sbt_sd,
    sst_sd,
    sic_sd, 
    log_chla_mean,
    depth,
    haul_id,
    year, 
    latitude,
    longitude
  ) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# Prepare data for all responses including fishing
fishing_data <- final_data %>%
  select(
    log_sum_fishing,
    kde_fric,
    kde_feve,       
    tric,
    teve,           
    sst_mean,       
    sbt_mean,       
    sic_mean, 
    sbt_sd,
    sst_sd,
    sic_sd, 
    log_chla_mean,
    depth,
    haul_id,
    year, 
    latitude,
    longitude
  ) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# Collinearity and distribution were already inspected in the biomass models

## ------------------------------------------------------------------------ ----
## Taxonomic evenness general model ----

# This is the model for teve including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(general_data$teve)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection(
  directory = "models/diversity",
  name = "general_teve",
  data = general_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = teve ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
  )

# Cleaning
rm(list = setdiff(ls(), keep))

## Taxonomic evenness fishing model ----

# This is the model for teve including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(fishing_data$teve)) # use gaussian distribution in modelling

# Model selection
fishing_model <- model_selection(
  directory = "models/diversity",
  name = "fishing_teve",
  data = fishing_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = teve ~
    log_sum_fishing +
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Taxonomic richness general model ----

# This is the model for tric including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(general_data$tric)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection(
  directory = "models/diversity",
  name = "general_tric",
  data = general_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = tric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## Taxonomic richness fishing model ----

# This is the model for tric including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(fishing_data$tric)) # use gaussian distribution in modelling

# Model selection
fishing_model <- model_selection(
  directory = "models/diversity",
  name = "fishing_tric",
  data = fishing_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = tric ~
    log_sum_fishing +
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Functional evenness general model ----

# This is the model for kde_feve including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(general_data$kde_feve)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection(
  directory = "models/diversity",
  name = "general_kde_feve",
  data = general_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = kde_feve ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## Functional evenness fishing model ----

# This is the model for kde_feve including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(fishing_data$kde_feve)) # use gaussian distribution in modelling

# Model selection
fishing_model <- model_selection(
  directory = "models/diversity",
  name = "fishing_kde_feve",
  data = fishing_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = kde_feve ~
    log_sum_fishing +
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Functional richness general model ----

# This is the model for kde_fric including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(general_data$kde_fric)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection(
  directory = "models/diversity",
  name = "general_kde_fric",
  data = general_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = kde_fric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## Functional richness fishing model ----

# This is the model for kde_fric including the most data (all areas and time, no
# fishing effort) 

# Response distribution
plot(density(fishing_data$kde_fric)) # use gaussian distribution in modelling

# Model selection
fishing_model <- model_selection(
  directory = "models/diversity",
  name = "fishing_kde_fric",
  data = fishing_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = kde_fric ~
    log_sum_fishing +
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Merge all coefficient files in one csv ----
csv_files <- list.files(
  path = "models/diversity",
  pattern = ".csv$",
  full.names = TRUE
)

combined_csv <- bind_rows(lapply(csv_files, read.csv)) |>
  filter(term !="(Intercept)")

write.csv(combined_csv,
          file = "models/diversity/diversity_models_coefficients.csv")

# Cleaning
file.remove(csv_files)

# End