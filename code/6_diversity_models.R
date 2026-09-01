## File: 6_diversity_models.R
## Purpose: model the response of diversity to warming and fishing
## Author: Filippomaria Cassarino
## Date: 31 Aug 2026
## Notes ----

## Library ---- 

# Import data, objects and functions
final_data <- read.csv("data/final/final_data.csv")

# Load functions
load("tools/install_load_packages_function.RData")
load("tools/model_selection_functions.RData")
load("tools/model_prediction_function.RData")

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
  "manual_model_selection",
  "model_prediction"
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
    log_chla_sd,
    depth,
    haul_id,
    haul_dur,
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
    log_chla_sd,
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
## Haul duration and species richness ----

# Haul duration may affect the number of species encountered, so here
# tric is models as response to haul duration, including spatiotemporal
# random fields for possible autocorrelations

# Data
haul_data <- final_data %>%
  select(      
    tric,
    haul_dur,
    year, 
    latitude,
    longitude
  ) %>%
  drop_na() 

# Response distribution
plot(density(haul_data$tric)) # use gaussian distribution in modelling

# Fit model
duration_model <- model_selection(
  directory = "models/haul_duration",
  name = "haul_duration",
  data = haul_data,  
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  fixed_formula = tric ~ haul_dur      
)

# Wald test - H0: coefficient = 0
w <- ((tidy(duration_model)[2, 2] / tidy(duration_model)[2, 3])[1, 1])^2
p <- 1 - pchisq(w, df = 1) # retain H0

save(p, file = "models/haul_duration/wald_test_p.RData")

# Predict
model_prediction(
  model = duration_model,
  vars = "haul_dur",
  name = "duration",
  directory = "models/haul_duration"
)

# Model validation shows that assumptions are largely met. 
# haul_dur has a minimal effect

rm(list = setdiff(ls(), keep)) # cleaning

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
  family = student(link = "identity"),
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
  family = student(link = "identity"),
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
  family = student(link = "identity"),
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
  family = student(link = "identity"),
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

csv_files <- csv_files[
  csv_files != "models/diversity/diversity_models_coefficients.csv"
  ]

combined_csv <- bind_rows(lapply(csv_files, read.csv)) |>
  filter(term !="(Intercept)")

write.csv(combined_csv,
          file = "models/diversity/diversity_models_coefficients.csv")

# Cleaning
file.remove(csv_files)

## ------------------------------------------------------------------------ ----
## Check models' underdispersion ----

# Most models show underdispersion, which may be a result of complex
# random fields capturing a lot of the variation in the residuals. 
# We refit them with a simpler structure to verify if underdisperison 
# persists.

# Model selection
general_model <- manual_model_selection(
  directory = "models/diversity",
  name = "general_tric_underdisperison_check",
  data = general_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  random_selection = "spatial_random_fields",
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

file.remove(
  "models/diversity/general_tric_underdisperison_check_model_coefficients.csv")

# The answer is no - underdispersion disappears (or it is greately reduced).
# However, the AIC difference is large (check AIC tables in the selection pdf),
# so the more complex structure must be preferred. 
# The model likely captures a large amount of real structure.

# End
