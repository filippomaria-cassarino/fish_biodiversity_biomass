## File: 7_diversity_models.R
## Purpose: model the response of diversity to warming and fishing
## Author: Filippomaria Cassarino
## Date: 
## ------------------------------------------------------------------------ ----
## Notes ----

## Model checks

# Teve general good (dispersion)
# Teve fishing good (res vs fitted, dispersion)

# Tric general good (dispersion)
# Tric fishing good (res vs fitted, dispersion)

# Tdiv general good (dispersion)
# Tdiv fishing good (dispersion)

# Feve general good (dispersion)
# Feve fishing good (dispersion, res vs pred, a bit too high spt cor)

# Fric general good (dispersion)
# Fric fishing good (res vs fitted maybe nonlinear)

# Fdis general good
# Fdis fishing good (positive trend in res vs fitted, dispersion)

## Library ---- 

# Import data, objects and functions
start_data <- read.csv("data/final/final_data.csv")
load("tools/model_selection_function.RData")
load("tools/model_validation_function.RData")
load("tools/model_prediction_function.RData")
load("tools/vif_function.RData")
load("tools/install_load_function.RData")

# Load required packages
required_pakages <- c("dplyr",
                      "ggplot2",
                      "tidyr",   
                      "sdmTMB", # model fitting
                      "DHARMa", # model validation
                      "sf") 

install_load_function(required_pakages)

# Objects to keep when cleaning
keep <- c("install_load_function",
          "manual_model_selection",
          "model_prediction",          
          "model_selection",
          "model_selection_only_random",
          "model_validation",
          "start_data",
          "vif",
          "general_data",
          "fishing_data",
          "keep")

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## General and fishing data ----

# Prepare data for all responses excluding fishing
general_data <- start_data %>%
  select(c(Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, year, latitude, longitude),
                base::scale))

# Collinearity X (based on biomass model)
general_data %>% select(
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth) %>% vif() # All VIF values < 2

# Prepare data for all responses including fishing
fishing_data <- start_data %>%
  select(c(Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           log_fishing_effort,      # fishing
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, year, latitude, longitude),
                base::scale))

# Collinearity X (based on biomass model)
fishing_data %>% select(
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth,
  log_fishing_effort) %>% vif() # All VIF values < 2

## ------------------------------------------------------------------------ ----
## Teve general model ----

# This is the model for Teve including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(general_data$Teve)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection_only_random(
  directory = "models/Teve",
  name = "general_Teve",
  data = general_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Teve ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth)

# Validation
model_validation(model = general_model,
                 data = general_data,
                 name = "general_Teve",
                 directory = "models/Teve")

# Prediction
#model_prediction(model = general_model, name = "general_Teve", directory = "models/Teve")

# Cleaning
rm(list = setdiff(ls(), keep))

## Teve fishing model ----

# This is the model for Teve including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(fishing_data$Teve)) # use gaussian distribution in modelling

# Model selection - it requires manual intervention. iid is selected
fishing_model <- manual_model_selection(
  directory = "models/Teve",
  name = "fishing_Teve",
  data = fishing_data,
  mesh_cutoff = 65, 
  random_selection = "spatiotemporal_random_fields_iid",
  family = gaussian(link = "identity"),
  fixed_formula = Teve ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_Teve",
                 directory = "models/Teve")

# Prediction
#model_prediction(model = fishing_model, name = "fishing_Teve", directory = "models/Teve")

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Tric general model ----

# This is the model for Tric including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(general_data$Tric)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection_only_random(
  directory = "models/Tric",
  name = "general_Tric",
  data = general_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Tric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth)

# Validation
model_validation(model = general_model,
                 data = general_data,
                 name = "general_Tric",
                 directory = "models/Tric")

# Prediction
#model_prediction(model = general_model, name = "general_Tric", directory = "models/Tric")

# Cleaning
rm(list = setdiff(ls(), keep))

## Tric fishing model ----

# This is the model for Tric including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(fishing_data$Tric)) # use gaussian distribution in modelling

# Model selection 
fishing_model <- model_selection_only_random(
  directory = "models/Tric",
  name = "fishing_Tric",
  data = fishing_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Tric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_Tric",
                 directory = "models/Tric")

# Prediction
#model_prediction(model = fishing_model, name = "fishing_Tric", directory = "models/Tric")

# Cleaning
rm(list = setdiff(ls(), keep))


## ------------------------------------------------------------------------ ----
## Tdiv general model ----

# This is the model for Tdiv including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(general_data$Tdiv)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection_only_random(
  directory = "models/Tdiv",
  name = "general_Tdiv",
  data = general_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Tdiv ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth)

# Validation
model_validation(model = general_model,
                 data = general_data,
                 name = "general_Tdiv",
                 directory = "models/Tdiv")

# Prediction
#model_prediction(model = general_model, name = "general_Tdiv", directory = "models/Tdiv")

# Cleaning
rm(list = setdiff(ls(), keep))

## Tdiv fishing model ----

# This is the model for Tdiv including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(fishing_data$Tdiv)) # use gaussian distribution in modelling

# Model selection 
fishing_model <- model_selection_only_random(
  directory = "models/Tdiv",
  name = "fishing_Tdiv",
  data = fishing_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Tdiv ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_Tdiv",
                 directory = "models/Tdiv")

# Prediction
#model_prediction(model = fishing_model, name = "fishing_Tdiv", directory = "models/Tdiv")

# Cleaning
rm(list = setdiff(ls(), keep))
## ------------------------------------------------------------------------ ----
## Feve general model ----

# This is the model for Feve including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(general_data$Feve)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection_only_random(
  directory = "models/Feve",
  name = "general_Feve",
  data = general_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Feve ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth)

# Validation
model_validation(model = general_model,
                 data = general_data,
                 name = "general_Feve",
                 directory = "models/Feve")

# Prediction
#model_prediction(model = general_model, name = "general_Feve", directory = "models/Feve")

# Cleaning
rm(list = setdiff(ls(), keep))

## Feve fishing model ----

# This is the model for Feve including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(fishing_data$Feve)) # use gaussian distribution in modelling

# Model selection - required manual intervention
fishing_model <- manual_model_selection(
  directory = "models/Feve",
  name = "fishing_Feve",
  data = fishing_data,
  mesh_cutoff = 65, 
  random_selection = "spatiotemporal_random_fields_iid_with_anisotropy",
  family = gaussian(link = "identity"),
  fixed_formula = Feve ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_Feve",
                 directory = "models/Feve")

# Prediction
#model_prediction(model = fishing_model, name = "fishing_Feve", directory = "models/Feve")

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Fric general model ----

# This is the model for Fric including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(general_data$Fric)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection_only_random(
  directory = "models/Fric",
  name = "general_Fric",
  data = general_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Fric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth)

# Validation
model_validation(model = general_model,
                 data = general_data,
                 name = "general_Fric",
                 directory = "models/Fric")

# Prediction
#model_prediction(model = general_model, name = "general_Fric", directory = "models/Fric")

# Cleaning
rm(list = setdiff(ls(), keep))

## Fric fishing model ----

# This is the model for Fric including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(fishing_data$Fric)) # use gaussian distribution in modelling

# Model selection 
fishing_model <- model_selection_only_random(
  directory = "models/Fric",
  name = "fishing_Fric",
  data = fishing_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Fric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_Fric",
                 directory = "models/Fric")

# Prediction
#model_prediction(model = fishing_model, name = "fishing_Fric", directory = "models/Fric")

# Cleaning
rm(list = setdiff(ls(), keep))
## ------------------------------------------------------------------------ ----
## Fdis general model ----

# This is the model for Fdis including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(general_data$Fdis)) # use gaussian distribution in modelling

# Model selection
general_model <- model_selection_only_random(
  directory = "models/Fdis",
  name = "general_Fdis",
  data = general_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Fdis ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth)

# Validation
model_validation(model = general_model,
                 data = general_data,
                 name = "general_Fdis",
                 directory = "models/Fdis")

# Prediction
#model_prediction(model = general_model, name = "general_Fdis", directory = "models/Fdis")

# Cleaning
rm(list = setdiff(ls(), keep))

## Fdis fishing model ----

# This is the model for Fdis including the most data (all areas and time, no
# fishing effort) and selected for the random effect only

# Response distribution
plot(density(fishing_data$Fdis)) # use gaussian distribution in modelling

# Model selection 
fishing_model <- model_selection_only_random(
  directory = "models/Fdis",
  name = "fishing_Fdis",
  data = fishing_data,
  mesh_cutoff = 65, 
  family = gaussian(link = "identity"),
  fixed_formula = Fdis ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_Fdis",
                 directory = "models/Fdis")

# Prediction
#model_prediction(model = fishing_model, name = "fishing_Fdis", directory = "models/Fdis")

# Cleaning
rm(list = setdiff(ls(), keep))
## ------------------------------------------------------------------------ ----
