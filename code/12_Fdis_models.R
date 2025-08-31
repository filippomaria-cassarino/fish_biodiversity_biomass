## File: 12_Fdis_models.R
## Purpose: model the response of taxonomic richness to warming and fishing
## Author: Filippomaria Cassarino
## Date: 
## ------------------------------------------------------------------------ ----
## Notes ----

## Library ---- 

# Load required packages
library(dplyr)       # data manipulation
library(tidyr)       # pivot to wide format 
library(ggplot2)     # plotting
library(sdmTMB)      # model
library(performance) # VIF
library(DHARMa)      # validation
library(writexl)     # save coefficients as .xlsx

# Import data and functions
final_data <- read.csv("data/final/final_data.csv")
load("tools/model_selection_function.RData")
load("tools/model_validation_function.RData")
load("tools/model_prediction_function.RData")

#
## ------------------------------------------------------------------------ ----
## All model no fishing - data exploration ----

# Prepare data of the whole study area (all) for exploration
check_data <- final_data %>%
  select(c(Fdis,                   # response 
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           log_chla_mean,          # productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() 

# Variables distribution and out-liers were assessed in 6_biomass_models.R
boxplot(check_data$Fdis, horizontal = TRUE) # no majour outliers
plot(density(check_data$Fdis)) # use gaussian distribution in modelling

# Collinearity X was assessed in 6_biomass_models.R, these were the correlated
#  variables:
# sst_mean - sbt_mean -> remove sst_mean
# sst_mean - sic_mean -> remove sst_mean
# sic_mean - sic_sd   -> remove sic_sd

# Data for modelling
all_data <- check_data  %>%
  mutate(across(-c(haul_id, year, latitude, longitude),
                base::scale))

## All model no fishing - selection, validation, prediction ----

# Selection
all_model <- model_selection(directory = "models/Fdis/all/no_fishing",
                             name = "all",
                             data = all_data,
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
model_validation(model = all_model,
                 data = all_data,
                 directory = "models/Fdis/all/no_fishing")

# Prediction
model_prediction(model = all_model,
                 directory = "models/Fdis/all/no_fishing")

## Arctic model no fishing - data exploration ----

# Identify arctic data
divide_arc <- quantile(all_data$sst_mean, 0.25, na.rm = TRUE)

# Prepare data for exploration
arctic_data <- all_data[all_data$sst_mean <= divide_arc, ]

# Variables distribution and out-liers were assessed in 6_biomass_models.R
boxplot(arctic_data$Fdis, horizontal = TRUE) # no majour outliers
plot(density(arctic_data$Fdis)) # use gaussian distribution in modelling

# Collinearity of X  was assessed in 6_biomass_models.R, these were the correlated
#  variables:
# sst_mean - sst_sd   -> remove sst_mean
# sst_mean - sic_mean -> remove sst_mean

## Arctic model no fishing - selection, validation, prediction ----

# Selection
arctic_model <- model_selection(directory = "models/biomass/arctic/no_fishing",
                                name = "arctic",
                                data = arctic_data,
                                mesh_cutoff = 65, 
                                family = gaussian(link = "identity"),
                                fixed_formula = Fdis ~ 
                                  sbt_mean +
                                  sbt_sd +
                                  #sst_mean +
                                  sst_sd +
                                  sic_mean + 
                                  sic_sd +
                                  log_chla_mean +
                                  depth)

# Validation
model_validation(model = arctic_model,
                 data = arctic_data,
                 directory = "models/Fdis/arctic/no_fishing")

# Prediction
model_prediction(model = arctic_model,
                 directory = "models/Fdis/arctic/no_fishing")

## Boreal model no fishing - data exploration ----

# Identify boreal data
divide_bor <- quantile(all_data$sst_mean, 0.75, na.rm = TRUE)

# Prepare data
boreal_data <- all_data[all_data$sst_mean >= divide_bor, ]

# Variables distribution and out-liers were assessed in 6_biomass_models.R
boxplot(boreal_data$Fdis, horizontal = TRUE) # no majour outliers
plot(density(boreal_data$Fdis)) # use gaussian distribution in modelling

# Collinearity X was assessed in 6_biomass_models.R, these were the correlated
#  variables:
# sst_mean - sbt_mean -> remove sst_mean
# sbt_mean - depth    -> remove depth
# sic variables can be removed because they are almost all 0s

## Boreal model no fishing  - selection, validation, prediction ----

# Selection
boreal_model <- model_selection(directory = "models/Fdis/boreal/no_fishing",
                                name = "boreal",
                                data = boreal_data,
                                mesh_cutoff = 65, 
                                family = gaussian(link = "identity"),
                                fixed_formula = Fdis ~ 
                                  sbt_mean +
                                  sbt_sd +
                                  #sst_mean +
                                  sst_sd +
                                  #sic_mean + 
                                  #sic_sd +
                                  #depth +
                                  log_chla_mean)

# Validation
model_validation(model = boreal_model,
                 data = boreal_data,
                 directory = "models/Fdis/boreal/no_fishing")

# Prediction
model_prediction(model = boreal_model,
                 directory = "models/Fdis/boreal/no_fishing")


# Cleaning before the model with fishing
rm(list = setdiff(ls(), c("final_data",
                          "model_prediction",
                          "model_selection",
                          "model_validation")))

## ------------------------------------------------------------------------ ----
## All model with fishing - data exploration ----

# Prepare data of the whole study area (all) for exploration
check_data <- final_data %>%
  select(c(Fdis,                   # response
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           log_chla_mean,          # productivity
           log_fishing_effort,     # fishing
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() 

# Variables distribution and outliers were assessed in 6_biomass_models.R
boxplot(check_data$Fdis, horizontal = TRUE) # no majour outliers
plot(density(check_data$Fdis)) # use gaussina distribution in modelling

# Collinearity of X was assessed in 6_biomass_models.R, these were the correlated
#  variables:
# sst_mean - sbt_mean -> remove sst_mean
# sst_mean - sic_mean -> remove sst_mean
# sst_mean - sic_sd   -> remove sst_mean
# sic_mean - sic_sd   -> remove sic_sd

# Data for modelling
all_data <- check_data  %>%
  mutate(across(-c(haul_id, year, latitude, longitude),
                base::scale))

## All model with fishing - selection, validation, prediction ----

# Selection
all_model <- model_selection(directory = "models/Fdis/all/with_fishing",
                             name = "all",
                             data = all_data,
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
                               log_fishing_effort +
                               depth)

# Validation
model_validation(model = all_model,
                 data = all_data,
                 directory = "models/Fdis/all/with_fishing")

# Prediction
model_prediction(model = all_model,
                 directory = "models/Fdis/all/with_fishing")

## Arctic model with fishing - data exploration ----

# Identify arctic data
divide_arc <- quantile(all_data$sst_mean, 0.25, na.rm = TRUE)

# Prepare data for exploration
arctic_data <- all_data[all_data$sst_mean <= divide_arc, ]

# Variables distribution and out-liers were assessed in 6_biomass_models.R
boxplot(arctic_data$Fdis, horizontal = TRUE) # no majour outliers
plot(density(arctic_data$Fdis)) # use gaussian distribution in modelling


# Collinearity of X was assessed in 6_biomass_models.R, these were the correlated
# variables:
# sst_mean - sst_sd   -> remove sst_mean
# sic_mean - sic_sd   -> remove sic_sd

## Arctic model with fishing - selection, validation, prediction ----

# Selection
arctic_model <- model_selection(directory = "models/Fdis/arctic/with_fishing",
                                name = "arctic",
                                data = arctic_data,
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
                                  log_fishing_effort +
                                  depth)

# Validation
model_validation(model = arctic_model,
                 data = arctic_data,
                 directory = "models/Fdis/arctic/with_fishing")

# Prediction
model_prediction(model = arctic_model,
                 directory = "models/Fdis/arctic/with_fishing")

## Boreal model with fishing - data exploration ----

# Identify boreal data
divide_bor <- quantile(all_data$sst_mean, 0.75, na.rm = TRUE)

# Prepare data
boreal_data <- all_data[all_data$sst_mean >= divide_bor, ]

# Variables distribution and out-liers were assessed in 6_biomass_models.R
boxplot(boreal_data$Fdis, horizontal = TRUE) # no majour outliers
plot(density(boreal_data$Fdis)) # use gaussian distribution in modelling 

# Collinearity of X was assessed in 6_biomass_models.R, these were the correlated
# variables:
# sst_mean - sbt_mean -> remove sst_mean
# sbt_mean - depth    -> remove depth
# sic variables can be removed because they are almost all 0s

## Boreal model  with fishing - selection, validation, prediction ----

# Selection
boreal_model <- model_selection(directory = "models/Fdis/boreal/with_fishing",
                                name = "boreal",
                                data = boreal_data,
                                mesh_cutoff = 65, 
                                family = gaussian(link = "identity"),
                                fixed_formula = Fdis ~ 
                                  sbt_mean +
                                  sbt_sd +
                                  #sst_mean +
                                  sst_sd +
                                  #sic_mean + 
                                  #sic_sd +
                                  log_chla_mean +
                                  #depth +
                                  log_fishing_effort)

# Validation
model_validation(model = boreal_model,
                 data = boreal_data,
                 directory = "models/Fdis/boreal/with_fishing")

# Prediction
model_prediction(model = boreal_model,
                 directory = "models/Fdis/boreal/with_fishing")

# Cleaning 
rm(list = setdiff(ls(), c("final_data",
                          "model_prediction",
                          "model_selection",
                          "model_validation")))

## End