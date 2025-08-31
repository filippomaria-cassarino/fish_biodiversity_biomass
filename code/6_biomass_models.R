## File: 6_biomass_models.R
## Purpose: model the response of biomass to diversity, warming, and fishing
## Author: Filippomaria Cassarino
## Date: 
## ------------------------------------------------------------------------ ----
## Notes ----

# All model is good but has some outlier problems, remove outlier for Fric maybe

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

## ------------------------------------------------------------------------ ----
## All model no fishing - data exploration ----

# Prepare data of the whole study area (all) for exploration
check_data <- final_data %>%
  select(c(total_biomass,          # response
           Fric, Feve, Fdis, # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() 

# Variables distribution and outliers
vars <- setdiff(colnames(check_data), c("year", "latitude",
                                           "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(check_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i])  # chla_mean and biomass are very skewed
}
par(mfrow = c(1, 1))  

# Transform 
boxplot(log(check_data$total_biomass), horizontal = TRUE) # no majour outliers
plot(density(log(check_data$total_biomass))) # use log-normal distribution in modelling

boxplot(log(check_data$chla_mean), horizontal = TRUE) # log_chla_mean is better


# Collinearity X
cor_m_pearson <- check_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis     -> remove Fdis
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sbt_mean -> remove sst_mean
# sic_mean - sic_sd   -> remove sic_sd
# latitude must be retained to account for spatial dependencies

# Data for modelling
all_data <- check_data  %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale)) 

## All model no fishing - selection, validation, prediction ----

# Selection
all_model <- model_selection(directory = "models/biomass/all/no_fishing",
                             name = "all",
                             data = all_data,
                             mesh_cutoff = 65, 
                             family = lognormal(link = "log"),
                             fixed_formula = total_biomass ~  
                               Fric +
                               Feve +
                               #Fdis +
                               Tric +
                               Teve + 
                               #Tdiv +
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
                 directory = "models/biomass/all/no_fishing")

# Prediction
model_prediction(model = all_model,
                 directory = "models/biomass/all/no_fishing")

## Arctic model no fishing - data exploration ----

# Identify arctic data
divide_arc <- quantile(all_data$sst_mean, 0.25, na.rm = TRUE)

# Prepare data for exploration
arctic_data <- all_data[all_data$sst_mean <= divide_arc, ]

# Variables distribution and outliers
vars <- setdiff(colnames(arctic_data), c("year", "latitude",
                                        "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(arctic_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i]) # chla_mean and biomass are very skewed, do same as all_data
}
par(mfrow = c(1, 1))  

# Collinearity X
cor_m_pearson <- arctic_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis         -> remove Fdis
# Feve - Teve         -> remove Feve
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sic_mean -> remove sst_mean
# latitude must be retained to account for spatial dependencies

## Arctic model no fishing - selection, validation, prediction ----

# Selection
arctic_model <- model_selection(directory = "models/biomass/arctic/no_fishing",
                                name = "arctic",
                                data = arctic_data,
                                mesh_cutoff = 65, 
                                family = lognormal(link = "log"),
                                fixed_formula = total_biomass ~  
                                  Fric +
                                  #Feve +
                                  #Fdis +
                                  Tric +
                                  Teve + 
                                  #Tdiv +
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
                 directory = "models/biomass/arctic/no_fishing")

# Prediction
model_prediction(model = arctic_model,
                 directory = "models/biomass/arctic/no_fishing")

## Boreal model no fishing - data exploration ----

# Identify boreal data
divide_bor <- quantile(all_data$sst_mean, 0.75, na.rm = TRUE)

# Prepare data
boreal_data <- all_data[all_data$sst_mean >= divide_bor, ]

# Variables distribution and outliers
vars <- setdiff(colnames(boreal_data), c("year", "latitude",
                                         "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(boreal_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i]) # chla_mean and biomass are very skewed, do same as all_data
}
par(mfrow = c(1, 1))  

# Collinearity X
cor_m_pearson <- boreal_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis     -> remove Fdis
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sbt_mean -> remove sst_mean
# sbt_mean - depth    -> remove depth
# sic variables can be removed because they are almost all 0s
# latitude must be retained to account for spatial dependencies

## Boreal model no fishing  - selection, validation, prediction ----

# Selection
boreal_model <- model_selection(directory = "models/biomass/boreal/no_fishing",
                                name = "boreal",
                                data = boreal_data,
                                mesh_cutoff = 65, 
                                family = lognormal(link = "log"),
                                fixed_formula = total_biomass ~  
                                  Fric +
                                  Feve +
                                  #Fdis +
                                  Tric +
                                  Teve + 
                                  #Tdiv +
                                  sbt_mean +
                                  sbt_sd +
                                  #sst_mean +
                                  sst_sd +
                                  #sic_mean + 
                                  #sic_sd +
                                  log_chla_mean #+
                                  #depth
                                )

# Validation
model_validation(model = boreal_model,
                 data = boreal_data,
                 directory = "models/biomass/boreal/no_fishing")

# Prediction
model_prediction(model = boreal_model,
                 directory = "models/biomass/boreal/no_fishing")


# Cleaning before the model with fishing
rm(list = setdiff(ls(), c("final_data",
                          "model_prediction",
                          "model_selection",
                          "model_validation")))

## ------------------------------------------------------------------------ ----
## All model with fishing - data exploration ----

# Prepare data of the whole study area (all) for exploration
check_data <- final_data %>%
  select(c(total_biomass,          # response
           Fric, Feve, Fdis, # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           log_chla_mean,          # productivity
           log_fishing_effort,     # fishing
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() 

# Variables distribution and outliers
vars <- setdiff(colnames(check_data), c("year", "latitude",
                                        "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(check_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i])  # chla_mean, biomass and fishing are very skewed
}
par(mfrow = c(1, 1))  

# Transform 
boxplot(log(check_data$total_biomass), horizontal = TRUE) # no majour outliers
plot(density(log(check_data$total_biomass))) # use log-normal distribution in modelling

boxplot(check_data$log_chla_mean, horizontal = TRUE) # better
boxplot(check_data$log_fishing_effort, horizontal = TRUE) # better

# Collinearity X
cor_m_pearson <- check_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis         -> remove Fdis
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sbt_mean -> remove sst_mean
# sic_mean - sic_sd   -> remove sic_sd
# latitude must be retained to account for spatial dependencies

# Data for modelling
all_data <- check_data  %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

## All model with fishing - selection, validation, prediction ----

# Selection
all_model <- model_selection(directory = "models/biomass/all/with_fishing",
                             name = "all",
                             data = all_data,
                             mesh_cutoff = 65, 
                             family = lognormal(link = "log"),
                             fixed_formula = total_biomass ~  
                               Fric +
                               Feve +
                               #Fdis +
                               Tric +
                               Teve + 
                               #Tdiv +
                               sbt_mean * log_fishing_effort +
                               sbt_sd * log_fishing_effort +
                               #sst_mean +
                               sst_sd +
                               sic_mean + 
                               #sic_sd +
                               log_chla_mean +
                               depth)

# Validation
model_validation(model = all_model,
                 data = all_data,
                 directory = "models/biomass/all/with_fishing")

# Prediction
model_prediction(model = all_model,
                 directory = "models/biomass/all/with_fishing")

## Arctic model with fishing - data exploration ----

# Identify arctic data
divide_arc <- quantile(all_data$sst_mean, 0.25, na.rm = TRUE)

# Prepare data for exploration
arctic_data <- all_data[all_data$sst_mean <= divide_arc, ]

# Variables distribution and outliers
vars <- setdiff(colnames(arctic_data), c("year", "latitude",
                                         "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(arctic_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i]) # chla_mean and biomass are very skewed, do same as all_data
}
par(mfrow = c(1, 1))  

# Collinearity X
cor_m_pearson <- arctic_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis         -> remove Fdis
# Feve - Teve         -> remove Feve
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sst_sd   -> remove sst_mean
# sic_mean - sic_sd   -> remove sic_sd
# latitude must be retained to account for spatial dependencies

## Arctic model with fishing - selection, validation, prediction ----

# Selection
arctic_model <- model_selection(directory = "models/biomass/arctic/with_fishing",
                                name = "arctic",
                                data = arctic_data,
                                mesh_cutoff = 65, 
                                family = lognormal(link = "log"),
                                fixed_formula = total_biomass ~  
                                  Fric +
                                  #Feve +
                                  #Fdis +
                                  Tric +
                                  Teve + 
                                  #Tdiv +
                                  sbt_mean +
                                  sbt_sd +
                                  #sst_mean +
                                  sst_sd +
                                  sic_mean + 
                                  #sic_sd +
                                  log_chla_mean +
                                  depth)

# Validation
model_validation(model = arctic_model,
                 data = arctic_data,
                 directory = "models/biomass/arctic/with_fishing")

# Prediction
model_prediction(model = arctic_model,
                 directory = "models/biomass/arctic/with_fishing")

## Boreal model with fishing - data exploration ----

# Identify boreal data
divide_bor <- quantile(all_data$sst_mean, 0.75, na.rm = TRUE)

# Prepare data
boreal_data <- all_data[all_data$sst_mean >= divide_bor, ]

# Variables distribution and outliers
vars <- setdiff(colnames(boreal_data), c("year", "latitude",
                                         "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(boreal_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i]) # chla_mean and biomass are very skewed, do same as all_data
}
par(mfrow = c(1, 1))  

# Collinearity X
cor_m_pearson <- boreal_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis         -> remove Fdis
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sbt_mean -> remove sst_mean
# sbt_mean - depth    -> remove depth
# sic variables can be removed because they are almost all 0s
# latitude must be retained to account for spatial dependencies

## Boreal model  with fishing - selection, validation, prediction ----

# Selection
boreal_model <- model_selection(directory = "models/biomass/boreal/with_fishing",
                                name = "boreal",
                                data = boreal_data,
                                mesh_cutoff = 65, 
                                family = lognormal(link = "log"),
                                fixed_formula = total_biomass ~  
                                  Fric +
                                  Feve +
                                  #Fdis +
                                  Tric +
                                  Teve + 
                                  #Tdiv +
                                  sbt_mean +
                                  sbt_sd +
                                  #sst_mean +
                                  sst_sd +
                                  #sic_mean + 
                                  #sic_sd +
                                  log_chla_mean #+
                                  #depth
)

# Validation
model_validation(model = boreal_model,
                 data = boreal_data,
                 directory = "models/biomass/boreal/with_fishing")

# Prediction
model_prediction(model = boreal_model,
                 directory = "models/biomass/boreal/with_fishing")

# Cleaning 
rm(list = setdiff(ls(), c("final_data",
                          "model_prediction",
                          "model_selection",
                          "model_validation")))

## End