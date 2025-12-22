## File: 6_biomass_models.R
## Purpose: model the response of biomass to diversity, warming, and fishing
## Author: Filippomaria Cassarino
## Date: 
## ------------------------------------------------------------------------ ----
## Notes ----

# Fix Fric outlier (keep or remove? now removed).
# check sdmTMB::simulate (not all converge)

## Models checks (good or bad):

# general good
# fishing good (positive trend in res vs fitted)

# Arctic good (non-linear trend in res vs fitted)
# Boreal good (insanity of the spatial model, which would not be selected)

# 2004-2011 good (minor patterns in rs vs fitted)
# 2012-2016 good
# 2017-2022 good (positive trend in res vs fitted, dispersion issues)

# remove_1 good
# remove_2 good
# remove_3 good
# remove_4 good
# remove_5 good (insanity of the spatial model, which would not be selected)
# remove_6 good (insanity of the spatial model, which would not be selected)
# remove_7 good (insanity of the spatial model, which would not be selected)
# remove_8 good (insanity of the spatial model, which would not be selected)
# remove_9 good (insanity of the spatial model, which would not be selected)
# remove_10 good

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
                      "sdmTMB", 
                      "DHARMa") 

install_load_function(required_pakages)

# Objects to keep when cleaning
keep <- c("install_load_function",
          "manual_model_selection",
          "model_selection_only_random",
          "model_selection",
          "model_prediction",
          "model_validation",
          "start_data", 
          "vif",
          "general_data",
          "keep")

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## General model ----

# This is the model for biomass including the most data (all areas and time, no
# fishing effort) and selected only for the random term to 
# enable direct comparisons with other models

# Prepare data of the whole study area (all) for exploration
general_data <- start_data %>%
  select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

# Variables distribution and outliers
vars <- setdiff(colnames(general_data), c("year", "latitude",
                                          "longitude", "haul_id"))

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(general_data[[vars[i]]], # biomass, chla and fishing are skewed
          horizontal = TRUE,
          main = vars[i])  
}
par(mfrow = c(1, 1))  

# Transform 
boxplot(log(general_data$total_biomass), horizontal = TRUE) # no majour outliers
plot(density(log(general_data$total_biomass))) # use log-normal distribution in modelling

boxplot(general_data$log_chla_mean, horizontal = TRUE) # better

# Collinearity X
cor_m_pearson <- general_data %>%       # Pearson correlation
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

general_data %>% select(
  Fric,
  Feve,
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth) %>% vif() # Teve is barely above the treshold, since it is a central
                   # variable, it is maintained

# Selection
general_model <- model_selection_only_random(
  directory = "models/biomass/general_model",
  name = "general_biomass",
  data = general_data,  
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
model_validation(model = general_model,
                 data = general_data,
                 name = "general_biomass",
                 directory = "models/biomass/general_model")

# Prediction
#model_prediction(model = general_model,
#                 name = "general_biomass",
#                 directory = "models/biomass/general_model")

# Cleaning
rm(list = setdiff(ls(), keep))

#
## Fishing model ----

# This model investigates the effect of fishing, incorporating all areas and 
# time. It is selected for the fixed term to compare with the general model.

# Prepare data for exploration
fishing_data <- start_data %>%
  select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           fishing_effort, log_fishing_effort,     # fishing
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

# Variables distribution and outliers
vars <- setdiff(colnames(fishing_data), c("year", "latitude",
                                        "longitude", "haul_id"))

par(mfrow = c(4, 5), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(fishing_data[[vars[i]]], # chla and fishing should be log-transformed
          horizontal = TRUE,
          main = vars[i])  
}
par(mfrow = c(1, 1))  

# Collinearity X
cor_m_pearson <- fishing_data %>%       # Pearson correlation
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

fishing_data %>% select(
  Fric,
  Feve,
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth,
  log_fishing_effort) %>% vif() # sic_mean is barely above the treshold, 
                                # it is maintained to compare with general

# Selection
fishing_model <- model_selection_only_random(
  directory = "models/biomass/fishing_model",
  name = "fishing_biomass",
  data = fishing_data,
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
    depth +
    log_fishing_effort)

# Validation
model_validation(model = fishing_model,
                 data = fishing_data,
                 name = "fishing_biomass",
                 directory = "models/biomass/fishing_model")

# Prediction
# model_prediction(model = fishing_model, name = "fishing_biomass", directory = "models/biomass/fishing_model")

# Cleaning
rm(list = setdiff(ls(), keep))

#
## Arctic and Boreal models ----

# These models explore differences between Arctic and Boreal areas, here divided
# using sst_mean. # They are selected only for the random term.

# Identify Arctic and Boreal data
divide_arc <- quantile(start_data$sst_mean, 0.25, na.rm = TRUE)
divide_bor <- quantile(start_data$sst_mean, 0.75, na.rm = TRUE)

# Prepare data for exploration
arctic_data <- start_data %>%
  filter(sst_mean <= divide_arc) %>%
  select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

boreal_data <- start_data %>%
  filter(sst_mean >= divide_bor) %>%
  select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

# Variables distribution and outliers
vars <- setdiff(colnames(arctic_data), c("year", "latitude",
                                         "longitude", "haul_id"))

# Arctic data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(arctic_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i]) # chla_mean and biomass are skewed
}
par(mfrow = c(1, 1))  

# Boreal data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))  
for(i in seq_along(vars)) {
  boxplot(boreal_data[[vars[i]]],
          horizontal = TRUE,
          main = vars[i]) # chla_mean is not so skewed, but will use log to keep comparison
}
par(mfrow = c(1, 1)) 

# Collinearity X in Arctic
cor_m_pearson <- arctic_data %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson

# Correlated variables:
# Fric - Fdis         -> remove Fdis
# Feve - Teve         -> should remove Feve
# Teve - Tdiv         -> remove Tdiv
# sst_mean - sic_mean -> remove sst_mean
# latitude must be retained to account for spatial dependencies

arctic_data %>% select(
  Fric,
  #Feve, 
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd, # remove to improve comparison
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth) %>% vif() # with Feve removed, all VIF < 2

# Collinearity X in Boreal data
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
# sbt_mean - depth    -> should remove depth
# latitude must be retained to account for spatial dependencies

boreal_data %>% select(
  Fric,
  Feve,
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  #sic_mean,
  #sic_sd,
  log_chla_mean,
  #depth
  ) %>% vif() # removing depth, all VIF < 2

# Arctic model selection
arctic_model <- model_selection_only_random(
  directory = "models/biomass/arctic_boreal_model",
  name = "arctic_biomass",
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
    depth) # removing depth does not change much compared to the boreal model

# Boreal model selection 
boreal_model <- model_selection_only_random(
  directory = "models/biomass/arctic_boreal_model",
  name = "boreal_biomass",
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
    sic_mean + 
    #sic_sd +   # removed to improve comparisons 
    #depth +   # removing depth has a strong effect on sbt_mean
    log_chla_mean)

# Arctic and Boreal models validation
model_validation(model = arctic_model,
                 data = arctic_data,
                 name = "arctic_biomass",
                 directory = "models/biomass/arctic_boreal_model")

model_validation(model = boreal_model,
                 data = boreal_data,
                 name = "boreal_biomass",
                 directory = "models/biomass/arctic_boreal_model")

# Arctic and Boreal models prediction
#model_prediction(model = arctic_model, name = "arctic_biomass",
#                 directory = "models/biomass/arctic_boreal_model")

#model_prediction(model = boreal_model, name = "boreal_biomass",
#                 directory = "models/biomass/arctic_boreal_model")

# Cleaning
rm(list = setdiff(ls(), keep))

#
## Partial timeline models ----

# These models investigate differences between the three time periods
# identified by Eriksen et al., 2025. They are selected only for the
# random term.

# Divide data based on Eriksen et al., 2025
data_2004_2011 <- filter(start_data, year < 2012) %>% # coldest period
 select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

data_2012_2016 <- filter(start_data, year > 2011 & year < 2017) %>% # warmer period
 select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

data_2017_2022 <- filter(start_data, year > 2016) %>% # cooler period
  select(c(total_biomass,          # response
           Fric, Feve, Fdis,       # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(-c(haul_id, total_biomass, year, latitude, longitude),
                base::scale))

# Collinearity X 
data_2004_2011 %>% select(
  Fric,
  Feve,
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth) %>% vif() # Teve is slightly above threshold, but it is kept

# Collinearity X 
data_2012_2016 %>% select(
  Fric,
  Feve,
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  sic_mean,
  #sic_sd,
  log_chla_mean,
  depth) %>% vif() # All < 2

# Collinearity X (more problematic)

# Collinearity X 
cor_m_pearson <- data_2017_2022 %>%       # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson") 

cor_m_pearson[lower.tri(cor_m_pearson, diag = TRUE)] <- NA

cor_df_pearson <- as.data.frame(as.table(cor_m_pearson)) %>%
  filter(!is.na(Freq), abs(Freq) > 0.65) 

cor_df_pearson # sbt_mean and sic_mean are collinear

data_2017_2022 %>% select(
  Fric,
  Feve,
  #Fdis,
  Tric,
  Teve,
  #Tdiv, 
  sbt_mean, 
  sbt_sd,
  #sst_mean, 
  sst_sd,
  #sic_mean,
  #sic_sd,
  log_chla_mean,
  depth) %>% vif() # Teve above, its effect on others will be checked

# Select 2004_2011 model
model_2004_2011 <- model_selection_only_random(
  directory = "models/biomass/time_period_model",
  name = "2004_2011_biomass",
  data = data_2004_2011,
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

# Select 2012_2016 model
model_2012_2016 <- model_selection_only_random(
  directory = "models/biomass/time_period_model",
  name = "2012_2016_biomass",
  data = data_2012_2016,
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

# Select 2017_2022 model
model_2017_2022 <- model_selection_only_random(
  directory = "models/biomass/time_period_model",
  name = "2017_2022_biomass",
  data = data_2017_2022,
  mesh_cutoff = 65, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    Fric +
    Feve +
    #Fdis +
    Tric +
    Teve +# removing it makes the effect of Feve and Fric mildly negative. Keep
    #Tdiv +
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    #sic_mean +  # removed to avoid collinearity
    #sic_sd +       
    log_chla_mean +
    depth)

# Models' validation
model_validation(model = model_2004_2011,
                 data = data_2004_2011,
                 name = "2004_2011_biomass",
                 directory = "models/biomass/time_period_model")

model_validation(model = model_2012_2016,
                 data = data_2012_2016,
                 name = "2012_2016_biomass",
                 directory = "models/biomass/time_period_model")

model_validation(model = model_2017_2022,
                 data = data_2017_2022,
                 name = "2017_2022_biomass",
                 directory = "models/biomass/time_period_model")

# Model's prediction
#model_prediction(model = model_2004_2011, name = "2004_2011_biomass",
#                 directory = "models/biomass/time_period_model")

#model_prediction(model = model_2012_2016, name = "2012_2016_biomass",
#                 directory = "models/biomass/time_period_model")

# Model's prediction
#model_prediction(model = model_2017_2022, name = "2017_2022_biomass",
#                 directory = "models/biomass/time_period_model")

# Cleaning
rm(list = setdiff(ls(), keep))

#
## Partial biomass models ----

# These models investigate the effect of progressively removing dominant species
# from the total_biomass. They ate selected only for the random term. 
# Note that data exploration is the same as general_model

# Prepare data
b_data <- start_data

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
           Fric, Feve, Fdis, # FD
           Tric, Teve, Tdiv,       # TD
           sst_mean, sst_sd,       # surface temperature
           sbt_mean, sbt_sd,       # bottom temperature
           sic_mean, sic_sd,       # sea ice
           chla_mean, log_chla_mean,# productivity
           year, depth, longitude, latitude, haul_id)) %>%
  drop_na() %>%
  mutate(across(starts_with("remove_"), ~ .x + 0.0001), # a few observations are slightly negative due to rounding
         across(-c(haul_id, total_biomass, all_of(new_res), year, latitude, longitude),
                base::scale))

# Models
for (i in 1:length(new_res)) {
  
  name <- paste0(new_res[i], "_biomass")
  
  # Fitting
  model <- model_selection_only_random(
    directory = "models/biomass/species_removal_model",
    name = name,
    data = biomass_data,
    mesh_cutoff = 65, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i],
      "  ~  
      Fric +
      Feve +
      Tric +
      Teve + 
      sbt_mean +
      sbt_sd +
      sst_sd +
      sic_mean +        
      log_chla_mean +
      depth"))) # Fdis, Tdiv, sst_mean, sic_sd are removed
  
  # Validation
  model_validation(model = model,
                   data = biomass_data,
                   name = paste0(new_res[i], "_biomass"),
                   directory = "models/biomass/species_removal_model")
  
}

# Cleaning
rm(list = setdiff(ls(), keep))

#
## ------------------------------------------------------------------------ ----


