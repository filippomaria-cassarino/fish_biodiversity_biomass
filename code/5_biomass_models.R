## File: 5_biomass_models.R
## Purpose: model the response of biomass to diversity, warming, and fishing
## Author: Filippomaria Cassarino
## Date: 2 Apr 2026
## Notes ----

# All the correlation tables must be updated once the final set of
# variables is selected

# It would be nice to merge all .txt and .pdf at the end of the script 
# to have them as attachments. It could also be just a pdf file

## Model validation
# general_biomass OK
# fishing_biomass NOT OK

## Library ---- 

# Import data, objects and functions
final_data <- read.csv("data/final/final_data.csv")

# Load functions
load("tools/model_selection_functions.RData")
load("tools/model_prediction_function.RData")
load("tools/vif_function.RData")
load("tools/install_load_packages_function.RData")

# Load required packages
required_pakages <- c("dplyr",
                      "ggplot2",
                      "tidyr",   
                      "sdmTMB", 
                      "DHARMa",
                      "pdftools",
                      "readr") 

install_load_packages(required_pakages)

# Recurring variables
recurring_vars <- c(
  "total_biomass",
  "kde_fric",
  "kde_feve",       
  "tric",
  "teve",       
  "sst_mean",       
  "sbt_mean",       
  "sic_mean", 
  "sbt_sd",
  "sst_sd",
  "sic_sd", 
  "log_chla_mean",
  "depth",
  "year",
  "longitude",
  "latitude",
  "haul_id"
)

# Variables to inspect during data exploration
check_vars <- setdiff(recurring_vars,
                c("year", "latitude", "longitude", "haul_id"))

# Objects to keep when cleaning
keep <- c(
  "keep",
  "final_data",
  "vif",
  "model_selection",
  "manual_model_selection", 
  "model_prediction",
  "recurring_vars",
  "check_vars",
  "start",
  "end"
)

rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----

## General model ----

# This is the model for biomass including the most data including 
# all areas and time periods, but excluding fishing effort 

# Prepare data of the whole study area (all) for exploration
general_data <- final_data %>%
  select(all_of(recurring_vars)) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.)))) # mean = 0; sd = 1 

# Variables distribution and outliers
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))

for(i in seq_along(check_vars)) {
  boxplot(general_data[[check_vars[i]]], # biomass, chla and fishing are skewed
          horizontal = TRUE,
          main = check_vars[i])  
}

par(mfrow = c(1, 1))  

# Transform 
boxplot(log(general_data$total_biomass), horizontal = TRUE) # no major outliers
plot(density(log(general_data$total_biomass))) # use log-normal distribution 

# Collinearity X
general_data %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.8112112 sbt_mean
#2 sst_mean sic_mean -0.8251594 sic_mean
#3 sst_mean   sic_sd -0.8267052
#4 sic_mean   sic_sd  0.8740684 sic_mean

general_data %>% select(
  kde_fric,
  kde_feve,       
  tric,
  teve,           
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
  ) %>% vif() # All VIF values < 2

# Selection
general_model <- model_selection(
  directory = "models/biomass",
  name = "general_biomass",
  data = general_data,  
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +      
    tric +
    teve +       
    #sst_mean +      
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    depth
  )

# Prediction
model_prediction(model = general_model,
                 vars = c("teve", "sbt_mean"),
                 name = "general_biomass",
                 directory = "figures/model_figures")

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----

## Fishing model ----

# This model investigates the effect of fishing, incorporating all areas and 
# time periods.

# Prepare data for exploration
fishing_data  <- final_data %>%
  select(
    log_sum_fishing,
    all_of(recurring_vars)
  ) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.)))) 

# Variables distribution and outliers
par(mfrow = c(4, 5), mar = c(2, 4, 2, 1)) 

for(i in seq_along(check_vars)) {
  boxplot(fishing_data[[check_vars[i]]], 
          horizontal = TRUE,
          main = check_vars[i])  
}

par(mfrow = c(1, 1))  

boxplot(fishing_data$log_sum_fishing, 
        horizontal = TRUE, main = "log_sum_fishing")

# Collinearity X
fishing_data %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.8566891 sbt_mean
#2 sst_mean sic_mean -0.8188343 sic_mean
#3 sst_mean   sic_sd -0.8322784
#4 sic_mean   sic_sd  0.9059084 sic_mean

fishing_data %>% select(
  log_sum_fishing,
  kde_fric,
  kde_feve,      
  tric,
  teve,       
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
  ) %>% vif() # All VIF values < 2

# Selection
fishing_model <- model_selection(
  directory = "models/biomass",
  name = "fishing_biomass",
  data = fishing_data,
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~ 
    log_sum_fishing +
    kde_fric +
    kde_feve +     
    tric +
    teve +
    #sst_mean +      
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    depth
  )

# The validation has some warnings, but when simplifying the model, the 
# coefficients remain stable, so the model should be trusted.

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----

## Arctic model ----

# This model focuses on Arctic areas

# Identify Arctic data
divide_arc <- quantile(final_data$sbt_mean, 0.25, na.rm = TRUE)

# Prepare Arctic data for exploration
arctic_data <- final_data %>%
  filter(sbt_mean <= divide_arc) %>%
  select(all_of(recurring_vars)) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# Arctic data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1)) 

for(i in seq_along(check_vars)) {
  boxplot(arctic_data[[check_vars[i]]],
          horizontal = TRUE,
          main = check_vars[i]) 
}

par(mfrow = c(1, 1))  

# Collinearity X in Arctic data
arctic_data %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1      Var2       Freq       Choice
#1 sst_mean  sic_mean -0.8677860 sic_mean
#2 sst_mean    sic_sd -0.8081397 
#3 sic_mean    sic_sd  0.8431559 sic_mean

arctic_data %>% select( # VIF
  kde_fric,
  kde_feve,       
  tric,
  teve, 
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # All VIF values < 2

# Arctic model selection
arctic_model <- model_selection(
  directory = "models/biomass",
  name = "arctic_biomass",
  data = arctic_data,
  mesh_cutoff = 60,  
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +      
    tric +
    teve +     
    #sst_mean +      
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    depth
)

## ------------------------------------------------------------------------ ----

## Boreal model ----

# This model focuses on Boreal areas

# Identify Boreal data
divide_bor <- quantile(final_data$sbt_mean, 0.75, na.rm = TRUE)

# Prepare Boreal data for exploration
boreal_data <- final_data %>%
  filter(sbt_mean >= divide_bor) %>%
  select(all_of(recurring_vars)) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# Boreal data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1)) 

for(i in seq_along(check_vars)) {
  boxplot(boreal_data[[check_vars[i]]],
          horizontal = TRUE,
          main = check_vars[i]) 
}

par(mfrow = c(1, 1)) 

# Collinearity X in Boreal data
boreal_data %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.7911535 sbt_mean
#2 sic_mean   sic_sd  0.9070993 neither, very few observations

boreal_data %>% select( # VIF
  kde_fric,
  kde_feve,       
  tric,
  teve,      
  #sst_mean,       
  sbt_mean,       
  #sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # all VIF < 2

# Boreal model selection 
boreal_model <- model_selection(
  directory = "models/biomass",
  name = "boreal_biomass",
  data = boreal_data,
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +      
    tric +
    teve +      
    #sst_mean +      
    sbt_mean +       
    #sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    depth
)

# Selection failed for spatial_random_fields_with_anisotropy,
# but it would probably not be selected as the best structure given that 
# spatiotemporal random fields are favoured in all the other models

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----

## 2004_2011 model ----

# This model focuses on the first (coldest) time period identifies by
# Eriksen et al., 2025. https://doi.org/10.1038/s41598-025-96964-x

# Divide data based on Eriksen et al., 2025
data_2004_2011 <- filter(final_data, year < 2012) %>% # coldest period
  select(all_of(recurring_vars)) %>%
  drop_na()  %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# 2004_2011 data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1)) 

for(i in seq_along(check_vars)) {
  
  boxplot(data_2004_2011[[check_vars[i]]],
          horizontal = TRUE,
          main = check_vars[i]) 
}

par(mfrow = c(1, 1)) 

# Collinearity X in 2004-2011 data
data_2004_2011 %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.7685743 sbt_mean
#2 sst_mean sic_mean -0.8326126 sic_mean
#3 sst_mean   sic_sd -0.8256690 
#4 sic_mean   sic_sd  0.8605019 sic_mean

data_2004_2011 %>% select( # VIF
  kde_fric,
  kde_feve,       
  tric,
  teve,       
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # sic_mean is just slightly above the treshold, so it is maintained

# Select 2004_2011 model
model_2004_2011 <- model_selection(
  directory = "models/biomass",
  name = "2004_2011_biomass",
  data = data_2004_2011,
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +
    tric +
    teve +
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +       
    log_chla_mean +
    depth
  )

## ------------------------------------------------------------------------ ----

## 2012_2016 model ----

# This model focuses on the second (warmest) time period identifies by
# Eriksen et al., 2025. https://doi.org/10.1038/s41598-025-96964-x

data_2012_2016 <- filter(final_data, year > 2011 & year < 2017) %>% # warmer period
  select(all_of(recurring_vars)) %>%
  drop_na()  %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# 2012_2016 data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1)) 

for(i in seq_along(check_vars)) {
  
  boxplot(data_2012_2016[[check_vars[i]]],
          horizontal = TRUE,
          main = check_vars[i]) 
}

par(mfrow = c(1, 1)) 

# Collinearity X in 2004-2011 data
data_2012_2016 %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.8349426 sbt_mean
#2 sst_mean sic_mean -0.7548561 sic_mean
#3 sst_mean   sic_sd -0.7666632 
#4 sic_mean   sic_sd  0.8682476 sic_mean

data_2012_2016 %>% select( # VIF
  kde_fric,
  kde_feve,       
  tric,
  teve,     
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # All VIF values < 2

# Select 2012_2016 model
model_2012_2016 <- model_selection(
  directory = "models/biomass",
  name = "2012_2016_biomass",
  data = data_2012_2016,
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +
    tric +
    teve + 
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +       
    log_chla_mean +
    depth
)

## ------------------------------------------------------------------------ ----

## 2017_2022 model ----

# This model focuses on the third (cooler) time period identifies by
# Eriksen et al., 2025. https://doi.org/10.1038/s41598-025-96964-x

data_2017_2022 <- filter(final_data, year > 2016) %>% # cooler period
  select(all_of(recurring_vars)) %>%
  drop_na()  %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

# 2017_2022 data distribution
par(mfrow = c(4, 4), mar = c(2, 4, 2, 1)) 

for(i in seq_along(check_vars)) {
  
  boxplot(data_2017_2022[[check_vars[i]]],
          horizontal = TRUE,
          main = check_vars[i]) 
}

par(mfrow = c(1, 1)) 

# Collinearity X in 2004-2011 data
data_2017_2022 %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.8728749 sbt_mean
#2 sst_mean sic_mean -0.8686514 
#3 sst_mean   sic_sd -0.8842051
#4 sbt_mean   sic_sd -0.6711975 sbt_mean
#5 sic_mean   sic_sd  0.9257741 sic_mean

data_2017_2022 %>% select( # VIF
  kde_fric,
  kde_feve,      
  tric,
  teve,       
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # sic_mean is above two, but not my much. Keep to favour comparisons

# Select 2017_2022 model
model_2017_2022 <- model_selection(
  directory = "models/biomass",
  name = "2017_2022_biomass",
  data = data_2017_2022,
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +
    tric +
    teve + 
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

## Partial biomass models ----

# These models investigate the effect of the progressiveremoval of 
# dominant species from the total_biomass. 
# Note that data exploration is the same as general_model

# Prepare data
b_data <- final_data

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
  select(  
    all_of(new_res),         
    all_of(recurring_vars)
  ) %>%
  drop_na() %>%
  mutate(across(
    starts_with("remove_"),
    ~ .x + 0.0001), # few observations are slightly negative due to rounding
         across(
    -c(
      starts_with("remove_"),
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.)))) 

# Responses distribution
par(mfrow = c(3, 4), mar = c(2, 4, 2, 1)) 

for(i in seq_along(new_res)) {
  
  plot(density(log(biomass_data[[new_res[i]]])),
       main = new_res[i])
}

par(mfrow = c(1, 1)) 

# Collinearity is the same as general model

# Models
for (i in 1:length(new_res)) {
  
  name <- paste0(new_res[i], "_biomass")
  
  # Fitting
  model <- model_selection(
    directory = "models/biomass",
    name = name,
    data = biomass_data,
    mesh_cutoff = 60, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i],
      "  ~  
      kde_fric +
      kde_feve +
      tric +
      teve + 
      sbt_mean +
      sbt_sd +
      sst_sd +
      sic_mean +        
      log_chla_mean +
      depth"
      )))
  
  # Return remove_10_model to visualize predictions
  if(formula(model)[[2]] == "remove_10") {
    
    remove_10_model <- model
    
    return(remove_10_model)
  }
  
  message(paste0(name, " model fitted and inspected"))
  
}

# Check remove_10 model predictions for tric and teve
model_prediction(model = remove_10_model,
                 name = "remove_10_biomass_model",
                 vars = c("tric", "teve"),
                 directory = "figures/model_figures")

# Merge all pdfs
pdf_combine(
  input = paste0("models/biomass/", new_res, "_biomass_model_selection.pdf"),
  output = "models/biomass/partial_biomass_models_selection.pdf"
)

# Cleaning
file.remove(
  paste0("models/biomass/", new_res, "_biomass_model_selection.pdf")
)

rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----

## General model with Pielou's evenness ----

# This model uses Pielou's taxonomic evenness rather than Simpson's to 
# investigate  whole-community effects

# Prepare data of the whole study area (all) for exploration
pielou_data <- final_data %>%
  select(all_of(recurring_vars), p_teve) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.)))) # mean = 0; sd = 1 

# New variable distribution and outliers
  boxplot(pielou_data$p_teve,
          horizontal = TRUE,
          main = "Pielou's evenness")  

# Collinearity X
pielou_data %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq       Choice
#1 sst_mean sbt_mean  0.8112112 sbt_mean
#2 sst_mean sic_mean -0.8251594 
#3 sst_mean   sic_sd -0.8267052 
#4 sic_mean   sic_sd  0.8740684 sic_mean

pielou_data %>% select(
  kde_fric,
  kde_feve,      
  tric,
  p_teve,           
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # All VIF values < 2

# Selection
pielou_model <- model_selection(
  directory = "models/biomass",
  name = "pielou_biomass",
  data = pielou_data,  
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +     
    tric +
    p_teve +       
    #sst_mean +      
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----

## Partial biomass models with Pielou's evenness ----

# These models use Pielou's taxonomic evenness rather than Simpson's to 
# investigate  whole-community effects while removing top species

b_data <- final_data

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
  select(  
    all_of(new_res),         
    all_of(recurring_vars),
    p_teve
  ) %>%
  drop_na() %>%
  mutate(across(
    starts_with("remove_"),
    ~ .x + 0.0001), # few observations are slightly negative due to rounding
    across(
      -c(
        starts_with("remove_"),
        haul_id,
        total_biomass,
        year, 
        latitude,
        longitude
      ), ~ as.numeric(scale(.)))) 

# Collinearity is the same as pielou model

# Models
for (i in 1:length(new_res)) {
  
  name <- paste0("pielou_", new_res[i], "_biomass")
  
  # Fitting
  model <- model_selection(
    directory = "models/biomass",
    name = name,
    data = biomass_data,
    mesh_cutoff = 60, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i],
      "  ~  
      kde_fric +
      kde_feve +
      tric +
      p_teve + 
      sbt_mean +
      sbt_sd +
      sst_sd +
      sic_mean +        
      log_chla_mean +
      depth"
    )))
  
  message(paste0(name, " model fitted and inspected"))
  
}

# Merge all validation pdfs
pdftools::pdf_combine(
  input = c(paste0("models/biomass/pielou_", new_res,
                   "_biomass_model_selection.pdf")),
  output = "models/biomass/pielou_partial_biomass_models_selection.pdf"
)

# Cleaning
file.remove(
  paste0("models/biomass/pielou_", new_res, 
         "_biomass_model_selection.pdf")
)

rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Other functional diversity model ----

# This is the model for biomass using convex hull diversity metrics

# Prepare data of the whole study area (all) for exploration
other_fd_data <- final_data %>%
  select(all_of(recurring_vars), starts_with("ch"), -starts_with("kde")) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.)))) # mean = 0; sd = 1 

# Variables distribution and outliers
other_vars <- setdiff(colnames(other_fd_data),
                      c("year", "latitude", "longitude", "haul_id")) 

par(mfrow = c(4, 4), mar = c(2, 4, 2, 1))

for(i in seq_along(other_vars)) {
  boxplot(other_fd_data[[other_vars[i]]], # biomass, chla and fishing are skewed
          horizontal = TRUE,
          main = other_vars[i])  
}

par(mfrow = c(1, 1))  

# Collinearity X
other_fd_data %>% # Pearson correlation
  select(where(is.numeric)) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq), abs(Freq) > 0.65)

# Correlated variables:
#Var1     Var2       Freq             Choice
#1      sst_mean  sbt_mean  0.8129230 sbt_mean
#2      sst_mean  sic_mean -0.8254482 
#3      sst_mean    sic_sd -0.8258865
#4      sic_mean    sic_sd  0.8743660 sic_mean
#8          tric   ch_fric  0.6975578 ch_fric

other_fd_data %>% select(
  ch_fric,
  ch_feve,
  tric,    
  teve,         
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  depth
) %>% vif() # tric shares a small collinearity with ch_fric

# Selection
other_fd_model <- model_selection(
  directory = "models/biomass",
  name = "other_fd_biomass",
  data = other_fd_data,  
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = total_biomass ~  
    ch_fric +
    ch_feve +      
    tric +
    teve +     
    #sst_mean +      
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Merge all coefficient files in one csv ----

csv_files <- list.files(
  path = "models/biomass",
  pattern = ".csv$",
  full.names = TRUE
)

combined_csv <- bind_rows(lapply(csv_files, read.csv)) |>
  filter(term !="(Intercept)")

write.csv(combined_csv,
          file = "models/biomass/biomass_models_coefficients.csv")

# Cleaning
file.remove(csv_files)

# End

## ------------------------------------------------------------------------ ----