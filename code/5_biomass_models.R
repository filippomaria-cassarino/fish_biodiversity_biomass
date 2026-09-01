## File: 5_biomass_models.R
## Purpose: model the response of biomass to diversity, warming, and fishing
## Author: Filippomaria Cassarino
## Date: 31 Aug 2026
## Notes ----

## Library ---- 

# Import data, objects and functions
final_data <- read.csv("data/final/final_data.csv")

# Load functions
load("tools/model_selection_functions.RData")
load("tools/model_prediction_function.RData")
load("tools/vif_function.RData")
load("tools/install_load_packages_function.RData")

# Load required packages
install_load_packages(c(
  "dplyr",
  "ggplot2",
  "tidyr",   
  "sdmTMB", 
  "DHARMa",
  "pdftools",
  "readr"
))

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
  "log_chla_mean",
  "sbt_sd",
  "sst_sd",
  "sic_sd", 
  "log_chla_sd",
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

# Correlated variables (arrow represents choice):

# sst_mean   ->   sbt_mean
# sic_mean   ->   sic_sd
# log_chla_mean -> log_chla_sd 

# Check diversity correlations
general_data %>% 
  select(teve, tric, kde_feve, kde_fric) %>%
  cor(use = "pairwise.complete.obs") %>%
  replace(lower.tri(., diag = TRUE), NA) %>%
  as.data.frame.table(responseName = "Freq") %>%
  filter(!is.na(Freq))

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
  #log_chla_sd,
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
    #log_chla_sd +
    depth
  )

# Prediction
model_prediction(model = general_model,
                 vars = c("teve", "tric", "sbt_mean"),
                 name = "general_biomass",
                 directory = "figures/model_figures")

# Refit with only one diversity metric at a time

vars <- c("kde_fric", "kde_feve", "teve", "tric")

for (var in vars) {
  
  model <- model_selection(
    directory = "models/biomass",
    name = paste0("only_", var, "_biomass"),
    data = general_data,  
    mesh_cutoff = 60, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      "total_biomass  ~ ", var, " + 
      sbt_mean +
      sic_mean +
      sbt_sd +
      sst_sd +        
      log_chla_mean +
      depth"
      )))
  
}

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

# sst_mean  ->  sbt_mean
# sic_mean  <-   sic_sd 
# log_chla_mean -> log_chla_sd

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
  #log_chla_sd,
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
    #log_chla_sd +
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
divide_arc <- quantile(final_data$sbt_mean, 0.3, na.rm = TRUE)

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

# sst_mean  ->  sic_mean
#  sic_mean   <-   sic_sd
# log_chla_mean <- log_chla_sd

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
  #log_chla_sd,
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
    #log_chla_sd +
    depth
)

# The residuals suggest potential non-linearities, but adding a smoothers to
# non-linear drivers does not solve the issues.

## ------------------------------------------------------------------------ ----
## Boreal model ----

# This model focuses on Boreal areas

# Identify Boreal data
divide_bor <- quantile(final_data$sbt_mean, 0.7, na.rm = TRUE)

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

# sst_mean  ->  sbt_mean
#  sic_mean   /   sic_sd - sic_mean
# log_chla_mean <- log_chla_sd 

boreal_data %>% select( # VIF
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
  #log_chla_sd,
  depth
) %>% vif() # all VIF < 2

# How many sic values?
table(boreal_data$sic_mean > 0)
prop.table(table(boreal_data$sic_mean > 0))

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
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    #log_chla_sd +
    depth
)

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

# sst_mean  ->  sbt_mean 
# sic_mean   <-   sic_sd
# log_chla_mean <- log_chla_sd

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
  #log_chla_sd,
  depth
) %>% vif() 

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
    #log_chla_sd +
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

# sst_mean  ->  sbt_mean
# sic_mean   <-   sic_sd
# log_chla_mean <- log_chla_sd

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
  #log_chla_sd,
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
    #log_chla_sd +
    depth
)

## ------------------------------------------------------------------------ ----
## 2017_2022 model ----

# This model focuses on the third (cooler) time period identifies by
# Eriksen et al., 2025. https://doi.org/10.1038/s41598-025-96964-x

data_2017_2022 <- filter(final_data, year > 2016) %>% # cooler period
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

# sst_mean  ->  sbt_mean
# sic_mean   <-   sic_sd 
# log_chla_mean <- log_chla_sd

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
  #log_chla_sd,
  depth
) %>% vif() # sic_mean is above two

# Select 2017_2022 model (manual because of an error)
model_2017_2022 <- manual_model_selection( 
  directory = "models/biomass",
  name = "2017_2022_biomass",
  data = data_2017_2022,
  mesh_cutoff = 55, # sanity check fails at 60, so we use 55
  family = lognormal(link = "log"),
  random_selection = "spatiotemporal_random_fields_ar1",
  fixed_formula = total_biomass ~  
    kde_fric +
    kde_feve +
    tric +
    teve + 
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    #sic_mean + 
    #sic_sd +       
    log_chla_mean +
    #log_chla_sd +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Partial biomass models ----

# These models investigate the effect of the progressive removal of 
# dominant species from the total_biomass. 
# Note that data exploration is the same as general_model

new_res <- paste0("remove_", 1:10)

# Prepare data
biomass_data <- final_data %>%
  select(         
    all_of(recurring_vars),
    starts_with("remove"),
    log_sum_fishing
  ) 

# Responses distribution
par(mfrow = c(3, 4), mar = c(2, 4, 2, 1)) 

for (i in 1:length(new_res)) { 
  
  plot(density(log(biomass_data[[paste0(new_res[i], "_biomass")]])),
       main = paste0(new_res[i], "_biomass"))
}

par(mfrow = c(1, 1)) 

# Collinearity (vif only)
for (i in 1:length(new_res)) { 
  
  message(new_res[i])
  
  x <- final_data %>%
    select(
      paste0(new_res[i], "_tric"),
      paste0(new_res[i], "_teve"),
      paste0(new_res[i], "_kde_fric"),
      paste0(new_res[i], "_kde_feve"),
      "sbt_mean",
      "sbt_sd",
      "sst_sd",
      #"sic_mean",
      "depth",
      "log_chla_mean",
      "log_sum_fishing"
    ) %>%
    drop_na() %>%
    vif()
  
  print(x)
  
}

# In some data sets sic_mean and/or sbt_mean surpass the VIF threshold 
# We remove sic_mean from all models to keep them comparable and remove all
# excessive correlations. sic_mean was not a driver in the general model and 
# does not become an important one in the partial biomass models if kept.

# Models
for (i in 1:length(new_res)) { 
  
  # Identify new response
  response <- paste0(new_res[i], "_biomass")
  
  data <- biomass_data %>%
    select(
      paste0(new_res[i], c(
        "_biomass",
        "_kde_fric",
        "_kde_feve",
        "_tric",
        "_teve")),
      sbt_mean,
      sbt_sd,
      sst_sd,
      #sic_mean,        
      log_chla_mean,
      depth,
      log_sum_fishing, 
      haul_id,
      year, 
      latitude,
      longitude) %>%
    drop_na() %>%
    filter(.data[[response]] > 0) %>%
    mutate(
      across(
        -c(
          ends_with("biomass"),
          haul_id,
          year, 
          latitude,
          longitude
        ), ~ as.numeric(scale(.)))) 
  
  # Fitting 
  model <- model_selection(
    directory = "models/biomass",
    name = response,
    data = data,
    mesh_cutoff = 60, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i], "_biomass  ~", 
      new_res[i], "_kde_fric +",
      new_res[i], "_kde_feve +",
      new_res[i], "_tric +",
      new_res[i], "_teve + 
      sbt_mean +
      sbt_sd +
      sst_sd +        
      log_chla_mean +
      depth +
      log_sum_fishing"
      )))
  
  # Return remove_10_model to visualize predictions
  if(formula(model)[[2]] == "remove_10_biomass") {
    
    remove_10_biomass_model <- model
    
    return(remove_10_biomass_model)
  }
  
  message(
    "\n======================================================================\n",
    paste0(new_res[i], " model fitted and inspected"),
    "\n======================================================================\n"
    )
  
}

# Check remove_10 model predictions for tric, teve, and sbt_mean
model_prediction(model = remove_10_biomass_model,
                 name = "remove_10_biomass_model",
                 vars = c("remove_10_teve", "remove_10_tric", "sbt_mean"),
                 directory = "figures/model_figures")

# Residuals indicate the presence of some non-linear relationships. This seems
# to be potentially problematic only for teve and tric, so we refit the models
# with smoothers and check 

remove_10_biomass_gamm_model <- model_selection(
  directory = "models/biomass",
  name = "remove_10_biomass_gamm",
  data = remove_10_biomass_model$data,
  mesh_cutoff = 60, 
  family = lognormal(link = "log"),
  fixed_formula = remove_10_biomass  ~ 
    remove_10_kde_fric +
    remove_10_kde_feve +
    s(remove_10_tric) +
    s(remove_10_teve) + 
      sbt_mean +
      sbt_sd +
      sst_sd +        
      log_chla_mean +
      depth +
      log_sum_fishing
  )

# Check remove_10 model (gamm) predictions for tric and teve
model_prediction(model = remove_10_biomass_gamm_model,
                 name = "remove_10_biomass_gamm_model",
                 vars = c("remove_10_teve", "remove_10_tric"),
                 directory = "figures/model_figures")

# Overall, the gamm has a better fit, but the difference in the predictions
# is minor. Given that the aim is to compare across models, we can accept a
# slightly worse fit (due to minor non-linearities) to maintain high 
# comparabiliy

# Merge all pdfs
pdf_combine(
  input = c(
    paste0("models/biomass/", new_res, "_biomass_model_selection.pdf"),
    "models/biomass/remove_10_biomass_gamm_model_selection.pdf"),
  output = "models/biomass/partial_biomass_models_selection.pdf"
)

# Cleaning
file.remove(
  c(
    paste0("models/biomass/", new_res, "_biomass_model_selection.pdf"),
    "models/biomass/remove_10_biomass_gamm_model_selection.pdf",
    "models/biomass/remove_10_biomass_gamm_model_coefficients.csv")
)

# Partial biomass models using the same diversity metrics  
partial_data <- final_data %>%
  select(all_of(recurring_vars),
         log_sum_fishing,
         ends_with("biomass")) %>%
  mutate(across(
    ends_with("biomass"),
    ~ replace(., . < 1e-5, NA)
  )) %>%
  drop_na() %>%
  mutate(
    across(
    -c(
      ends_with("biomass"),
      haul_id,
      total_biomass,
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

for (i in 1:length(new_res)) { 
  
  # Fitting 
  model <- model_selection(
    directory = "models/biomass",
    name = paste0(new_res[i], "_samediv_biomass"),
    data = partial_data,
    mesh_cutoff = 60, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i], "_biomass  ~
      kde_fric +
      kde_feve +
      tric +
      teve + 
      sbt_mean +
      sbt_sd +
      sst_sd +        
      log_chla_mean +
      depth +
      log_sum_fishing"
    )))
  
  message(
    "\n======================================================================\n",
    paste0(new_res[i], " (same diversity) model fitted and inspected"),
    "\n======================================================================\n"
  )
  
}

# Merge all pdfs
pdf_combine(
  input = paste0("models/biomass/", new_res, "_samediv_biomass_model_selection.pdf"),
  output = "models/biomass/partial_biomass_same_diversity_models_selection.pdf"
)

# Cleaning
file.remove(
  paste0("models/biomass/", new_res, "_samediv_biomass_model_selection.pdf"))


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

# sst_mean   -> sbt_mean
# sic_mean   <-   sic_sd
#  log_chla_mean <- log_chla_sd

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
  #log_chla_sd,
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
    #log_chla_sd +
    depth
)

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Other functional diversity model ----

# This is the model for biomass using convex hull diversity metrics

# Prepare data of the whole study area (all) for exploration
other_fd_data <- final_data %>%
  select(all_of(recurring_vars), starts_with("ch_"), -starts_with("kde")) %>%
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

# sst_mean   ->  sbt_mean
# sic_mean   <-   sic_sd
# log_chla_mean <- log_chla_sd
# tric  ->   ch_fric 

other_fd_data %>% select(
  ch_fric,
  ch_feve,
  #tric,    
  teve,         
  #sst_mean,       
  sbt_mean,       
  sic_mean, 
  sbt_sd,
  sst_sd,
  #sic_sd, 
  log_chla_mean,
  #log_chla_sd,
  depth
) %>% vif() # tric shares a small collinearity with ch_fric, but they are 
            # both maintained to favor comparisons

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
    #tric +
    teve +     
    #sst_mean +      
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    #sic_sd + 
    log_chla_mean +
    #log_chla_sd +
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

csv_files <- csv_files[
  csv_files != "models/biomass/biomass_models_coefficients.csv"
]

combined_csv <- bind_rows(lapply(csv_files, read.csv)) |>
  filter(term !="(Intercept)")

write.csv(combined_csv,
          file = "models/biomass/biomass_models_coefficients.csv")

# Cleaning
file.remove(csv_files)

rm(list = ls())

# End

## ------------------------------------------------------------------------ ----