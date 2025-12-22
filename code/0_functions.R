## File: 0_functions.R
## Purpose: create useful functions to streamline analysis and code
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

#

## Install and/or load packages function ----
install_load_function <- function(pkg){
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg))
  install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}

save(install_load_function, file = "tools/install_load_function.RData")

## Functional diversity function ----

# Description

# Required packages 
library(BAT)

# Function 
Fdiversity <- function(species_matrix, trait_matrix) {
  hv <<- BAT::kernel.build(comm = species_matrix,   
                           trait = trait_matrix,
                           method.hv = "gaussian", 
                           abund = TRUE,    
                           cores = 1) 
  
  Fric <- kernel.alpha(hv)
  Fdis <- kernel.dispersion(hv, func = "divergence", frac = 0.1)
  Feve <- kernel.evenness(hv)                     
  
  functional_diversity <- data.frame(haul_id = names(Fric), Fric, Feve, Fdis)
  save(functional_diversity, file = "data/intermediate/functional_diversity.RData")
}

save(Fdiversity, file = "tools/Fdiversity_function.RData")

## Model selection functions ----

# Description

# Required packages 
library(sdmTMB)
library(dplyr)
  
# Function for complete model selection (both fixed and random term)
model_selection <- function(directory,
                            name,
                            data,
                            mesh_cutoff,
                            family,
                            fixed_formula) {

vars <- all.vars(fixed_formula)

mesh <- make_mesh(data,
                  xy_cols = c("longitude", "latitude"),
                  cutoff = mesh_cutoff)

# Random-term selection ----

# 1) random error only
error_only <- sdmTMB(
  formula = fixed_formula, 
  data = data, 
  mesh = mesh,
  spatial = "off",
  family = family,
  reml = TRUE) 

# 2) 1 + spatial random fields 
spatial_random_fields <- update(error_only , spatial = "on")

# 3) 2 + anisotropy
spatial_random_fields_with_anisotropy <- update(spatial_random_fields,
                                                anisotropy = TRUE) 

# 4) 2 + spatiotemporal random fields ("iid" = temporally independent)
spatiotemporal_random_fields_iid <- update(spatial_random_fields,
                                           time = "year", spatiotemporal = "iid") 

# 5) 4 + anisotropy
spatiotemporal_random_fields_iid_with_anisotropy <- update(spatiotemporal_random_fields_iid,
                                                           anisotropy = TRUE) 

# 6) 2 + spatiotemporal random fields ("ar1" = temporally dependent) # !!! check
spatiotemporal_random_fields_ar1 <- update(spatial_random_fields, 
                                           time = "year", spatiotemporal = "ar1")

# 7) 6 + anisotropy
spatiotemporal_random_fields_ar1_with_anisotropy <- update(spatiotemporal_random_fields_ar1,
                                                           anisotropy = TRUE)

# List models
random_list <- list(
  error_only = error_only,
  spatial_random_fields = spatial_random_fields,
  spatial_random_fields_with_anisotropy = spatial_random_fields_with_anisotropy,
  spatiotemporal_random_fields_iid = spatiotemporal_random_fields_iid,
  spatiotemporal_random_fields_iid_with_anisotropy = spatiotemporal_random_fields_iid_with_anisotropy,
  spatiotemporal_random_fields_ar1 = spatiotemporal_random_fields_ar1,
  spatiotemporal_random_fields_ar1_with_anisotropy = spatiotemporal_random_fields_ar1_with_anisotropy
)

# Check sanity
random_sanity_check <- sapply(random_list, sanity)
bad_random <- names(random_list)[!as.logical(random_sanity_check["all_ok", ])]


# AIC
random_aic <- AIC(error_only,
                  spatial_random_fields ,
                  spatial_random_fields_with_anisotropy,
                  spatiotemporal_random_fields_iid,
                  spatiotemporal_random_fields_iid_with_anisotropy,
                  spatiotemporal_random_fields_ar1,
                  spatiotemporal_random_fields_ar1_with_anisotropy)

# Order
random_aic <- random_aic[order(random_aic$AIC), ] 

# Differences between each model and the next better one
random_aic$delta_AIC <- c(diff(random_aic$AIC), NA)

# Print for possible manual selection
print(random_aic)

# Select the best model which is at least > 2 units better than the following    
if (random_aic$delta_AIC[1] > 2) {  
  selected <- rownames(random_aic)[1]
} else if (random_aic$delta_AIC[2] > 2) {
  lowest_df <- min(random_aic[1:2, ]$df) # lowest degrees of freedom between the two
  selected <- rownames(random_aic)[1:2][random_aic$df[1:2] == lowest_df][1] # choose the lowest (or the first if equal)
} else {
  stop("Manual random term selectio required; use 'random_aic' for comparisons") 
} 

# Fixed-term selection ----

# Refit model with maximum likelihood (Zuur et al, 2009)
fit_full <- update(get(selected), reml = FALSE) 

# Remove each covariate to assess model changes
covariates <- vars[-1]

fixed_list <- list()
for (var in covariates) {
  to_remove <- as.formula(paste(". ~ . -", var))
  formula_new <- update(formula(fit_full), to_remove)
  fixed_list[[paste0("model_without_", var)]] <- sdmTMB( # update() doesn't work for this
    formula = formula_new,
    data = data,
    mesh = mesh,
    family = family,
    spatial = fit_full$spatial,
    spatiotemporal = fit_full$spatiotemporal,
    anisotropy = if (is.null(fit_full$anisotropy)) FALSE else fit_full$anisotropy,
    time = fit_full$time,
    reml = FALSE
  )
}

# Check sanity
fixed_sanity_check <- sapply(fixed_list, sanity)
bad_fixed <- names(fixed_list)[!as.logical(fixed_sanity_check["all_ok", ])]

# Compute AIC of all models
fixed_aic <- data.frame(
  model = names(fixed_list),
  AIC = sapply(fixed_list, AIC)
) |>
  mutate(
    delta_AIC = AIC - AIC(fit_full),
    var_removed = gsub("model_without_", "", model)
  ) |>
  arrange(desc(AIC))

# Keep variables that reduce AIC by > 2
vars_to_keep <- covariates[!covariates %in% fixed_aic$var_removed[fixed_aic$delta_AIC <= 2]]
response <- vars[1]

if (length(vars_to_keep) > 0) {
  final_fixed_formula <- formula(paste(response, " ~ ", paste(vars_to_keep, collapse = " + ")))
} else {
  final_fixed_formula <- formula(paste0(response, " ~ 1"))
}

# Final model
model <- sdmTMB( 
  formula = final_fixed_formula,
  data = data,
  mesh = mesh,
  family = family,
  spatial = fit_full$spatial,
  spatiotemporal = fit_full$spatiotemporal,
  anisotropy = if (is.null(fit_full$anisotropy)) FALSE else fit_full$anisotropy,
  time = fit_full$time,
  reml = TRUE
)

# Objects to return or save ----

# Start logging to file
log_file <- file.path(directory, paste0(name, "_model_selection_checks.txt"))
sink(log_file)

cat(paste0("\nMODEL SELECTION CHECKS (", name, " model)\n"))
cat("\n---------------------------------------------------------------------\n")
cat("\nRANDOM TERM\n")

# Report random-term sanity
if (length(bad_random) > 0) {
  cat("\nSanity check failed for the following random-term models:\n")
  cat(paste(bad_random, collapse = ", "), "\n")
} else {
  cat("\nAll models in the random-term selection passed sanity checks.\n")
}

# Selected random structure
cat("\nSelected random structure:", selected, "\n")
cat("\n---------------------------------------------------------------------\n")
cat("\nFIXED TERM\n")

# Report fixed-term sanity
if (length(bad_fixed) > 0) {
  cat("\nSanity check failed for the following fixed-term models:\n")
  cat(paste(bad_fixed, collapse = ", "), "\n")
} else {
  cat("\nAll models in the fixed-term selection passed sanity checks.\n")
}

# Selected fixed structure
cat("\nSelected fixed structure:\n")
cat(paste(deparse(final_fixed_formula), collapse = ""), "\n")
cat("\n---------------------------------------------------------------------\n")
cat("\nCOLLINEARITY AMONG COVARIATES\n")

# Print AIC tables
cat("\nRandom-term selection AIC table:\n")
print(random_aic)

cat("\nFixed-term selection AIC table:\n")
print(fixed_aic)

cat("\n---------------------------------------------------------------------\n")
cat("\nMODEL\n")

# Print model
cat(paste0("\n", name,  " model\n"))
print(model)
cat("\n---------------------------------------------------------------------\n")

# Stop logging
sink()

# Save its coefficients as excell
write.table(data.frame(tidy(model),
                      model = name,
                      response = response,
           paste0(directory, "/", name, "_model_coefficients.xlsx")))

# Return final model
return(model)

}

# Function for partial model selection (only random term)
model_selection_only_random <- function(directory, 
                                        name,
                                        data,
                                        mesh_cutoff,
                                        family,
                                        fixed_formula) {
  
  vars <- all.vars(fixed_formula)
  
  mesh <- make_mesh(data,
                    xy_cols = c("longitude", "latitude"),
                    cutoff = mesh_cutoff)
  
  # Random-term selection ----
  
  # 1) random error only
  error_only <- sdmTMB(
    formula = fixed_formula, 
    data = data, 
    mesh = mesh,
    spatial = "off",
    family = family,
    reml = TRUE) 
  
  # 2) 1 + spatial random fields 
  spatial_random_fields <- update(error_only , spatial = "on")
  
  # 3) 2 + anisotropy
  spatial_random_fields_with_anisotropy <- update(spatial_random_fields,
                                                  anisotropy = TRUE) 
  
  # 4) 2 + spatiotemporal random fields ("iid" = temporally independent)
  spatiotemporal_random_fields_iid <- update(spatial_random_fields,
                                             time = "year", spatiotemporal = "iid") 
  
  # 5) 4 + anisotropy
  spatiotemporal_random_fields_iid_with_anisotropy <- update(spatiotemporal_random_fields_iid,
                                                             anisotropy = TRUE) 
  
  # 6) 2 + spatiotemporal random fields ("ar1" = temporally dependent) 
  spatiotemporal_random_fields_ar1 <- update(spatial_random_fields, 
                                             time = "year", spatiotemporal = "ar1")
  
  # 7) 6 + anisotropy
  spatiotemporal_random_fields_ar1_with_anisotropy <- update(spatiotemporal_random_fields_ar1,
                                                             anisotropy = TRUE)
  
  # List models
  random_list <- list(
    error_only = error_only,
    spatial_random_fields = spatial_random_fields,
    spatial_random_fields_with_anisotropy = spatial_random_fields_with_anisotropy,
    spatiotemporal_random_fields_iid = spatiotemporal_random_fields_iid,
    spatiotemporal_random_fields_iid_with_anisotropy = spatiotemporal_random_fields_iid_with_anisotropy,
    spatiotemporal_random_fields_ar1 = spatiotemporal_random_fields_ar1,
    spatiotemporal_random_fields_ar1_with_anisotropy = spatiotemporal_random_fields_ar1_with_anisotropy
  )
  
  # Check sanity
  random_sanity_check <- sapply(random_list, sanity)
  bad_random <- names(random_list)[!as.logical(random_sanity_check["all_ok", ])]
  
  
  # AIC
  random_aic <- AIC(error_only,
                    spatial_random_fields,
                    spatial_random_fields_with_anisotropy,
                    spatiotemporal_random_fields_iid,
                    spatiotemporal_random_fields_iid_with_anisotropy,
                    spatiotemporal_random_fields_ar1,
                    spatiotemporal_random_fields_ar1_with_anisotropy)
  
  # Order
  random_aic <- random_aic[order(random_aic$AIC), ] 
  
  # Differences between each model and the next better one
  random_aic$delta_AIC <- c(diff(random_aic$AIC), NA)
  
  # Print for possible manual selection
  print(random_aic)
  
  # Select the best model which is at least > 2 units better than the following    
  if (random_aic$delta_AIC[1] > 2) {  
    selected <- rownames(random_aic)[1]
  } else if (random_aic$delta_AIC[2] > 2) {
    lowest_df <- min(random_aic[1:2, ]$df) # lowest degrees of freedom between the two
    selected <- rownames(random_aic)[1:2][random_aic$df[1:2] == lowest_df][1] # choose the lowest (or the first if equal)
  } else {
    stop("Manual random term selectio required; use 'random_aic' for comparisons") 
  } 
  
  # Final_model
  model <- get(selected)
  
  # Objects to return or save ----
  
  # Start logging to file
  log_file <- file.path(directory, paste0(name, "_model_selection_checks.txt"))
  sink(log_file)
  
  cat(paste0("\nMODEL SELECTION CHECKS (", name, " model)\n"))
  cat("\n---------------------------------------------------------------------\n")
  cat("\nRANDOM TERM\n")
  
  # Report random-term sanity
  if (length(bad_random) > 0) {
    cat("\nWARNING: model selection should not be trusted\n")
    
    cat("\nSanity check failed for the following random structures:\n")
    cat(paste(bad_random, collapse = ", "), "\n")
  } else {
    cat("\nAll models in the random-term selection passed sanity checks.\n")
  }
  
  # Selected random structure
  cat("\nSelected random structure:", selected, "\n")
  cat("\n---------------------------------------------------------------------\n")
  cat("\nFIXED TERM\n")
  
  # Print AIC tables
  cat("\nRandom-term selection AIC table:\n")
  print(random_aic)
  
  cat("\n---------------------------------------------------------------------\n")
  cat("\nMODEL\n")
  
  # Print model
  cat(paste0("\n", name,  " model\n"))
  print(model)
  cat("\n---------------------------------------------------------------------\n")
  
  # Stop logging
  sink()
  
  # Save its coefficients as excell
  write.table(data.frame(tidy(model),
                         model = name),
                        file = paste0(directory, "/", name, "_model_coefficients.csv"))
 
  # Return final model
  return(model)
}

# Function for manual model selection
manual_model_selection <- function(directory, 
                                   name,
                                   data,
                                   mesh_cutoff,
                                   family,
                                   random_selection,
                                   fixed_formula) {
  
  vars <- all.vars(fixed_formula)
  
  mesh <- make_mesh(data,
                    xy_cols = c("longitude", "latitude"),
                    cutoff = mesh_cutoff)
  
  # Random-term selection ----

  if (random_selection == "error_only") {
    
    # 1) random error only
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "off",
      family = family,
      reml = TRUE)
  }
   
  if (random_selection == "spatial_random_fields") {
    
    # 2) 1 + spatial random fields 
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "on",
      family = family,
      reml = TRUE)
  }
  
  if (random_selection == "spatial_random_fields_with_anisotropy") {
    
    # 3) 2 + anisotropy
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "on",
      family = family,
      anisotropy = TRUE,
      reml = TRUE)
  } 
  
  if (random_selection == "spatiotemporal_random_fields_iid") {
    
    # 4) 2 + spatiotemporal random fields ("iid" = temporally independent)
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "on",
      family = family,
      anisotropy = TRUE,
      time = "year",
      spatiotemporal = "iid",
      reml = TRUE)
  } 
  
  if (random_selection == "spatiotemporal_random_fields_iid_with_anisotropy") {
    
    # 5) 4 + anisotropy
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "on",
      family = family,
      time = "year",
      spatiotemporal = "iid",
      anisotropy = TRUE,
      reml = TRUE)
  } 
  
  if (random_selection == "spatiotemporal_random_fields_ar1") {
    
    # 6) 2 + spatiotemporal random fields ("ar1" = temporally dependent) 
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "on",
      family = family,
      anisotropy = TRUE,
      time = "year",
      spatiotemporal = "ar1",
      reml = TRUE)
  } 
  
  if (random_selection == "spatiotemporal_random_fields_ar1_with_anisotropy") {
    
    # 7) 6 + anisotropy 
    model <- sdmTMB(
      formula = fixed_formula, 
      data = data, 
      mesh = mesh,
      spatial = "on",
      family = family,
      anisotropy = TRUE,
      time = "year",
      spatiotemporal = "ar1",
      anisotropy = TRUE,
      reml = TRUE)
  } 
  
  
  # Check sanity
  random_sanity <- sanity(model)
  
  # Objects to return or save ----
  
  # Start logging to file
  log_file <- file.path(directory, paste0(name, "_model_selection_checks.txt"))
  sink(log_file)
  
  cat(paste0("\nMODEL SELECTION CHECKS (", name, " model)\n"))
  cat("\n---------------------------------------------------------------------\n")
  
  # Selected random structure
  cat("\nManually selected random structure:", random_selection, "\n")
  cat("\n---------------------------------------------------------------------\n")
  
  # Model sanity
  cat("\nModel passed sanity check:", random_sanity$all_ok, "\n")
  
  cat("\n---------------------------------------------------------------------\n")
  cat("\nMODEL\n")
  
  # Print model
  cat(paste0("\n", name,  " model\n"))
  print(model)
  cat("\n---------------------------------------------------------------------\n")
  
  # Stop logging
  sink()
  
  # Save its coefficients as excell
  write.table(data.frame(tidy(model),
                         model = name),
              file = paste0(directory, "/", name, "_model_coefficients.csv"))
  
  # Return final model
  return(model)
}

# Save
save(model_selection,
     model_selection_only_random,
     manual_model_selection,
     file = "tools/model_selection_function.RData")

## Model validation function ----

# Description

# Required packages
library(DHARMa)
library(dplyr)

# Function
model_validation <- function(model, data, name, directory) {
  set.seed(8)
  
  ## Residual simulation
  simulation <- simulate(model, nsim = 500, type = "mle-mvn") 
  residuals <- dharma_residuals(simulation, model, return_DHARMa = TRUE)
  
  ## Open pdf
  pdf(file = paste0(directory, "/", name, "_model_validation.pdf"),
      width = 8.5, height = 11)  # standard page size
  
  ## Basic information ----
  plot.new()
  
  text(
    x = 0, y = 1,
    labels = paste( 
      "Model validation report",
      "",
      paste("Model:", name),
      paste("Observations:", nrow(data)),
      paste("Date:", Sys.Date()),
      sep = "\n"
    ),
    adj = c(0, 1),
    cex = 1
  )
  
  
  ## Main diagnostics ----
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  DHARMa::testDispersion(residuals)
  DHARMa::testOutliers(residuals)
  DHARMa::plotQQunif(residuals)
  DHARMa::plotResiduals(residuals)
  
  
  ## Spatial autocorrelation ----
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  
  spatial_test <- DHARMa::testSpatialAutocorrelation(
    residuals,    
    x = data$longitude,
    y = data$latitude,
    plot = TRUE
  )
  
  pval <- round(spatial_test$p.value, 3)
  obs  <- round(spatial_test$statistic[1], 3)
  
  mtext(
    paste0(
      "Moran's I = ", obs,
      " | H1: Distance-based autocorrelation",
      " | p-value = ", pval
    ),
    side = 3, line = 0, cex = 0.8
  )
  
  ## Residuals vs predictors ----
  predictors <- setdiff(
    colnames(data),
    c(
      all.vars(formula(model))[1],
      "year", "depth", "longitude", "latitude", "haul_id"
    )
  )
  
  plots_per_page <- 16
  n_pages <- ceiling(length(predictors) / plots_per_page)
  
  for (p in seq_len(n_pages)) {
    idx <- ((p - 1) * plots_per_page + 1):
      min(p * plots_per_page, length(predictors))
    
    par(mfrow = c(4, 4), mar = c(4, 4, 4, 1))
    
    for (i in idx) {
      predictor <- data[[predictors[i]]]
      DHARMa::plotResiduals(residuals, form = predictor)
      mtext(paste0("Predictor: ", predictors[i]),
            side = 3, line = 0.5, cex = 0.7)
    }
  }
  
  par(mfrow = c(1, 1))
  
  ## Close PDF
  dev.off()
}

# Save
save(model_validation, file = "tools/model_validation_function.RData")


## Model prediction function ----

# Description

# Required packages
library(sdmTMB)
library(dplyr)
library(ggplot2)

# Function
model_prediction <- function(model, name, directory) {
  
  # Model variables
  fixed_vars <- all.vars(formula(model))[-1]
  response <- all.vars(formula(model))[1]
  
  # Model data
  data <- model.frame(model)
  
  # Log transform if the response is biomass
  if (family(model)$link == "log") {
    data[[response]] <- log(data[[response]])
  }
  
  # Loop through each variable in the fixed formula
  for (var in fixed_vars) {
    
    # Create new data: vary current var, others fixed at their mean
    nd <- data.frame(matrix(ncol = length(fixed_vars), nrow = 100))
    colnames(nd) <- fixed_vars
    
    for (v in fixed_vars) {
      if (v == var) {
        nd[[v]] <- seq(min(data[[v]], na.rm = TRUE),
                       max(data[[v]], na.rm = TRUE), length.out = 100)
      } else {
        nd[[v]] <- mean(data[[v]], na.rm = TRUE)
      }
    }
    
    # Add year column (required)
    nd$year <- 2014
    
    # Predict
    pred <- predict(model, newdata = nd, se_fit = TRUE, re_form = NA)
    
    # Create ggplot
    p <- ggplot() +
      geom_point(data = data, aes_string(x = var, y =  response), color = "gray", size = 0.5) +
      geom_line(data = pred, aes_string(x = var, y = "est"), color = "darkred", linewidth = 0.5) +
      geom_ribbon(data = pred,
                  aes_string(x = var,
                             ymin = "est - 1.96 * est_se",
                             ymax = "est + 1.96 * est_se"),
                  alpha = 0.4, fill = "darkred") +
      labs(title = paste(var, "effect on",  response),
           x = var, y = response) +
      theme_bw() +
      theme(axis.title = element_text(size = 14))
    
    # Save plot
    ggsave(paste0(directory, "/", name, "_prediction_", var, ".png"),
           plot = p, width = 8, height = 8, units = "cm")
    
    # Model coefficients 
    coefs <- tidy(model, conf.int = TRUE)[-1, ]
    
    # Add color column
    coefs$color <- ifelse(coefs$estimate >= 0, "Positive", "Negative")
    
    # Order terms for clean plotting
    coefs$term <- factor(coefs$term, levels = coefs$term[order(coefs$estimate)])
    
    # Plot coefficients
    p <- ggplot(coefs, aes(x = estimate, y = term, fill = color)) +
      geom_col(width = 0.6) +  # bar from 0 to estimate
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
      scale_fill_manual(values = c("Positive" = "forestgreen", "Negative" = "firebrick")) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      labs(x = "Estimate (± 95% CI)", y = "Covariate", title = "Model Coefficients") +
      theme_bw() +
      theme(legend.position = "none")
    
    # Save plot
    ggsave(paste0(directory, "/", name, "_coefficients_plot.png"),
           plot = p, width = 8, height = 8, units = "cm")
  }
}

# Save
save(model_prediction, file = "tools/model_prediction_function.RData")

## Grid plotting function ----

# Description

# Required packages 
c("dplyr",
   "ggplot2",
   "tidyr", # pivot to wide format
   "terra",
   "sf",
   "rnaturalearth")

# Required objects
land <- st_crop(ne_countries(scale = "medium", returnclass = "sf"),
                xmin = -10, xmax = 65, 
                ymin = 55, ymax = 90) |>
  st_transform(
    crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") 

theme_custom <- function() {
  theme(
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray", color = "gray45", linewidth = 1),
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 8)
  )
}

save(land, theme_custom, file = "tools/mapping_objects.RData")

# Function
grid_plot <- function(data, variable, unit) {
  
  # plotting devices
  load("tools/mapping_objects.RData")
  
  # remove extreme values (top and bottom 0.5%) from data
  lower <- quantile(data[[variable]], 0.005, na.rm = TRUE)  # 0.5th percentile
  upper <- quantile(data[[variable]], 0.995, na.rm = TRUE)  # 99.5th percentile
  data <- data %>% 
    filter(
      !is.na(.data[[variable]]),
      .data[[variable]] >= lower,
      .data[[variable]] <= upper
    ) 
  
  # create grid
  g <- data %>%
    st_make_grid(cellsize = 64.82,  # 35 nautical miles in km
                 what = "polygons",
                 square = TRUE, 
                 offset = c(-483, -573)) %>% # center the cells
    st_sf()
  
  # Spatial join
  grid_r <- st_join(g, data, left = FALSE)
  
  # Plot limits
  xlim <- st_bbox(data)[c(1, 3)] + c(-65, 65)
  ylim <- st_bbox(data)[c(2, 4)] + c(-65, 65)
  
  # Mean per cell
  grid_mean <- grid_r %>%
    group_by(geometry) %>%
    summarize(x = mean(.data[[variable]], na.rm = TRUE)) %>%
    st_as_sf()
  
  # Range
  max_val <- ceiling(max(grid_mean$x, na.rm = TRUE) * 100)/100
  min_val <- floor(min(grid_mean$x, na.rm = TRUE) * 100)/100
  
  # Plot
  ggplot() +
    geom_sf(data = grid_mean, aes(fill = x), color = NA) +
    geom_sf(data = land, fill = "gray45", color = "gray45") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                         values = scales::rescale(c(min_val, max_val)),
                         limits = c(min_val, max_val),
                         name = unit,
                         breaks = seq(min_val, max_val, length.out = 3)) +
    labs(x = "\nLongitude", y = "Latitude\n", title = variable) +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 10, ticks.linewidth = .1)) +
    theme_custom()
}

save(grid_plot, file = "tools/grid_plot_function.RData")



## VIF function ----

## Description
# The function computes the variance inflation factor of all the variables 
# in a data frame, returning those above 2. This is a strict threshold 
# that avoids the hiding of weak ecological effects. 
# See Zuur et al., 2010
# https://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2009.00001.x

## Function
vif <- function(df) {
  vif_values <- numeric(ncol(df))
  names(vif_values) <- colnames(df)
  
  for (i in seq_along(df)) {
    y <- df[[i]]
    x <- df[, -i, drop = FALSE]
    r2 <- summary(lm(y ~ ., data = x))$r.squared
    vif_values[i] <- 1 / (1 - r2)
  }
  vif_vals <- vif_values[vif_values >= 2]
  if (length(vif_vals > 0)) {
    message("VIF values >= 2 detected:")
    return(vif_vals)
    } else {message("All VIF values < 2")}

}

save(vif, file = "tools/vif_function.RData")
