## File: 0_functions.R
## Purpose: create useful functions to streamline analysis and code
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

#
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

## Model selection function ----

# Description

# Required packages 
library(sdmTMB)
library(performance)
library(writexl)
library(dplyr)
  
# Function
model_selection <- function(directory, name, data, mesh_cutoff, family, fixed_formula) {

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
  formula_new <- as.formula(paste("~ . -", var))
  fixed_list[[paste0("model_without_", var)]] <- update(fit_full, formula = update(formula(fit_full), formula_new))
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

# Check collinearity of the final variables in the model
collinearity_model <- update(fit_full, formula = final_fixed_formula)

# Check multicollinearity using performance (works only on ML models)
vif <- performance::check_collinearity(collinearity_model) 

problematic <- vif[vif$VIF > 2, ] # threshold from Zuur et al., 2010

# Final_model
model <- update(collinearity_model, reml = TRUE)

# Objects to return or save ----

# Save its coefficients as excel file 
write_xlsx(
  data.frame(tidy(model),
             area = name,
             response = response),
  path = paste0(directory, "/", name, "_model_coefficients.xlsx")
)

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

# Collinearity
if (length(rownames(problematic)) > 0) {
  cat("\nPotential collinearity issues (VIF >= 2) with:\n")
  cat(paste(problematic$Term, collapse = ", "), "\n")
} else {
  cat("\nNo collinearity detected among covariates (VIF < 2).\n")
}

cat("\n---------------------------------------------------------------------\n")
cat("\nAIC SCORES\n")

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

# Save fnal model
save(model,
     file = paste0(directory, "/", name, "_model.RData"))

# Return final model
return(model)
}

# Save
save(model_selection, file = "tools/model_selection_function.RData")

## Model validation function ----

# Description

# Required packages
library(DHARMa)
library(dplyr)

# Function
model_validation <- function(model, data, directory) {
  set.seed(8)
  
  # Residual
  simulation <- simulate(model, nsim = 1000, type = "mle-mvn") 
  residuals <- dharma_residuals(simulation, model, return_DHARMa = TRUE)
  
  # Spatial autocorrelation detection with Moran's I test
  png(paste0(directory, "/residuals_vs_space.png"),
      width = 12, height = 16, unit = "cm", res = 400)
  
  spatial_test <- DHARMa::testSpatialAutocorrelation(residuals,    
                                                     x = data$longitude,
                                                     y = data$latitude,
                                                     plot = TRUE)
  
  pval <- round(spatial_test$p.value * 1000)/1000
  mtext(paste0("p-value = ", pval, "    H1: Distance-based autocorrelation"),
        side = 3, line = .7, cex = .7)
  
  dev.off()
  
  # Dispersion
  png(paste0(directory, "/dispersion.png"),
      width = 15, height = 10, unit = "cm", res = 400)
  DHARMa::testDispersion(residuals)
  dev.off()
  
  # Outliers
  png(paste0(directory, "/outliers.png"),
      width = 15, height = 10, unit = "cm", res = 400)
  DHARMa::testOutliers(residuals) # here outliers are all values outside the simulation envelope (depends on nsim)
  dev.off()
  
  # Uniformity
  png(paste0(directory, "/uniformity.png"),
      width = 10, height = 12, unit = "cm", res = 400)
  DHARMa::plotQQunif(residuals) 
  dev.off()
  
  # Residuals vs fitted 
  png(paste0(directory, "/residuals_vs_fitted.png"),
      width = 10, height = 12, unit = "cm", res = 400)
  DHARMa::plotResiduals(residuals) 
  dev.off()
  
  # Residuals vs predictors
  predictors <- setdiff(colnames(data),
                        c(all.vars(formula(model))[1],
                          "year", "depth", "longitude", "latitude", "haul_id"))
  
  index <- paste0(1:length(predictors), ")")
  
  png(paste0(directory, "/residuals_vs_predictors.png"),
      width = 17, height = 35, unit = "cm", res = 400)
  par(mfrow = c(5, 3), mar = c(4, 4, 8, 4))
  for(i in 1:length(predictors)){
    predictor <- data[[predictors[i]]]
    DHARMa::plotResiduals(residuals, form = predictor) 
    mtext(paste0("Pred.: ", predictors[i]),
          side = 3, line = 1, cex = .7)
    mtext(index[i],
          side = 3, line = 1, cex = .7, adj = 0)
  }
  par(mfrow = c(1, 1))
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
model_prediction <- function(model, directory) {
  
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
    ggsave(paste0(directory, "/prediction_", var, ".png"),
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
    ggsave(paste0(directory, "/coefficients_plot.png"),
           plot = p, width = 8, height = 8, units = "cm")
  }
}

# Save
save(model_prediction, file = "tools/model_prediction_function.RData")
