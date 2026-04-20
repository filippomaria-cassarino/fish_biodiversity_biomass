## Name: 0_functions.R
## Purpose: create functions to streamline code
## Author: Filippomaria Cassarino
## Date: 16 Apr 2026

## Notes ----

# Grid plotting function is not finished yet

## Install and/or load packages function ----

# Description: installs missing packages and loads required ones

# Function
install_load_packages <- function(pkg){
  
  # identify missing packages
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  
  # install missing packages
  if (length(new.pkg) > 0) {
    
    install.packages(new.pkg, dependencies = TRUE)
  }
    
  # load required packages
  sapply(pkg, require, character.only = TRUE)
}

save(install_load_packages, file = "tools/install_load_packages_function.RData")

## Variance inflation factor (VIF) function ----

# Description: computes the VIF of all the variables in a data frame,
# returning those above 2. This is a strict threshold 
# that avoids the hiding of weak ecological effects. 
# See Zuur et al., 2010
# https://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2009.00001.x

# Function
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
    
  } else {
    
    message("All VIF values < 2")
    
    return(vif_values)}
}

save(vif, file = "tools/vif_function.RData")

## Model selection functions ----

# Description: random_model_selection builds a mesh to reduce 
# spatial complexity, fits glmms with different random structures,
# selects best random structure based of AIC 
# and opens a pdf detailing the selection process, model summary,
# and model validation.

# manual_model_selection does the same but allows to select random term

# Required packages 
c(
  "dplyr",
  "sdmTMB",
  "DHARMa"
)

# Selection function 
model_selection <- function(directory,     # directory to save results
                            name,          # model name
                            data,          # data
                            mesh_cutoff,   # minimum size for mesh in km
                            family,        # distribution family
                            fixed_formula) # fixed effects formula
{
  
  # Build mesh
  mesh <- sdmTMB::make_mesh(
    data,
    xy_cols = c("longitude", "latitude"),
    cutoff = mesh_cutoff
  )
  
  # Fit models with different random terms ----
  
  # 1) random error only
  error_only <- sdmTMB(
    formula = fixed_formula, 
    data = data, 
    mesh = mesh,
    spatial = "off",
    family = family,
    reml = TRUE
  ) 
  
  # 2) 1 + spatial random fields 
  spatial_random_fields <- update(
    error_only,
    spatial = "on"
  )
  
  # 3) 2 + anisotropy
  spatial_random_fields_with_anisotropy <- update(
    spatial_random_fields,
    anisotropy = TRUE
  ) 
  
  # 4) 2 + spatiotemporal random fields ("iid" = temporally independent)
  spatiotemporal_random_fields_iid <- update(
    spatial_random_fields,
    time = "year", 
    spatial = "off",
    spatiotemporal = "iid"
  ) 
  
  # 5) 4 + anisotropy
  spatiotemporal_random_fields_iid_with_anisotropy <- update(
    spatiotemporal_random_fields_iid,
    anisotropy = TRUE
  ) 
  
  # 6) 2 + spatiotemporal random fields ("ar1" = temporally dependent) 
  spatiotemporal_random_fields_ar1 <- update(
    spatial_random_fields, 
    spatial = "off",
    time = "year",
    spatiotemporal = "ar1"
  )
  
  # 7) 6 + anisotropy
  spatiotemporal_random_fields_ar1_with_anisotropy <- update(
    spatiotemporal_random_fields_ar1,
    anisotropy = TRUE
  )
  
  # list models
  random_list <- list(
    error_only = 
      error_only,
    spatial_random_fields = 
      spatial_random_fields,
    spatial_random_fields_with_anisotropy = 
      spatial_random_fields_with_anisotropy,
    spatiotemporal_random_fields_iid = 
      spatiotemporal_random_fields_iid,
    spatiotemporal_random_fields_iid_with_anisotropy = 
      spatiotemporal_random_fields_iid_with_anisotropy,
    spatiotemporal_random_fields_ar1 = 
      spatiotemporal_random_fields_ar1,
    spatiotemporal_random_fields_ar1_with_anisotropy = 
      spatiotemporal_random_fields_ar1_with_anisotropy
  )
  
  # check sanity
  random_sanity_check <- sapply(random_list, sanity)
  
  # identify potential bad models
  bad_random <- names(random_list)[!as.logical(random_sanity_check["all_ok", ])]
  
  # Calculate AIC and ΔAIC ----
  
  # calculate AIC
  random_aic <- AIC(error_only,
                    spatial_random_fields,
                    spatial_random_fields_with_anisotropy,
                    spatiotemporal_random_fields_iid,
                    spatiotemporal_random_fields_iid_with_anisotropy,
                    spatiotemporal_random_fields_ar1,
                    spatiotemporal_random_fields_ar1_with_anisotropy)
  
  # order 
  random_aic <- random_aic[order(random_aic$AIC), ] 
  
  # compute ΔAIC between each model and the next better one
  random_aic$delta_AIC <- c(diff(random_aic$AIC), NA)
  
  # print for potential manual selection
  print(random_aic)
  
  # Select ----
  
  #   if the 1st and 2nd models have ΔAIC ≥ 2, select 1st  
  if (random_aic$delta_AIC[1] >= 2) { 
    
    selected <- rownames(random_aic)[1]
    
    # if the 1st and 2nd have ΔAIC < 2 and 2nd and 3rd have ΔAIC ≥ 2,
    # select among the first 2 the one with the fewest degrees of freedom 
    # or 1st if they have the same degrees of freedom 
  } else if (random_aic$delta_AIC[2] > 2) {
    
    lowest_df <- min(random_aic[1:2, ]$df) 
    
    selected <- rownames(random_aic)[1:2][random_aic$df[1:2] == lowest_df][1] 
    
    # if the 2nd and 3rd have ΔAIC < 2, manual selection is required
  } else {
    
    stop(
      "Manual random term selectio required: use 'random_aic' for comparisons,
    then 'manual_random_model_selection' to fit the model, specifying the 
    structure through the argument 'random_selection' and the structure names
    from 'random_aic'."
    ) 
  } 
  
  # selected model
  model <- get(selected)
  
  # Report general information ----
  
  # model data
  data <- model.frame(model)
  
  # some calculations have a random component
  set.seed(8)
  
  # residual simulation for validation
  simulation <- simulate(
    model, 
    nsim = 500, 
    type = "mle-mvn"
  ) 
  
  residuals <- sdmTMB::dharma_residuals(
    simulation,
    model,
    return_DHARMa = TRUE
  )
  
  # open pdf
  pdf(file = paste0(
    directory,
    "/", 
    name, 
    "_model_selection.pdf"),
    width = 8.5, height = 11)  # standard page size
  
  # model information
  text_output <- capture.output({
    
    cat("-------------------------------------------------------------------\n")
    cat("MODEL SELECTION AND VALIDATION REPORT (", name, "model )\n")
    cat("-------------------------------------------------------------------\n\n")
    
    cat(paste0("## OBSERVATIONS:", nrow(data), "\n\n"))
    
    if (length(bad_random) > 0) {
      cat("## SANITY: warning - sanity check failed for the following models:\n")
      cat(paste(bad_random, collapse = ",\n "), "\n")
    } else {
      cat("## SANITY: OK - all models in the selection process passed sanity checks.\n")
    }
    
    cat("\n## SELECTED RANDOM TERM:", selected, "\n\n")
    
    cat("## MODEL SUMMARY (", name, "model )\n\n")
    print(model)
    
    cat("\n## SELECTION STATISTICS\n")
    print(random_aic)
    
  })
  
  par(mar = c(1, 1, 1, 1), oma = c(2, 2, 2, 2))
  
  plot.new()
  text(
    x = 0, y = 1,
    labels = paste(text_output, collapse = "\n"),
    adj = c(0, 1),
    family = "mono",
    cex = 0.8
  )
  
  # Validation ----
  # diagnostics plots
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 2, 2, 2))
  
  DHARMa::testDispersion(residuals)
  DHARMa::testOutliers(residuals)
  DHARMa::plotQQunif(residuals)
  DHARMa::plotResiduals(residuals)
  
  # spatial autocorrelation
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1), oma = c(2, 2, 2, 2))
  
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
  
  # residuals vs predictors
  predictors <- setdiff(
    colnames(data),
    c(
      all.vars(formula(model))[1],
      "year",
      "depth",
      "longitude",
      "latitude",
      "haul_id",
      "_sdmTMB_time"
    )
  )
  
  plots_per_page <- 16
  n_pages <- ceiling(length(predictors) / plots_per_page)
  
  for (p in seq_len(n_pages)) {
    idx <- ((p - 1) * plots_per_page + 1):
      min(p * plots_per_page, length(predictors))
    
    par(mfrow = c(4, 4), mar = c(4, 4, 4, 1), oma = c(2, 2, 2, 2))
    
    for (i in idx) {
      predictor <- data[[predictors[i]]]
      DHARMa::plotResiduals(residuals, form = predictor)
      mtext(paste0("Predictor: ", predictors[i]),
            side = 3, line = 0.5, cex = 0.7)
    }
  }
  
  par(mfrow = c(1, 1))
  
  # Objects to return ----
  
  # close PDF
  dev.off()
  
  # save coefficients
  write.table(
    data.frame(tidy(model), model = name),
    file = paste0(directory, "/", name, "_model_coefficients.csv"),
    sep = ","
  )
  
  # return selected model
  return(model)
  
}

# Manual selection function
manual_model_selection <- function(directory, 
                                   name,
                                   data,
                                   mesh_cutoff,
                                   family,
                                   random_selection, # select random structure
                                   fixed_formula) 
{
  
  mesh <- make_mesh(
    data,
    xy_cols = c("longitude", "latitude"),
    cutoff = mesh_cutoff
  )
  
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
      spatial = "off",
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
      spatial = "off",
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
      spatial = "off",
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
      spatial = "off",
      family = family,
      anisotropy = TRUE,
      time = "year",
      spatiotemporal = "ar1",
      anisotropy = TRUE,
      reml = TRUE)
  } 
  
  # check sanity
  random_sanity <- sanity(model)
  
  if (random_sanity["all_ok"] == FALSE) {
    stop("Sanity check failed for model.", call. = FALSE)
  }
  
  # Report general information ----
  
  # model data
  data <- model.frame(model)
  
  # some calculations have a random component
  set.seed(8)
  
  # residual simulation for validation
  simulation <- simulate(
    model, 
    nsim = 500, 
    type = "mle-mvn"
  ) 
  
  residuals <- sdmTMB::dharma_residuals(
    simulation,
    model,
    return_DHARMa = TRUE
  )
  
  # open pdf
  pdf(file = paste0(
    directory,
    "/", 
    name, 
    "_model_selection.pdf"),
    width = 8.5, height = 11)  # standard page size
  
  # model information
  text_output <- capture.output({
    
    cat("-------------------------------------------------------------------\n")
    cat("MODEL SELECTION AND VALIDATION REPORT (", name, "model )\n")
    cat("-------------------------------------------------------------------\n\n")
    
    cat(paste0("## OBSERVATIONS:", nrow(data), "\n\n"))
    
    cat("## SANITY: OK - the model passed sanity check.\n\n")
    
    cat("## MANUALLY SELECTED RANDOM TERM:", random_selection, "\n\n")
    
    cat("## MODEL SUMMARY (", name, "model )\n\n")
    print(model)
    
  })
  
  par(mar = c(1, 1, 1, 1), oma = c(2, 2, 2, 2))
  
  plot.new()
  text(
    x = 0, y = 1,
    labels = paste(text_output, collapse = "\n"),
    adj = c(0, 1),
    family = "mono",
    cex = 0.8
  )
  
  # Validation ----
  # diagnostics plots
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(2, 2, 2, 2))
  
  DHARMa::testDispersion(residuals)
  DHARMa::testOutliers(residuals)
  DHARMa::plotQQunif(residuals)
  DHARMa::plotResiduals(residuals)
  
  # spatial autocorrelation
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1), oma = c(2, 2, 2, 2))
  
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
  
  # residuals vs predictors
  predictors <- setdiff(
    colnames(data),
    c(
      all.vars(formula(model))[1],
      "year",
      "depth",
      "longitude",
      "latitude",
      "haul_id",
      "_sdmTMB_time"
    )
  )
  
  plots_per_page <- 16
  n_pages <- ceiling(length(predictors) / plots_per_page)
  
  for (p in seq_len(n_pages)) {
    idx <- ((p - 1) * plots_per_page + 1):
      min(p * plots_per_page, length(predictors))
    
    par(mfrow = c(4, 4), mar = c(4, 4, 4, 1), oma = c(2, 2, 2, 2))
    
    for (i in idx) {
      predictor <- data[[predictors[i]]]
      DHARMa::plotResiduals(residuals, form = predictor)
      mtext(paste0("Predictor: ", predictors[i]),
            side = 3, line = 0.5, cex = 0.7)
    }
  }
  
  par(mfrow = c(1, 1))
  
  # Objects to return ----
  
  # close PDF
  dev.off()
  
  # save coefficients
  write.table(
    data.frame(tidy(model), model = name),
    file = paste0(directory, "/", name, "_model_coefficients.csv"),
    sep = ","
  )
  
  # return selected model
  return(model)
}


# Save
save(model_selection,
     manual_model_selection,
     file = "tools/model_selection_functions.RData")

## Model prediction function ----

# Description: plots model predictions for each variable in the model

# Required packages
c(
  "dplyr",
  "ggplot2",
  "sdmTMB"
  )

# Function
model_prediction <- function(model, name, directory, vars) {
  
  # model variables
  fixed_vars <- all.vars(formula(model))[-1]
  response <- all.vars(formula(model))[1]
  
  # model data
  data <- model.frame(model)
  
  # log transform if the response is biomass
  if (family(model)$link == "log") {
    data[[response]] <- log(data[[response]])
  }
  
  # loop through each variable in the fixed formula
  for (var in vars) {
    
    # create new data: vary current var, others fixed at their mean
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
    
    # add year column (required)
    nd$year <- 2014
    
    # predict
    pred <- predict(
      model,
      newdata = nd,
      se_fit = TRUE,
      re_form = NA
      )
    
    # compute confidence intervals
    pred$lwr <- pred$est - 1.96 * pred$est_se
    pred$upr <- pred$est + 1.96 * pred$est_se
    
    # create ggplot
    p <- ggplot() +
      geom_point(
        data = data,
        aes(x = .data[[var]], y = .data[[response]]),
        pch = 19,
        fill = "gray40",
        alpha = 0.1
        ) +
      geom_line(
        data = pred, 
        aes(x = .data[[var]], y = est),
        color = "darkred",
        linewidth = 0.5
        ) +
      geom_ribbon(
        data = pred,
        aes(x = .data[[var]], 
            ymin = lwr, 
            ymax = upr),
        alpha = 0.4,
        fill = "darkred"
      ) + 
      labs(title = paste(var, "effect on",  response, " ± 95% CI"),
           x = var, y = response) +
      theme_bw() +
      theme(axis.title = element_text(size = 14))
    
    # save plot
    ggsave(paste0(
      directory,
      "/", name,
      "_prediction_",
      var,
      ".png"),
      plot = p, width = 9, height = 9, units = "cm")
  }
}

# Save
save(model_prediction, file = "tools/model_prediction_function.RData")

## Grid plotting function (TO BE FIXED) ----

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

rm(list = ls())

## End