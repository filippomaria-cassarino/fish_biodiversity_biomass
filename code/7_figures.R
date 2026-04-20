## File: 8_figures.R
## Purpose: to produce plots that compare models
## Author: Filippomaria Cassarino
## Date: 
## Notes ----

## Fix read.csv
## Library ---- 

# Import data, objects and functions
load("tools/install_load_packages_function.RData")

# Load required packages
required_pakages <- c("dplyr",
                      "ggplot2",
                      "tidyr") 

install_load_packages(required_pakages)

keep <- c(
  "biomass_coeff",
  "diversity_coeff",
  "keep"
  )

## Prepare biomass coefficient data ----

# Read biomass coefficient data
biomass_coeff <- read.csv(
  file = "models/biomass/biomass_models_coefficients.csv",
  sep = ",")

remove <- paste0("remove_", 1:10, "_biomass")

# Identify levels
biomass_coeff$model <- factor(biomass_coeff$model, 
                               levels = c(
                                 "general_biomass", 
                                 "fishing_biomass",      
                                 "arctic_biomass",
                                 "boreal_biomass",   
                                 "2004_2011_biomass",
                                 "2012_2016_biomass",
                                 "2017_2022_biomass",
                                 remove,
                                 "pielou_biomass",
                                 paste0("pielou_", remove),
                                 "other_fd_biomass"
                               ))

# Rename levels
levels(biomass_coeff$model) <- c(
  "General model",
  "Model with fishing",
  "Arctic model",
  "Boreal model",
  "2004 to 2011 (coldest)",
  "2012 to 2016 (warmest)",
  "2017 to 2022 (cooler)",
  paste0("Top ", 1:10, " removed"),
  "Pielou model", 
  paste0("Top ", 1:10, " removed (Pielou)"),
  "Other FD method"
)

# Identify levels
biomass_coeff$term <- factor(biomass_coeff$term, levels = c(
  "depth", 
  "log_chla_mean",
  "sst_sd",
  "sbt_sd",
  "sic_mean",
  "sbt_mean",
  "kde_fric",
  "kde_feve",
  "ch_fric",
  "ch_feve",
  "tric",
  "teve",
  "p_teve",
  "log_sum_fishing"
))

# Rename levels
levels(biomass_coeff$term) <- c(
  "Depth", 
  "Chla mean",
  "SST SD",
  "SBT SD",
  "SIC mean",
  "SBT mean",
  "Fric",
  "Feve",
  "Fric",
  "Feve",
  "Tric",
  "Teve",
  "Teve",
  "Fish. eff."
)

## General and fishing comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("General model",
                      "Model with fishing")) 

# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#000000",
                                "#D55E00")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 2,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± 95% CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_general_fishing_models.png",
       unit = "cm", height = 20, width = 25)

rm(list = setdiff(ls(), keep))

## Biogeographical area comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("General model",
                      "Arctic model",
                      "Boreal model")) 

# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#000000",
                                "#E69F00",
                                "#56B4E9")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 2,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± 95% CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_area_effect.png",
       unit = "cm", height = 20, width = 25)

rm(list = setdiff(ls(), keep))

## Time period comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("General model",
                      "2004 to 2011 (coldest)",
                      "2012 to 2016 (warmest)",
                      "2017 to 2022 (cooler)")) 

# Plot
cols <- colorRampPalette(c("darkblue", "lightblue"))(3)

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#000000", cols)) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 2,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± 95% CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_period_effect.png",
       unit = "cm", height = 20, width = 25)

rm(list = setdiff(ls(), keep))

## Species removal comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("General model",
                      paste0("Top ", 1:10, " removed"))) 

# Plot
cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = cols) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, 1),
                     breaks = seq(-.5, .75, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 1.75,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± 95% CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_species_effect.png",
       unit = "cm", height = 20, width = 26) 

rm(list = setdiff(ls(), keep))

## Pielou's evenness comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("Pielou model",
                      paste0("Top ", 1:10, " removed (Pielou)"))) 

# Plot
cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = cols) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, 1),
                     breaks = seq(-.5, .75, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 1.75,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± 95% CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_species_effect_pielou.png",
       unit = "cm", height = 20, width = 26) 

rm(list = setdiff(ls(), keep))

## Other functional diversity comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("General model",
                      "Other FD method")) 
# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#000000",
                                "forestgreen")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     
    aspect.ratio = 2,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± 95% CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/a_different_fd_effect.png",
       unit = "cm", height = 20, width = 25)

rm(list = setdiff(ls(), keep))

## Prepare diversity coefficient data ----

# Read biomass coefficient data
diversity_coeff <- read.csv(
  file = "models/diversity/diversity_models_coefficients.csv",
  sep = ",")

# Identify levels
diversity_coeff$model <- factor(diversity_coeff$model, 
                              levels = c(
                                "general_teve",
                                         "general_tric",
                                         "general_kde_feve",
                                         "general_kde_fric",
                                         "fishing_teve",
                                         "fishing_tric",
                                         "fishing_kde_feve",
                                         "fishing_kde_fric"
                                         ))

# Rename levels
levels(diversity_coeff$model) <- c(
  "Teve model",
  "Tric model",
  "Feve model",
  "Fric model",
  "Teve model (fishing)",
  "Tric model (fishing)",
  "Feve model (fishing)",
  "Fric model (fishing)"
)

# Identify levels
diversity_coeff$term <- factor(diversity_coeff$term, levels = c(
  "depth", 
  "log_chla_mean",
  "sst_sd","sbt_sd",
  "sic_mean",
  "sbt_mean",
  "log_sum_fishing"
))

# Rename levels
levels(diversity_coeff$term) <- c(
  "Depth", 
  "Chla mean",
  "SST SD",
  "SBT SD",
  "SIC mean",
  "SBT mean",
  "Fish. eff."
)

## Diversity models (no fishing) ----

# Divide data based on fishing
no_f_data <- diversity_coeff %>% filter(!grepl("(fishing)", model))

f_data <- diversity_coeff %>% filter(grepl("(fishing)", model))
       
# Plot for models without fishing
p <- ggplot(data = no_f_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#780000", "#c1121f", "#003049", "#669bbc")) +
  coord_flip() +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     aspect.ratio = 2,
             panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/diversity_models_no_fishing.png",
       unit = "cm", height = 20, width = 25)

# Plot for models with fishing
p <- ggplot(data = f_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#780000", "#c1121f", "#003049", "#669bbc")) +
  coord_flip() +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     aspect.ratio = 2,
             panel.grid.minor = element_blank()) +
  labs(
    x = "Model term\n",
    y = "\nCoefficient ± CI",
    color = "Model")

ggsave(p, file = "figures/model_figures/diversity_models_with_fishing.png",
       unit = "cm", height = 20, width = 25)

rm(list = setdiff(ls(), keep))

## End