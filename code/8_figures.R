## File: 8_figures.R
## Purpose: to produce plots that compare models
## Author: Filippomaria Cassarino
## Date: 
## ------------------------------------------------------------------------ ----
## Notes ----

## Library ---- 

# Import data, objects and functions
load("tools/install_load_function.RData")

# Load required packages
required_pakages <- c("dplyr",
                      "ggplot2",
                      "tidyr") 

install_load_function(required_pakages)

# Objects to keep when cleaning
keep <- c("install_load_function")

# Cleaning
rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Biomass components ----

# Import data
b_data <- read.csv("data/final/final_data.csv")

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

plot_data <- b_data %>%
  group_by(year) %>%
  summarise(
    across(all_of(sp), sum, na.rm = TRUE),
    total_biomass = sum(total_biomass, na.rm = TRUE)
  ) %>%
  mutate(
    Others = total_biomass - rowSums(across(all_of(sp)))
  ) %>%
  select(year, all_of(sp), Others) %>%
  pivot_longer(
    -year,
    names_to = "species",
    values_to = "biomass"
  )

plot_data$species <- factor(plot_data$species, levels = sp)

levels(plot_data$species) <- c(
  "G. morhua",
  "M. aeglefinus",
  "S. mentella",
  "M. poutassou",
  "H. platessoides",
  "M. villosus",
  "P. virens",
  "B. saida",
  "T. esmarkii",
  "R. hippoglossoides",
  "Others")

cols <- c(colorRampPalette(c("#2C3E50","#03A9F4", "#E2F4C7"))(10),
          "gray50")

p <- ggplot(
  plot_data,
  aes(
    x = year,
    y = biomass,
    fill = species
  )
) +
  geom_area(position = "fill") +
  geom_line(position = "fill", color = "gray30", linewidth = 0.1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Year",
    y = "Proportion of total biomass",
    fill = "Species"
  ) +
  theme_minimal(base_size = 14)

ggsave(p, file = "figures/biomass_plots/top_species_biomass.png",
       unit = "cm", height = 16, width = 25)

rm(plot_data)

## General and fishing model ----

# Import csv files of model coefficients

main_models <- read.csv2( # all biomass (general) model
  file = "models/biomass/general_model/general_biomass_model_coefficients.csv",
  sep = " ") 

f <- read.csv2(
  file = "models/biomass/fishing_model/fishing_biomass_model_coefficients.csv",
 sep = " ")

main_models <- bind_rows(main_models, f)
  
rm(f)


# Data to plot
plot_data <- main_models %>%
  filter(term != "(Intercept)") %>%
  drop_na() %>%
  mutate(across(-c(model, term), as.numeric))

plot_data$model <- factor(plot_data$model, 
                          levels = c("general_biomass", 
                                     "fishing_biomass"))

levels(plot_data$model) <- c("General model (2004-2022)",
                             "With fishing (2012-2022)")

plot_data$term <- factor(plot_data$term, 
                         levels = c("depth", "log_chla_mean",
                                    "sst_sd","sbt_sd",
                                    "sic_mean","sbt_mean",
                                    "Fric","Feve",
                                    "Tric","Teve","log_fishing_effort"))

# Plot

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.15,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("#000000", "#D55E00")) +
  coord_flip() +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     aspect.ratio = 2,
    panel.grid.minor = element_blank(),
    #legend.text = element_text(size = 10),
    #legend.title = element_text(size = 11),
    #legend.key.size = unit(0.4, "cm")
  ) +
  labs(
    x = "Model term",
    y = "\nCoefficient ± 1 SE",
    color = "Model")

ggsave(p, file = "figures/model_figures/general_fishing_models.png",
       unit = "cm", height = 20, width = 25)

rm(plot_data)

## Area comparison ----

# Import csv files of model coefficients
csv_name <- paste0(c("arctic_boreal_model/arctic",
                     "arctic_boreal_model/boreal"), "_biomass_model_coefficients.csv")

area_effect_data <- read.csv2( # all biomass (general) model
  file = "models/biomass/general_model/general_biomass_model_coefficients.csv",
  sep = " ") 

for (i in 1:length(csv_name)) { # bind all the others
  data <- read.csv2(file = paste0("models/biomass/",
                                  csv_name[i]), sep = " ")
  area_effect_data <- bind_rows(area_effect_data, data)
  
  rm(data)
}

# Data to plot
plot_data <- area_effect_data %>%
  filter(term != "(Intercept)") %>%
  drop_na() %>%
  mutate(across(-c(model, term), as.numeric))

plot_data$model <- factor(plot_data$model, 
                          levels = c("general_biomass", 
                                     "boreal_biomass",
                                     "arctic_biomass"))

levels(plot_data$model) <- c("All areas",
                             "Boreal areas",
                             "Arctic areas")

plot_data$term <- factor(plot_data$term, 
                         levels = c("depth", "log_chla_mean",
                                    "sst_sd","sbt_sd",
                                    "sic_mean","sbt_mean",
                                    "Fric","Feve",
                                    "Tric","Teve", "log_fishing_effort"))

# Plot

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.15,
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
    x = "Model term",
    y = "\nCoefficient ± 1 SE",
    color = "Model")

ggsave(p, file = "figures/model_figures/area_effect.png",
       unit = "cm", height = 20, width = 25)

rm(plot_data)

## Time period comparison ----

# Import csv files of model coefficients
csv_name <- paste0(c("2004_2011",
                     "2012_2016",
                     "2017_2022"), "_biomass_model_coefficients.csv")

time_effect_data <- read.csv2( # all biomass (general) model
  file = "models/biomass/general_model/general_biomass_model_coefficients.csv",
  sep = " ") 

for (i in 1:length(csv_name)) { # bind all the others
  data <- read.csv2(file = paste0("models/biomass/time_period_model/",
                                  csv_name[i]), sep = " ")
  time_effect_data <- bind_rows(time_effect_data, data)
  
  rm(data)
}

# Data to plot
plot_data <- time_effect_data %>%
  filter(term != "(Intercept)") %>%
  drop_na() %>%
  mutate(across(-c(model, term), as.numeric))

plot_data$model <- factor(plot_data$model, 
                          levels = c("general_biomass",
                                     "2004_2011_biomass",
                                     "2012_2016_biomass",
                                     "2017_2022_biomass"))

levels(plot_data$model) <- c("All time series",
                             "2004 to 2011 (coldest)",
                             "2012 to 2016 (warmest)",
                             "2017 to 2022 (cooler)")

plot_data$term <- factor(plot_data$term, 
                         levels = c("depth", "log_chla_mean",
                                    "sst_sd","sbt_sd",
                                    "sic_mean","sbt_mean",
                                    "Fric","Feve",
                                    "Tric","Teve", "log_fishing effort"))

# Plot
cols <- colorRampPalette(c("darkblue", "lightblue"))(3)

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.15,
                position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("black", cols)) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .7),
                     breaks = seq(-.5, .5, by = .25)) +
  theme_minimal(base_size = 14) +
  theme(     aspect.ratio = 2,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term",
    y = "\nCoefficient ± 1 SE",
    color = "Model")

ggsave(p, file = "figures/model_figures/time_period_effect.png",
       unit = "cm", height = 20, width = 25)

rm(plot_data)

## Partial biomass model coefficients ----

# Import csv files of model coefficients
csv_name <- paste0("remove_", 1:10, "_biomass_model_coefficients.csv")

species_effect_data <- read.csv2( # all biomass (general) model
  file = "models/biomass/general_model/general_biomass_model_coefficients.csv",
  sep = " ") 

for (i in 1:length(csv_name)) { # bind all the others
  data <- read.csv2(file = paste0("models/biomass/species_removal_model/",
                                  csv_name[i]), sep = " ")
  species_effect_data <- bind_rows(species_effect_data, data)
  
  rm(data)
}

# Data to plot
plot_data <- species_effect_data %>%
  filter(term != "(Intercept)") %>%
  drop_na() %>%
  mutate(across(-c(model, term), as.numeric))

plot_data$model <- factor(plot_data$model, 
                          levels = c("general_biomass", paste0("remove_", 1:10, "_biomass")))

levels(plot_data$model) <- c("All biomass", paste0("Top ", 1:10, " species removed"))

plot_data$term <- factor(plot_data$term, 
                         levels = c("depth", "log_chla_mean",
                                    "sst_sd","sbt_sd",
                                    "sic_mean","sbt_mean",
                                    "Fric","Feve",
                                    "Tric","Teve", "log_fishing_effort"))

# Plot
cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 2,
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = estimate - std.error,
                    ymax = estimate + std.error),
                width = 0.15,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = cols) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(-.7, .87),
                     breaks = c(seq(-.5, .5, by = .25), .75)) +
  theme_minimal(base_size = 14) +
  theme(     aspect.ratio = 2,
    panel.grid.minor = element_blank()) +
  labs(
    x = "Model term",
    y = "\nCoefficient ± 1 SE",
    color = "Model")

ggsave(p, file = "figures/model_figures/top_species_removal_effect.png",
       unit = "cm", height = 20, width = 25) 

rm(plot_data)

## Diversity models' coefficients ----

vars <- c("Teve", "Tric", "Tdiv", "Feve", "Fric", "Fdis")


for (i in 1:length(vars)) {
  
  # Import csv files of model coefficients
  
  main_models <- read.csv2( 
    file = paste0("models/", vars[i], "/general_",
    vars[i], "_model_coefficients.csv"),
    sep = " ")
  
  f <- read.csv2( 
    file = paste0("models/", vars[i], "/fishing_",
    vars[i], "_model_coefficients.csv"),
    sep = " ") 
  
  main_models <- bind_rows(main_models, f)
  
  rm(f)
  
  
  # Data to plot
  plot_data <- main_models %>%
    filter(term != "(Intercept)") %>%
    drop_na() %>%
    mutate(across(-c(model, term), as.numeric))
  
  plot_data$model <- factor(plot_data$model, 
                            levels = c(paste0("general_", vars[i]), 
                                       paste0("fishing_", vars[i])))
  
  levels(plot_data$model) <- c("General model (2004-2022)",
                               "With fishing (2012-2022)")
  
  plot_data$term <- factor(plot_data$term, 
                           levels = c("depth", "log_chla_mean",
                                      "sst_sd","sbt_sd",
                                      "sic_mean","sbt_mean","log_fishing_effort"))
  
  # Plot
  
  p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
    geom_point(size = 2,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = estimate - std.error,
                      ymax = estimate + std.error),
                  width = 0.15,
                  position = position_dodge(width = 0.3)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("#000000", "#D55E00")) +
    coord_flip() +
    scale_y_continuous(limits = c(-.7, .7),
                       breaks = seq(-.5, .5, by = .25)) +
    theme_minimal(base_size = 14) +
    theme(     aspect.ratio = 2,
      panel.grid.minor = element_blank()) +
    labs(
      x = "Model term",
      y = "\nCoefficient ± 1 SE",
      color = "Model")
  
  ggsave(p, file = paste0("figures/model_figures/a_", vars[i], "_models.png"),
         unit = "cm", height = 20, width = 25)
  
}

rm(plot_data)

## End