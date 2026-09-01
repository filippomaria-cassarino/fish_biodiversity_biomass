## File: 8_figures.R
## Purpose: to produce plots that compare models
## Author: Filippomaria Cassarino
## Date: 1 Sep 2026
## Notes ----

## Library ---- 

# Import data, objects and functions
load("tools/install_load_packages_function.RData")

# Install/load required packages
install_load_packages(c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "forcats", # na to level
  "ggspatial", # annotation_scale
  "sf",
  "terra",
  "rnaturalearth", # land shape
  "ggrepel", # add species names to the plots
  "ggcorrplot", # correlation plot
  "patchwork" # multiple plot figures
))

keep <- c(
  "biomass_coeff",
  "diversity_coeff",
  "keep",
  "coeff_theme",
  "map_theme"
  )

coeff_theme <- function(x) {
  theme(     
  aspect.ratio = 1.8,
  plot.margin = margin(0, 0, 0, 0),
  legend.position = "right",
  legend.justification = "top",
  legend.margin = margin(0, 0, 0, 0),
  legend.key.width = unit(0.2, "cm"),
  legend.text = element_text(margin = margin(l = 2, r = 2),
                             size = 10),
  legend.title = element_text(size = 12),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 10),
  axis.ticks = element_blank(),
  title = element_text(size = 13),
  panel.background = element_rect(fill = "white", color = NA),
  panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.3),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "gray92", linewidth = 0.3)
) 
}

map_theme <- function(x) {
  theme(
    plot.title = element_text(size = 13, margin = margin(b = 2)),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "bottom",
    legend.justification = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.8)
  ) 
}

## ------------------------------------------------------------------------ ----
## Study area maps ----

# depth data
depth_rast <- rast("data/original/depth_gebco_400.tiff")

# simplify depth layer and transform to data frame
d <- depth_rast %>%
  aggregate(fact = 20, fun = mean, na.rm = TRUE) %>% # EPSG:3996
  as.data.frame(xy = TRUE) %>%
  mutate(depth = ifelse(depth_gebco_400 < 0, 1, 0))

# plot of the arctic
p <- ggplot() +
  geom_raster(data = d, aes(x = x, y = y, fill = depth)) +
  scale_fill_gradientn(colors = c("gray45", "white"),
                       name = "Depth (m)") +
  labs(title = "Barents Sea: Bathymetry and Landmass",
       x = "Longitude", y = "Latitude") +
  guides(fill = guide_colorbar(barwidth = .5,
                               barheight = 5)) +
  theme_void()

ggsave(p, file = paste0("figures/maps/large_map.jpg"),
       width = 17, height = 15, units = "cm")

# common crs
crs_aeqd <- sf::st_crs(
  "+proj=aeqd +lat_0=76 +lon_0=27 +datum=WGS84 +units=km +no_defs"
  )

# simplify depth layer and transform to data frame
d <- depth_rast %>%
  aggregate(fact = 10, fun = mean, na.rm = TRUE) %>%
  terra::project(crs_aeqd$wkt) %>%
  as.data.frame(xy = TRUE) %>%
  mutate(
    depth = ifelse(depth_gebco_400 < 0, depth_gebco_400, NA),
    depth_cat = cut(
      abs(depth),
      breaks = c(0, 100, 300, 500, 1000, 2000, Inf),
      labels = c("0–100", "100–300", "300–500", "500–1000", "1000-2000", "2000+"),
      right = FALSE
    ),
    depth_cat = fct_na_value_to_level(depth_cat, level = "Land")
  )

# eez data
eez <- sf::st_read("data/original/EEZ/EEZ_land_union_v4_202410.shp") %>%
  filter(UNION %in% c("Norway", "Svalbard", "Russia", "Greenland")) %>%
  st_transform(crs_aeqd)

# study area plot
l <- st_crop(ne_countries(scale = "medium", returnclass = "sf"),
             xmin = -20, xmax = 80, # range for the study area (sf)
             ymin = 50, ymax = 90) |>
  st_transform(crs = crs_aeqd)

# community data
data <- read.csv("data/final/final_data.csv", sep = ",")

hauls <- data %>%
  sf::st_as_sf(
    coords = c("longitude", "latitude"), 
    crs = crs_aeqd) 

# plot
p <- ggplot() +
  geom_raster(data = d, aes(x = x, y = y, fill = depth_cat)) +
  geom_sf(data = eez, fill = NA, color = "white",
          linetype = "solid",
          linewidth = .3,
          size = .1) +
  geom_sf(data = hauls, size = .2) +
  scale_fill_manual(
    name = "Depth (m)",
    values = c(
      "0–100" = "#c6dbef",
      "100–300" = "#9ecae1",
      "300–500" = "#6baed6",
      "500–1000" = "#3182bd",
      "1000-2000" = "#08519c",
      "2000+" = "#06396d",
      "Land" = "gray45"
    )) +
  coord_sf(xlim = c(-600, 1200), ylim = c(-900, 800), expand = FALSE) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.7
  ) +
  theme_void() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.box.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.3),
    legend.box.margin = margin(2, 2, 2, 2))

ggsave(p, file = paste0("figures/maps/study_area.jpg"),
       width = 10, units = "cm", dpi = 1200)

## Boreal - Arctic distinction 
divide_bor <- quantile(data$sbt_mean, 0.7, na.rm = TRUE)
divide_arc <- quantile(data$sbt_mean, 0.3, na.rm = TRUE)

area_data <- hauls %>%
  mutate(
  area = case_when(
    sbt_mean <= divide_arc ~ "Arctic",
    sbt_mean >= divide_bor ~ "Boreal",
    sbt_mean > divide_arc & sbt_mean < divide_bor ~ "Excluded"
  )) %>%
  drop_na(area)

p <- ggplot() +
  geom_raster(data = d, aes(x = x, y = y, fill = depth_cat), alpha = .3) +
  geom_sf(data = area_data, aes(color = area), size = .5) +
  scale_fill_manual(
    name = "Depth (m)",
    values = c(
      "0–100" = "#c6dbef",
      "100–300" = "#9ecae1",
      "300–500" = "#6baed6",
      "500–1000" = "#3182bd",
      "1000-2000" = "#08519c",
      "2000+" = "#06396d",
      "Land" = "gray45"
    )) +
  scale_color_manual(
    name = "Area",
    values = c(
      "Arctic" = "#8EB1C7",
      "Boreal" = "#D5465E",
      "Excluded" = "black"
    )) +
  coord_sf(xlim = c(-600, 1200), ylim = c(-900, 800), expand = FALSE) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.7
  ) +
  theme_void() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.box.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.3),
    legend.box.margin = margin(2, 2, 2, 2))

ggsave(p, file = paste0("figures/maps/boreal_vs_arctic_area.jpg"),
       width = 12, height = 12, units = "cm")

## ------------------------------------------------------------------------ ----
## Biomass composition ----

# Load community data and top species
load("data/intermediate/community_filtered.RData")
load("data/intermediate/top_10_species.RData")

# Plot data
plot_data <- community_filtered %>%
  group_by(year, taxon) %>%
  summarise(biomass = sum(wgt_cpua, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    species = if_else(taxon %in% top_10, taxon, "Others")
  ) %>%
  group_by(year, species) %>%
  summarise(biomass = sum(biomass), .groups = "drop") 

plot_data$species <- factor(plot_data$species, levels = c(top_10, "Others"))

proportions <- plot_data %>% 
  group_by(species) %>%
  summarise(prop = round(sum(biomass)/sum(plot_data$biomass) * 100)) %>%
  left_join(data_frame(species = c(top_10, "Others")), by = "species")

levels(plot_data$species) <- paste0(
  proportions$prop,
  "% ",
  c(top_10, "Others"))

cols <- c(colorRampPalette(c("#2C3E50","#03A893", "#E2F4C7"))(10),
          "gray50")

p <- ggplot(plot_data, aes(x = year, y = biomass, fill = species)) +
  geom_area(position = "fill") +
  geom_line(position = "fill", color = "gray30", linewidth = 0.1) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(
    data = data.frame(year = seq(2004, 2022, by = 2)),
    aes(x = year, y = -0.05, label = year),
    angle = 90,
    vjust = 0.5,
    size = 2.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols) +
  labs(
    x = NULL,
    y = NULL,
    #fill = "Proportion of surveyed biomass"
    fill = NULL
  ) +
  geom_hline(
    yintercept = 0.5,
    color = "white",
    linetype = "dotted",
    linewidth = 0.2
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.margin = margin(1, 1, 1, 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length = unit(-0.40, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8, face = "italic"),
    legend.position = "right",
    legend.justification = "top",
    legend.key.width = unit(0.2, "cm"),
    legend.key.height  = unit(0.509, "cm"),
    legend.spacing.x = unit(0, "cm"),
    legend.box.margin = margin(6.5, 0, 0, -10)
  )

# Save plot
ggsave(p, file = "figures/figure_2.jpg",
       unit = "cm", height = 6.5, width = 8.5, 
       dpi = 1200)

## ------------------------------------------------------------------------ ----
## Functional space plot ----
load("data/intermediate/functional_diversity_matrixes.RData")
load("data/intermediate/top_10_species.RData")

# All species vs top 10 species trait combinations
all_species <- trait_matrixes[["total_biomass"]] %>%
  mutate(taxon = row.names(.)) 

top_10_species <- all_species %>%
  filter(taxon %in% top_10)

# loop to plot all dimension comparisons
axes <- combn(1:4, 2, simplify = FALSE)

for (ax in axes) {
  
  xvar <- paste0("Axis.", ax[1])
  yvar <- paste0("Axis.", ax[2])
  
  p <- ggplot(all_species, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(color = "Others"), alpha = 0.5,
               size = .3) +
    geom_point(
      data = top_10_species,
      aes(color = "Top 10"),
      size = .5
    ) +
    geom_text_repel(
      data = top_10_species,
      aes(label = taxon),
      size = 1.5,
      segment.size = .2
    ) +
    scale_color_manual(
      values = c(
        "Other" = "grey70",
        "Top 10" = "black"
      ),
      name = NULL
    ) +
    theme_bw(base_size = 8) +
    theme(legend.position = "none")
  
  # save
  ggsave(p,
         file = paste0("figures/functional_space_", ax[1], "_", ax[2], ".jpg"),
         width = 5, height = 5, units = "cm")
}

## ------------------------------------------------------------------------ ----
## Correlation plot ----

# Import data
final_data <- read.csv("data/final/final_data.csv")

# Prepare data of the whole study area (all) 
cor_data <- final_data %>%
  select(
    kde_fric,
    kde_feve,
    tric,
    teve,
    sst_mean,
    sbt_mean,
    sic_mean,
    log_chla_mean,
    sbt_sd,
    sst_sd,
    sic_sd,
    log_chla_sd,
    depth
  ) %>%
  rename(
    "Depth" = depth,
    "log chla mean" = log_chla_mean,
    "log chla SD" = log_chla_sd,
    "SST SD" = sst_sd,
    "SBT SD" = sbt_sd,
    "SIC SD" = sic_sd,
    "SST mean" = sst_mean,
    "SBT mean" = sbt_mean,
    "SIC mean" = sic_mean,
    "Fric" = kde_fric,
    "Feve" = kde_feve,
    "Tric" = tric,
    "Teve" = teve
  ) %>%
  drop_na() %>%
  mutate(across(everything(), scale)) %>%
  cor()

# Plot
p <- ggcorrplot(
  cor_data,
  type = "upper",
  hc.order = FALSE,
  lab = TRUE
  ) +
  scale_fill_gradientn(
    colours = c(
      "#003f5c",
      "lightblue",
      "gray85",
      "#e48646",
      "darkred"
    ))

ggsave(p, file = "figures/correlation_plot.jpg",
       unit = "cm", height = 17, width = 17)

## ------------------------------------------------------------------------ ----
## Prepare biomass coefficient data and save table ----

# Read biomass coefficient data
biomass_coeff <- read.csv(
  file = "models/biomass/biomass_models_coefficients.csv",
  sep = ",") 

new_res <- paste0("remove_", 1:10)

# Identify levels
biomass_coeff$model <- factor(biomass_coeff$model, 
                               levels = c(
                                 "general_biomass",
                                 paste0("only_", 
                                        c("teve",
                                          "kde_feve",
                                          "tric",
                                          "kde_fric"), "_biomass"),
                                 "fishing_biomass",
                                 "boreal_biomass",   
                                 "arctic_biomass",
                                 "2004_2011_biomass",
                                 "2012_2016_biomass",
                                 "2017_2022_biomass",
                                 paste0(new_res, "_biomass"),
                                 paste0(new_res, "_samediv_biomass"),
                                 "pielou_biomass",
                                 "other_fd_biomass"
                               ))

# Rename levels
levels(biomass_coeff$model) <- c(
  "Main",
  paste0("Only ", c("Teve",  "Feve", "Tric", "Fric")),
  "With fishing",
  "Boreal",
  "Arctic",
  "2004-2011",
  "2012-2016",
  "2017-2022",
  paste0("-Top ", 1:10),
  paste0("-Top ", 1:10, " (Same diversty)"),
  "Pielou", 
  "CH FD"
)

# Rename levels
biomass_coeff <- biomass_coeff %>%
  mutate(
    term = case_when(
      grepl("tric", term) ~ "tric",
      grepl("teve", term) ~ "teve",
      grepl("fric", term) ~ "fric",
      grepl("feve", term) ~ "feve",
      TRUE ~ term),
    # Add significance column
    significance = case_when(
      conf.low * conf.high > 0 ~ "Significant",
      conf.low * conf.high <= 0 ~ "Not significant"
    )) %>%
  select(-X)


# Identify levels
biomass_coeff$term <- factor(biomass_coeff$term, levels = c(
  "depth", 
  "log_chla_mean",
  "sst_sd",
  "sbt_sd",
  "sic_mean",
  "sbt_mean",
  "fric",
  "feve",
  "tric",
  "teve",
  "log_sum_fishing"
))

# Rename levels
levels(biomass_coeff$term) <- c(
  "Depth", 
  "log chla mean",
  "SST SD",
  "SBT SD",
  "SIC mean",
  "SBT mean",
  "Fric",
  "Feve",
  "Tric",
  "Teve",
  "log fish. eff."
)

# Save table
write.table(
  biomass_coeff,
  file = "tables/biomass_models_table.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)

## ------------------------------------------------------------------------ ----
## General and fishing comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c(
   "Main",
   "With fishing"
    )) 

# Plot
pa <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_color_manual(values = c("#000000",
                                "#7BB086")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
     limits = c(-.5, .5),
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "A)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme()

ggsave(pa, file = "figures/model_figures/A_general_fishing.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

rm(list = setdiff(ls(), c(keep, "pa")))

## Time period comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("Main",
                      "2004-2011",
                      "2012-2016",
                      "2017-2022")) 

# Plot
cols <- colorRampPalette(c("#334d80", "#cbcaa5"))(3)

pb <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_color_manual(values = c("#000000", cols)) +
  coord_flip(ylim = c(-.5, .5)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "B)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme()

ggsave(pb, file = "figures/model_figures/B_period.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

rm(list = setdiff(ls(), c(keep, "pa", "pb")))

## Biogeographical area comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("Main",
                      "Boreal",
                      "Arctic")) 

# Plot
pc <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_color_manual(values = c("#000000",
                                "#D5465E",
                                "#8EB1C7")) +
  coord_flip(ylim = c(-.5, .5)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "C)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme()

ggsave(pc, file = "figures/model_figures/C_area.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

rm(list = setdiff(ls(), c(keep, "pa", "pb", "pc")))

## Partial biomass comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("Main", paste0("-Top ", 1:10))) 

# Plot
cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

pd <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.8, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.8, reverse = TRUE)) +
  scale_color_manual(values = cols) +
  coord_flip(ylim = c(-.5, .75)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = seq(-.5, .75, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5", "0.75"),
    expand = c(0, 0)
  ) +
  labs(
    title = "D)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme() +
  theme(aspect.ratio = 1.45)

ggsave(pd, file = "figures/model_figures/D_partial_biomass.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

rm(list = setdiff(ls(), c(keep, "pa", "pb", "pc", "pd")))

## Figure 2 - main model comparisons ----

# Adjust themes
pA <- pa + 
  theme(
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

pB <- pb + 
  theme(
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

pC <- pc + 
  theme(
    axis.title.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0))

pD <- pd +
  theme(
    aspect.ratio = 1.45,
    plot.margin = margin(0, 0, 0, 0),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Space the plots
row1 <- pA + plot_spacer() + pB + plot_spacer() +
  plot_layout(widths = c(0.84, 0, 0.84, 0.1))

row2 <- pC + plot_spacer() + pD + plot_spacer() +
  plot_layout(widths = c(0.74, 0.15, 0.92, 0))

final_plot <- row1 / row2 +
  plot_annotation(
    theme = theme(plot.margin = margin(0, 0, 0, 0)))

# Save plot
ggsave(
  "figures/figure_3.jpg",
  plot = final_plot,
  width = 17,
  height = 19,
  units = "cm",
  dpi = 1200
)

## Pielou's evenness comparison ----

# Data to plot
plot_data <- biomass_coeff %>%
  filter(model %in% c("Main", "Pielou")) 

# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_color_manual(values = c("black", "gray")) +
  coord_flip(ylim = c(-.75, .5)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = seq(-.75, .5, by = .25),
    labels = c("-0.75", "-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme() +
  theme(aspect.ratio = 1.45)

ggsave(p, file = "figures/model_figures/supplementary_pielou.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

rm(list = setdiff(ls(), keep))

## Other models ----

## Partial biomass models using the same diversity metrics of the main model
plot_data <- biomass_coeff %>%
  filter(model %in% c("Main", paste0("-Top ", 1:10, " (Same diversty)"))) 

cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.8, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.8, reverse = TRUE)) +
  scale_color_manual(values = cols) +
  coord_flip(ylim = c(-.5, .75)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = seq(-.5, .75, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5", "0.75"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Supplemnetary X",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme() +
  theme(aspect.ratio = 1.45)

ggsave(p, file = "figures/model_figures/Supplementary_X_partial_biomass.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

## Main models using convex-hull estimated FD
plot_data <- biomass_coeff %>%
  filter(model %in% c("Main",
                      "CH FD")) 
# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.5, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.5, reverse = TRUE)) +
  scale_color_manual(values = c("#000000",
                                "#e67090")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
     limits = c(-.5, .5),
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Supplementary X",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme()

ggsave(p, file = "figures/model_figures/supplementary_other_fd.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

## Main models with only one diversity metric
plot_data <- biomass_coeff %>%
  filter(model %in% c(
    "Main", paste0("Only ", c("Teve", "Feve",  "Tric", "Fric")))) 

# Plot
p <- ggplot(data = plot_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6, reverse = TRUE)) +
  scale_color_manual(values = c("#000000",
                                "#326747", "#98CDAD",
                                "#67325F", "#CD98C5")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(-.5, .5),
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Supplemtary X",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme()

ggsave(p, file = "figures/model_figures/Supplementary_X_general_nofd_notd.jpg",
       unit = "cm", height = 9, width = 13, dpi = 1200)

rm(list = setdiff(ls(), keep))

## ------------------------------------------------------------------------ ----
## Prepare diversity coefficient data and save table ----

# Read biomass coefficient data
diversity_coeff <- read.csv(
  file = "models/diversity/diversity_models_coefficients.csv",
  sep = ",")

# Identify levels
diversity_coeff$model <- factor(diversity_coeff$model, 
                                levels = c(
                                  "general_teve",
                                  "general_kde_feve",
                                  "general_tric",
                                  "general_kde_fric",
                                  "fishing_teve",
                                  "fishing_kde_feve",
                                  "fishing_tric",
                                  "fishing_kde_fric"
                                ))

# Rename levels
levels(diversity_coeff$model) <- c(
  "Teve",
  "Feve",
  "Tric",
  "Fric",
  "Teve     (fishing)",
  "Feve     (fishing)",
  "Tric     (fishing)",
  "Fric     (fishing)"
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
  "log chla mean",
  "SST SD",
  "SBT SD",
  "SIC mean",
  "SBT mean",
  "log fish. eff."
)

diversity_coeff <- diversity_coeff %>%
mutate(
  significance = case_when(
    conf.low * conf.high > 0 ~ "Significant",
    conf.low * conf.high <= 0 ~ "Not significant"
  )) %>%
  select(-X)

# Save table
write.table(
  diversity_coeff,
  file = "tables/diversity_models_table.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)

## ------------------------------------------------------------------------ ----
## Diversity models ----

# Divide data based on fishing
no_f_data <- diversity_coeff %>% filter(!grepl("(fishing)", model))

f_data <- diversity_coeff %>% filter(grepl("(fishing)", model))
       
# Plot for models without fishing
p <- ggplot(data = no_f_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6, reverse = TRUE)) +
  
  scale_color_manual(values = c("#326747", "#98CDAD",
                                "#67325F", "#CD98C5")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(-.5, .5),
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "A)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme()  +
  theme(aspect.ratio = 1.15)

ggsave(p, file = "figures/model_figures/G_diversity_no_fishing.jpg",
       unit = "cm", height = 6.5, width = 13, dpi = 1200)

# Plot for models with fishing
p <- ggplot(data = f_data, aes(x = term, y = estimate, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  geom_point(size = 1,
             position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6, reverse = TRUE)) +
  scale_color_manual(values = c("#326747", "#98CDAD",
                                "#67325F", "#CD98C5")) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(-.5, .5),
    breaks = seq(-.5, .5, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5"),
    expand = c(0, 0)
  ) +
  labs(
    title = "B)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  coeff_theme() +
  theme(aspect.ratio = 1.15)

ggsave(p, file = "figures/model_figures/H_diversity_with_fishing.jpg",
       unit = "cm", height = 6.5, width = 13, dpi = 1200)

rm(list = setdiff(ls(), "map_theme"))

## ------------------------------------------------------------------------ ----
## Variable maps (grid) ----

# Load data
data <- read.csv("data/final/final_data.csv", sep = ",")

# Transform back to sf
plot_data <- data %>%
  select(haul_id, latitude,longitude, total_biomass, 
         kde_fric, kde_feve, tric, teve,
         log_sum_fishing,
         sbt_mean,
         gadus_morhua,         
         melanogrammus_aeglefinus,     
         sebastes_mentella,   
         pollachius_virens,         
         reinhardtius_hippoglossoides) %>%
  mutate(
    log_commercial_biomass = log(gadus_morhua + # demersal targets in west           
                                   melanogrammus_aeglefinus  +
                                   pollachius_virens +
                                   sebastes_mentella +
                                   reinhardtius_hippoglossoides + .1),
    log_total_biomass = log(total_biomass + .1),
         log_gadus_morhua = log(gadus_morhua + .1)) %>%
  sf::st_as_sf(
    coords = c("longitude", "latitude"), 
    crs = "+proj=aeqd +lat_0=76 +lon_0=27 +datum=WGS84 +units=km +no_defs") 

# Data limits
xlim <- st_bbox(plot_data)[c("xmin", "xmax")] + c(-5, 35)
ylim <- st_bbox(plot_data)[c("ymin", "ymax")] + c(0, 0) 

# Continuous grid to bin the observations (65 km cells)
g <- plot_data %>%
  sf::st_make_grid(cellsize = 65,  
                   what = "polygons",
                   square = TRUE, 
                   offset = st_bbox(plot_data)[c("xmin", "ymin")] +
                     c(-15, -30)) %>%
  sf::st_as_sf() 

# Spatial join (inner join) to match observations to cells
grid_r <- st_join(g, plot_data, left = FALSE) 

# Check observation number per cell (ideally = n of years)
cell_counts <- grid_r %>%
  group_by(x) %>%
  summarize(n_obs = n(), .groups = "drop")

prop.table(table((cell_counts$n_obs > 5))) 

# Responses to plot
res <- c(
  "log_total_biomass",
  "log_commercial_biomass",
  "log_gadus_morhua",
  "tric",
  "teve",
  "kde_fric",
  "kde_feve"
)

titles <- c(
  "Total biomass",
  "Commercial biomass",
  "Cod biomass",
  "Tric",
  "Teve",
  "Fric",
  "Feve"
  )

unit <- c(
  "log(kg/km2)",
  "log(kg/km2)",
  "log(kg/km2)",
  "N° species",
  rep(" ", 3)
)

# Variables' rounding 
round <- c(rep(10, 4), rep(100, 3))

# Land shape
land <- sf::st_crop(ne_countries(scale = "medium", returnclass = "sf"),
                    xmin = -10, xmax = 65, 
                    ymin = 55, ymax = 90) |>
  sf::st_transform(crs = st_crs(plot_data)) 

# Compute mean values over the study area, plot and save
for (i in 1:length(res)) { 
  
  grid_mean <- grid_r %>%
    group_by(x) %>%
    summarize(mean = mean(get(res[i]), na.rm = TRUE)) %>%
    st_as_sf()
  
  max <- ceiling(quantile(grid_mean$mean, .95, na.rm = TRUE) *
                   round[i])/round[i]
  min <- floor(quantile(grid_mean$mean, .05, na.rm = TRUE) *
                 round[i])/round[i]
  
  p <- ggplot() +
    geom_sf(data = grid_mean, aes(fill = mean), color = NA) +
    geom_sf(data = land, fill = "gray40", color = "gray40") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    scale_fill_gradientn(colors = c(
      "#003f5c", 
      "lightblue", 
      "gray85",
      "#e48646", 
      "darkred"),
      oob = scales::squish,
      limits = c(min, max), 
      name = paste0(unit[i], " \n"),
      breaks = seq(min, max, length = 2)) + 
    labs(x = "\nLongitude", y = "Latitude\n",
         title = titles[i]) + 
    guides(fill = guide_colorbar(barwidth = 1.5,
                                 barheight = .2,
                                 ticks = FALSE)) +
    map_theme() 
  
  ggsave(p, file = paste0("figures/maps/", res[i], "_grid.jpg"),
         width = 4, height = 6, units = "cm")
}

# Load environmental rasters
temperature_stack <- terra::rast("data/original/temperature_data.nc")
chla_stack <- terra::rast("data/original/chlorophyll_data.nc")
fishing_stack <- terra::rast("data/intermediate/fish_stack.tiff")

# Subset temperature data to compute each variable on its own
sst_stack <- temperature_stack[[grep("thetao", names(temperature_stack))]] # sst
sbt_stack <- temperature_stack[[grep("bottomT", names(temperature_stack))]]# sbt
sic_stack <- temperature_stack[[grep("siconc", names(temperature_stack))]] # sic

# Load community data for haul locations in space and time
load("data/intermediate/community_filtered.RData")

# Loop to compute mean and sd values, project, and plot
stacks <- list(
  SST  = "sst_stack",
  SBT  = "sbt_stack",
  SIC  = "sic_stack",
  Chla = "chla_stack"
)

unit_list <- c(
  SST  = "°C",
  SBT  = "°C",
  SIC  = "Cover",
  Chla = "mg/m3"
)

for (v in names(stacks)) {
  
  r <- get(stacks[[v]])
  
  # Keep only the time period used for modelling
  dates <- as.Date(time(r))
  
  keep <- dates >= as.Date("2003-08-01") & dates <= as.Date("2022-07-01")
  
  r <- r[[keep]]
  
  # Mean and SD across all months
  r_mean <- terra::app(r, mean, na.rm = TRUE)
  r_sd   <- terra::app(r, sd, na.rm = TRUE)
  
  # Match CRS of grid
  r_mean <- terra::project(r_mean, st_crs(g)$wkt)
  r_sd   <- terra::project(r_sd, st_crs(g)$wkt)
  
  # Extract mean raster value within each grid cell
  g_stats <- g
  
  g_stats$mean <- terra::extract(
    r_mean,
    terra::vect(g),
    fun = mean,
    na.rm = TRUE
  )[, 2]
  
  g_stats$SD <- terra::extract(
    r_sd,
    terra::vect(g),
    fun = mean,
    na.rm = TRUE
  )[, 2]
  
  for (stat in c("mean", "SD")) {
    
    vals <- g_stats[[stat]]
    
    max <- ceiling(quantile(vals, .95, na.rm = TRUE) *
                     10)/10
    min <- floor(quantile(vals, .05, na.rm = TRUE) *
                   10)/10
    
    p <- ggplot() +
      geom_sf(data = g_stats, aes(fill = .data[[stat]]), color = NA) +
    geom_sf(data = land, fill = "gray40", color = "gray40") +
      scale_x_continuous(limits = xlim) +
      scale_y_continuous(limits = ylim) +
      scale_fill_gradientn(colors = c(
        "#003f5c", 
        "lightblue", 
        "gray85",
        "#e48646", 
        "darkred"),
        na.value = NA,
        oob = scales::squish,
        limits = c(min, max), 
        name = paste0(unit_list[v], " \n"),
        breaks = seq(min, max, length = 2)) + 
      labs(x = "\nLongitude", y = "Latitude\n",
           title = paste(v, stat)) + 
      guides(fill = guide_colorbar(barwidth = 1.5,
                                   barheight = .2,
                                   ticks = FALSE)) +
      map_theme() 
    
    ggsave(paste0( "figures/maps/", v, "_", stat, "_grid.jpg"),
      p,
      width = 4,
      height = 6,
      units = "cm"
    )
    }
  }

# Same for fishing but using sum
r <- fishing_stack

dates <- as.Date(time(r))

keep <- dates >= as.Date("2011-08-01") & # minor bias (data starts in Jan 2012)
  dates <= as.Date("2022-07-01")

r <- r[[keep]]

# Sum
r_sum <- terra::app(r, sum, na.rm = TRUE)

# Match CRS of grid
r_sum <- terra::project(r_sum, st_crs(g)$wkt)

# Extract sum of raster value within each grid cell
g_stats <- g

g_stats$sum <- terra::extract(
  r_sum,
  terra::vect(g),
  fun = sum,
  na.rm = TRUE
)[, 2]

g_stats$sum <- log(g_stats$sum + 0.0001)
  
max <- ceiling(quantile(g_stats$sum , .95, na.rm = TRUE) *
                 10)/10
min <- floor(quantile(g_stats$sum , .05, na.rm = TRUE) *
               10)/10
  
  p <- ggplot() +
    geom_sf(data = g_stats, aes(fill = sum), color = NA) +
    geom_sf(data = land, fill = "gray40", color = "gray40") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    scale_fill_gradientn(colors = c(
      "#003f5c", 
      "lightblue", 
      "gray85",
      "#e48646", 
      "darkred"),
      na.value = NA,
      oob = scales::squish,
      limits = c(min, max), 
      name = "log(hours)\n",
      breaks = seq(min, max, length = 2)) + 
    labs(x = "\nLongitude", y = "Latitude\n",
         title = "Fishing effort") + 
    guides(fill = guide_colorbar(barwidth = 1.5,
                                 barheight = .2,
                                 ticks = FALSE)) +
    map_theme() 
  
  ggsave("figures/maps/fishing_grid.jpg",
         p, width = 4, height = 6, units = "cm")
  
# Depth
depth_rast <- rast("data/original/depth_gebco_400.tiff")

d <- depth_rast %>%
  aggregate(fact = 20, fun = mean, na.rm = TRUE) %>% # EPSG:3996
  terra::project(st_crs(g)$wkt) 

d[d > 0 | d < -500] <- NA

# Extract mean of raster value within each grid cell
g_stats <- g

g_stats$depth <- terra::extract(
  d,
  terra::vect(g),
  fun = mean,
  na.rm = TRUE
)[, 2]

max <- ceiling(quantile(g_stats$depth, .95, na.rm = TRUE))
min <- floor(quantile(g_stats$depth, .05, na.rm = TRUE))

p <- ggplot() +
  geom_sf(data = g_stats, aes(fill = depth), color = NA) +
  geom_sf(data = land, fill = "gray40", color = "gray40") +
  scale_x_continuous(limits = xlim) +
  scale_y_continuous(limits = ylim) +
  scale_fill_gradientn(colors = c(
    "#003f5c", 
    "lightblue", 
    "gray85",
    "#e48646", 
    "darkred"),
    na.value = NA,
    oob = scales::squish,
    limits = c(min, max), 
    name = "m\n",
    breaks = seq(min, max, length = 2)) + 
  labs(x = "\nLongitude", y = "Latitude\n",
       title = "Depth") + 
  guides(fill = guide_colorbar(barwidth = 1.5,
                               barheight = .2,
                               ticks = FALSE)) +
  map_theme() 

ggsave("figures/maps/depth_grid.jpg",
       p, width = 4, height = 6, units = "cm")

# Plot the number of hauls included in each cell and the hauls
p <- ggplot() +
  geom_sf(data = land, fill = "gray80", color = "gray80") + 
  geom_sf(data = plot_data, aes(color = "Hauls"), size = .5, alpha = 0.3) +
  geom_sf(data = cell_counts, fill = NA, color = "black", size = 0.1) + 
  geom_sf_text(data = cell_counts, aes(label = n_obs), 
               size = 1.5, color = "black") +
  scale_color_manual(
    name = NULL,
    values = c("Hauls" = "#CAAA98")
  ) +
  scale_x_continuous(limits = xlim) +
  scale_y_continuous(limits = ylim) +
  labs(
    x = "\nLongitude", y = "Latitude\n",
    title = "Hauls per Cell"
  ) +
  map_theme() 

ggsave(p, file = "figures/maps/haul_number.jpg",
       width = 4, height = 6, units = "cm")

rm(list = setdiff(ls(), c("plot_data", "xlim", "ylim", "map_theme")))

## Variable maps (point) ----

# Responses to plot
res <- c(
  "log_sum_fishing",
  "log_total_biomass",
  "log_commercial_biomass",
  "log_gadus_morhua",
  "sbt_mean",
  "tric",
  "teve",
  "kde_fric",
  "kde_feve"
)

titles <- c(
  "Fishing effort",
  "Total biomass",
  "Commercial biomass",
  "Cod biomass",
  "Mean SBT",
  "Tric",
  "Teve",
  "KDE Fric",
  "KDE Feve"
)

unit <- c(
  "log(hours/year)",
  "log(kg/km2)",
  "log(kg/km2)",
  "log(kg/km2)",
  "°C",
  "N° species",
  rep(" ", 3)
)

# Variables' rounding 
round <- c(rep(10, 6), rep(100, 3))

# Land shape
land <- sf::st_crop(ne_countries(scale = "medium", returnclass = "sf"),
                    xmin = -10, xmax = 65, 
                    ymin = 55, ymax = 90) |>
  sf::st_transform(crs = st_crs(plot_data)) 

# Compute mean values over the study area, plot and save
for (i in 1:length(res)) { 
  
  max <- ceiling(quantile(plot_data[[res[i]]], .95, na.rm = TRUE) * 
                   round[i])/round[i]
  min <- floor(quantile(plot_data[[res[i]]], .05, na.rm = TRUE) * 
                 round[i])/round[i]
  
  p <- ggplot() +
    geom_sf(data = plot_data,
            aes(color = get(res[i])),
            size = .5) +
    geom_sf(data = land, fill = "gray40", color = "gray40") +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    scale_color_gradientn(colors = c(
      "#003f5c", 
      "lightblue", 
      "gray85",
      "#e48646", 
      "darkred"),
      na.value = NA,
      oob = scales::squish,   # cap values outside limits
      limits = c(min, max), 
      name = paste0(unit[i], " \n"),
      breaks = seq(min, max, length = 2)) + 
    labs(x = "\nLongitude", y = "Latitude\n",
         title = titles[i]) + 
    guides(color = guide_colorbar(barwidth = 1.5,
                                 barheight = .2,
                                 ticks = FALSE)) +
    map_theme() 
  
  ggsave(p, file = paste0("figures/maps/point_version/", res[i], "_point.jpg"),
         width = 4, height = 6, units = "cm")
}

## End

