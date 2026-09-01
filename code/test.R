

## Partial biomass models withou FD or TD ----

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
    directory = "models/test",
    name = response,
    data = data,
    mesh_cutoff = 60, 
    family = lognormal(link = "log"),
    fixed_formula = as.formula(paste0( # this allows to loop the function
      new_res[i], "_biomass  ~", 
     # new_res[i], "_kde_fric +",
      #new_res[i], "_kde_feve +",
      new_res[i], "_tric +",
      new_res[i], "_teve + 
      sbt_mean +
      sbt_sd +
      sst_sd +        
      log_chla_mean +
      depth +
      log_sum_fishing"
    )))
  
  message(
    "\n======================================================================\n",
    paste0(new_res[i], " model fitted and inspected"),
    "\n======================================================================\n"
  )
  
}

csv_files <- list.files(
  path = "models/test",
  pattern = ".csv$",
  full.names = TRUE
)

combined_csv <- bind_rows(lapply(csv_files, read.csv)) |>
  filter(term !="(Intercept)")

write.csv(combined_csv,
          file = "models/test/biomass_models_coefficients.csv")


# Read biomass coefficient data
biomass_coeff <- read.csv(
  file = "models/test/biomass_models_coefficients.csv",
  sep = ",") 

remove <- paste0("remove_", 1:10, "_biomass")

# Identify levels
biomass_coeff$model <- factor(biomass_coeff$model, 
                              levels = paste0("remove_", 1:10, "_biomass"))

# Rename levels
levels(biomass_coeff$model) <- paste0("-Top ", 1:10)

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
  "Tric",
  "Teve",
  "log fish. eff."
)

# Plot
cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = biomass_coeff, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 1,
             position = position_dodge(width = 0.8, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.8, reverse = TRUE)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  scale_color_manual(values = cols[2:11]) +
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

ggsave(p, file = "models/test/D_partial_biomass_no_FD.png",
       unit = "cm", height = 9, width = 13)

rm(list = setdiff(ls(), keep))

## Plots of partial diversity ----
# Load data
data <- read.csv("data/final/final_data.csv", sep = ",")

# Transform back to sf
plot_data <- data %>%
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
  "total_biomass",
  paste0("remove_", 1:10, "_biomass")
)

titles <- res

unit <- rep("", 11)

# Variables' rounding 
round <- c(rep(100, 11))

# Land shape
land <- sf::st_crop(ne_countries(scale = "medium", returnclass = "sf"),
                    xmin = -10, xmax = 65, 
                    ymin = 55, ymax = 90) |>
  sf::st_transform(crs = st_crs(plot_data)) 

# Compute mean values over the study area, plot and save
for (i in 1:length(res)) { 
  
  grid_mean <- grid_r %>%
    group_by(x) %>%
    summarize(mean = mean(log(get(res[i]) + 0.0001), na.rm = TRUE)) %>%
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
  
  ggsave(p, file = paste0("figures/remove_biomass_maps/", res[i], "_grid.png"),
         width = 4, height = 6, units = "cm")
}

## check partial diversity effect ----
p <- read.csv(
  file = "models/biomass/different_diversity/biomass_models_coefficients.csv",
  sep = ",")  |> select(-X)

biomass_coeff <- read.csv(
  file = "models/biomass/biomass_models_coefficients.csv",
  sep = ",") |> select(-X)

plot_coeff <- rbind(p, general)

# Identify levels
plot_coeff$model <- factor(plot_coeff$model, 
                              levels = c(
                                "general_biomass",
                                paste0("remove_", c(1, 5, 10), "_biomass")
                              ))

# Rename levels
levels(plot_coeff$model) <- c(
  "General",
  paste0("-Top ", c(1, 5, 10))
)

plot_coeff <- plot_coeff %>%
  mutate(
    term = case_when(
      grepl("tric", term) ~ "tric",
      grepl("teve", term) ~ "teve",
      grepl("fric", term) ~ "fric",
      grepl("feve", term) ~ "feve",
      TRUE ~ term
    )
  )

# Identify levels
plot_coeff$term <- factor(plot_coeff$term, levels = c(
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
levels(plot_coeff$term) <- c(
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

cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = plot_coeff, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 1,
             position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  scale_color_manual(values = cols[c(1, 2, 6, 11)]) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(-.75, .75),
    breaks = seq(-.75, .75, by = .25),
    labels = c("-0.75", "-0.5", "-0.25", "0", "0.25", "0.5", "0.75"),
    expand = c(0, 0)
  ) +
  labs(
    title = "D) - partial diversity",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  theme(     
    aspect.ratio = 1.2,
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right",
    legend.justification = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(margin = margin(l = 2, r = 2), size = 10),
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

p

ggsave(p, file = "figures/model_figures/D_partial_biomass(partial_diversity).png",
       unit = "cm", height = 9, width = 13)

### Issues with boreal model ----
model <- manual_model_selection(
  directory = "models",
  name = "boreal_biomass",
  data = boreal_data,
  mesh_cutoff = 60,  
  random_selection = "spatiotemporal_random_fields_ar1",
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

boreal_data$latitude <- round(boreal_data$latitude)
boreal_data$longitude <- round(boreal_data$longitude)

### Issues with kde diversity models ----
model_gaus <- manual_model_selection(
  directory = "models",
  name = "kde_fric_gaus",
  data = general_data,
  mesh_cutoff = 60, 
  family = gaussian(link = "identity"),
  random_selection = "spatiotemporal_random_fields_ar1",
  fixed_formula = kde_fric ~
    sbt_mean +
    sbt_sd +
    #sst_mean +
    sst_sd +
    sic_mean + 
    #sic_sd +
    log_chla_mean +
    depth
)

sanity(model)$all_OK

## Partial biomass models excluding top species from diversity metrics ----

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

# read final data
final_data <- read.csv("data/final/final_data.csv")

# total community model
model_data <- final_data %>% 
  select(
    tric,
    teve, 
    kde_fric,
    kde_feve,
    sbt_mean,       
    sic_mean,
    sbt_sd,
    sst_sd,
    log_chla_mean,
    depth,
    haul_id,
    total_biomass,
    year, 
    latitude,
    longitude
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

mesh <- make_mesh(model_data, 
                  cutoff = 60,
                  xy_cols = c("longitude", "latitude"))

model <- sdmTMB(
  data = model_data,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  time = "year",
  reml = TRUE,
  formula = total_biomass ~  
    tric +
    teve + 
    kde_fric +
    kde_feve +
    sbt_mean +       
    sic_mean + 
    sbt_sd +
    sst_sd +
    log_chla_mean +
    depth
)

sanity(model)$all_ok

test_coeff <- data.frame(tidy(model), model = "General")

# Loops
new_res <- paste0("remove_", 1:10)

model_names <- paste0("-Top ", 1:10)

load("tools/vif_function.RData")

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
      "sic_mean",
      "depth",
      "log_chla_mean",
      "log_sum_fishing"
    ) %>%
    drop_na() %>%
    vif()
  
  print(x)
  
  }

for (i in 1:length(new_res)) {
  
new_data <- final_data %>%
  select(
    paste0(new_res[i], "_biomass"),
    paste0(new_res[i], "_tric"),
    paste0(new_res[i], "_teve"),
    paste0(new_res[i], "_kde_fric"),
    paste0(new_res[i], "_kde_feve"),
    "sbt_mean",
    "sbt_sd",
    "sst_sd",
    "sic_mean",
    "depth",
    "log_chla_mean",
    "log_sum_fishing",
    "year",
    "latitude",
    "longitude",
    "haul_id"
  ) %>%
  drop_na() %>%
  mutate(across(
    -c(
      haul_id,
      paste0(new_res[i], "_biomass"),
      year, 
      latitude,
      longitude
    ), ~ as.numeric(scale(.))))

mesh <- make_mesh(new_data, cutoff = 60, xy_cols = c("longitude", "latitude"))

model <- sdmTMB(
  data = new_data,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  time = "year",
  reml = TRUE,
  formula = as.formula(paste0( # this allows to loop the function
    new_res[i], "_biomass  ~",  
      new_res[i], "_tric +",
      new_res[i], "_teve +", 
      new_res[i], "_kde_fric +",
      new_res[i], "_kde_feve +
      sbt_mean +
      sbt_sd +
      sst_sd +
      sic_mean +        
      log_chla_mean +
      log_sum_fishing +
      depth"
  ))
)

print(new_res[i])
print(sanity(model)$all_ok)

coeff <- data.frame(tidy(model), model = model_names[i])

test_coeff <- rbind(test_coeff, coeff)

}

## Plot ----
plot_coeff <- drop_na(test_coeff) |> filter(term != "(Intercept)")

plot_coeff$model <- factor(plot_coeff$model, levels = c("General", model_names))

plot_coeff <- plot_coeff %>%
  mutate(
    term = case_when(
      grepl("tric", term) ~ "tric",
      grepl("teve", term) ~ "teve",
      grepl("fric", term) ~ "fric",
      grepl("feve", term) ~ "feve",
      TRUE ~ term
    )
  )

# Identify levels
plot_coeff$term <- factor(plot_coeff$term, levels = c(
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
levels(plot_coeff$term) <- c(
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

cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = plot_coeff, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 1,
             position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  scale_color_manual(values = cols) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(-.75, .75),
    breaks = seq(-.75, .75, by = .25),
    labels = c("-0.75", "-0.5", "-0.25", "0", "0.25", "0.5", "0.75"),
    expand = c(0, 0)
  ) +
  labs(
    title = "F)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  theme(     
    aspect.ratio = 1.2,
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right",
    legend.justification = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(margin = margin(l = 2, r = 2), size = 10),
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

p

ggsave(p, file = "figures/model_figures/partial_biomass_new_.png",
       unit = "cm", height = 9, width = 13)

## Old partial biomass exlusion ----

# Read final data
final_data <- read.csv("data/final/final_data.csv")

# Load filtered community data
load("data/intermediate/community_filtered.RData")

# Total community model
model_data <- final_data %>%
  select(
    total_biomass,
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
    depth,
   # log_sum_fishing,
    year,
    longitude,
    latitude,
    haul_id
  ) %>%
  drop_na() %>%
  mutate(
    across(
      -c(
        haul_id,
        total_biomass,
        year,
        latitude,
        longitude
      ),
      ~ as.numeric(scale(.))
    )
  )

mesh <- make_mesh(
  model_data,
  cutoff = 60,
  xy_cols = c("longitude", "latitude")
)

model <- sdmTMB(
  data = model_data,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial = "off",
  spatiotemporal = "ar1",
  time = "year",
  reml = TRUE,
  formula = total_biomass ~
   # log_sum_fishing +
    tric +
    teve +
    sbt_mean +
    sic_mean +
    sbt_sd +
    sst_sd +
    log_chla_mean +
    depth
)

sanity(model)$all_ok

test_coeff <- data.frame(
  tidy(model),
  model = "General"
)


# Sequential species removal
top_10 <- c(
  "Gadus morhua",
  "Melanogrammus aeglefinus",
  "Sebastes mentella",
  "Micromesistius poutassou",
  "Hippoglossoides platessoides",
  "Mallotus villosus",
  "Pollachius virens",
  "Boreogadus saida",
  "Trisopterus esmarkii",
  "Reinhardtius hippoglossoides"
)

new_res <- paste0("-Top ", 1:10)

other_vars <- final_data %>%
  select(
    sst_mean,
    sbt_mean,
    sic_mean,
    log_chla_mean,
    sbt_sd,
    sst_sd,
    sic_sd,
    log_chla_sd,
    #log_sum_fishing,
    depth,
    year,
    longitude,
    latitude,
    haul_id
  )

for (i in seq_along(new_res)) {
  
  new_data <- community_filtered %>%
    filter(!(taxon %in% top_10[1:i])) %>%
    group_by(haul_id) %>%
    summarize(
      total_biomass = sum(wgt_cpua, na.rm = TRUE),
      tric = length(num_cpua),
      teve = log(1 / sum((num_cpua / sum(num_cpua))^2)) /
        log(length(num_cpua)),
      .groups = "drop"
    ) %>%
    left_join(other_vars, by = "haul_id") %>%
    drop_na() %>%
    mutate(
      across(
        -c(
          haul_id,
          total_biomass,
          year,
          latitude,
          longitude
        ),
        ~ as.numeric(scale(.))
      )
    )
  
  mesh <- make_mesh(
    new_data,
    cutoff = 60,
    xy_cols = c("longitude", "latitude")
  )
  
  model <- sdmTMB(
    data = new_data,
    mesh = mesh,
    family = lognormal(link = "log"),
    spatial = "off",
    spatiotemporal = "ar1",
    time = "year",
    reml = TRUE,
    formula = total_biomass ~
     # log_sum_fishing +
      tric +
      teve +
      sbt_mean +
      sic_mean +
      sbt_sd +
      sst_sd +
      log_chla_mean +
      depth
  )
  
  print(new_res[i])
  print(sanity(model)$all_ok)
  
  coeff <- data.frame(
    tidy(model),
    model = new_res[i]
  )
  
  test_coeff <- rbind(test_coeff, coeff)
}

## Plot ----
plot_coeff <- drop_na(test_coeff) |> filter(term != "(Intercept)")

plot_coeff$model <- factor(plot_coeff$model, levels = c("General", model_names))

plot_coeff <- plot_coeff %>%
  mutate(
    term = case_when(
      grepl("tric", term) ~ "tric",
      grepl("teve", term) ~ "teve",
      TRUE ~ term
    )
  )

# Identify levels
plot_coeff$term <- factor(plot_coeff$term, levels = c(
  "depth", 
  "log_chla_mean",
  "sst_sd",
  "sbt_sd",
  "sic_mean",
  "sbt_mean",
  "Fric",
  "Feve",
  "tric",
  "teve",
  "log_sum_fishing"
))

# Rename levels
levels(plot_coeff$term) <- c(
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

cols <- colorRampPalette(c("black", "#2b9348", "#eeef20"))(11)

p <- ggplot(data = plot_coeff, aes(x = term, y = estimate, color = model)) +
  geom_point(size = 1,
             position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.6, reverse = TRUE)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) +
  scale_color_manual(values = cols) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    limits = c(-.5, 1),
    breaks = seq(-.5, .75, by = .25),
    labels = c("-0.5", "-0.25", "0", "0.25", "0.5", "0.75"),
    expand = c(0, 0)
  ) +
  labs(
    title = "F)",
    x = "Model term\n",
    y = "\nSlope ± CI",
    color = "Model") +
  theme(     
    aspect.ratio = 1.2,
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "right",
    legend.justification = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(0.2, "cm"),
    legend.text = element_text(margin = margin(l = 2, r = 2), size = 10),
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

p

ggsave(p, file = "figures/model_figures/partial_biomass_new_cpua.png",
       unit = "cm", height = 9, width = 13)

