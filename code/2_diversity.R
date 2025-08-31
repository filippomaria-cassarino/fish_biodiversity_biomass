## File: biomass.R
## Purpose: compute taxonomic and functional diversity metrics
## Author: Filippomaria Cassarino
## Date: 

## Notes ----

## Library ---- 

library(dplyr)
library(tidyr) # replace_na function
library(ggplot2)
library(vegan)
library(GGally) # collinearity plot
library(ape)
library(gawdis)
library(purrr)
library(tibble) # column and row names
library(hypervolume)
library(BAT)
library(fundiversity) # modern FD, compare with BAT
library(mFD) # how many pcoas to retain
library(terra)
library(sf)
library(rnaturalearth) # land shape, better and faster than mapdata
library(modi) # to compute cwv

library(forcats) # PCoA plot
library(ggrepel) # PCoA plot

library(janitor) # clean_names 

#
## Taxonomic diversity metrics ----

# Load filtered community data and trait data
load("data/intermediate/community_filtered.RData")

# Taxonomic richness, evenness (Pielou's), and diversity (Shannon-Wiener)
taxonomic_diversity <- community_filtered %>%
  group_by(haul_id) %>%
  summarize(Tric = vegan::specnumber(accepted_name),
            Teve = vegan::diversity(num_cpua)/log(specnumber(accepted_name)),
            Tdiv = vegan::diversity(num_cpua))

save(taxonomic_diversity, file = "data/intermediate/taxonomic_diversity.RData")

#
## Functional diversity metrics ----

# Load trait data
traits <- read.csv2("data/original/trait_beukhof.csv")

# Filter traits
trait <- filter(traits, LME == 20) %>%# retain only Barents, but remove "Icelus" 
  rbind(filter(traits, taxon == "Icelus")[1,]) %>% # add again "Icelus" 
  mutate(habitat = case_when(
    habitat == "non-pelagic" ~ NA,  # convert "non-pelagic" to NA
    TRUE ~ habitat)) 

# Merge traits with survey species and remove irrelevant columns
nor_trait <- data.frame(taxon = unique(community_filtered$accepted_name)) %>%
  left_join(trait, by = c("taxon" = "taxon")) %>% # joint to trait
  select(-which(grepl("reference", names(.)) |    # remove unwanted columns
                  grepl("level", names(.))),      
         -c(family,
            genus,
            species,
            taxonomic.rank,
            FAO,
            LME)) 

# Check how many NAs
missing_traits <- nor_trait[rowSums(is.na(nor_trait)) > 0, ]  
colSums(is.na(missing_traits)) # Ar is missing for 28 species - remove
rowSums(is.na(missing_traits)) # 2 species miss all traits (Liparis liparis and Zeugopterus norvegicus) - remove


# Exclude species
sum(community_filtered[community_filtered$accepted_name == "Somniosus microcephalus", ]$num) # extreme traits but only 11
sum(community_filtered[community_filtered$accepted_name == "Scomber scombrus", ]$num) # only pelagic species, just 3
 
exclude <- c("Liparis liparis",  "Zeugopterus norvegicus",  # no trait data
             "Somniosus microcephalus")                     # extreme traits             

# Categorical traits
table(as.factor(nor_trait$habitat)) # highly unbalanced (problems for gawdis) and inaccurate
table(as.factor(nor_trait$feeding.mode))
table(as.factor(nor_trait$body.shape)) # unbalanced for gawdis
table(as.factor(nor_trait$fin.shape))
table(as.factor(nor_trait$spawning.type)) # not meaningful for the study

# Species matrix: row = haul, column = taxon
sm <- community_filtered %>% 
  filter(!accepted_name %in% exclude) %>% # exclude unwanted species
  select(haul_id, taxon = accepted_name, count = num_cpue) %>%
  spread(taxon, count) %>% # spread to wide format
  filter(rowSums(!is.na(.)) > 2) %>%  # exclude hauls with less than 3 species (Mammola et al., 2024)
  mutate(across(everything(), ~ replace_na(., 0))) %>% # NAs are actually 0s
  column_to_rownames(var = "haul_id") # set haul_id as row name 

# Traits matrix: row = taxon, column = trait 
tm <- nor_trait %>%
  filter(!taxon %in% exclude) %>%   # exclude unwanted species
  select(!c(habitat,       # its levels are messy: demersal is unclear
            spawning.type, # only three levels, does not find any simple connection to climate or fishing 
            AR             # 28 species miss it
            )) %>% # exlude unwanted traits
  column_to_rownames(var = "taxon") %>% # set taxon as row name
  mutate(across(where(is.character), as.factor)) %>%    # set to factors to compute gawdis
  BAT::standard(trait = ., method = "standard",         # standardize 
                convert = which(sapply(., is.numeric))) # continuous variables

# Check that tm's row names match sm's column names (must be TRUE)
all(sort(rownames(tm)) == sort(colnames(sm))) 

# Check collinearity among continuous traits
GGally::ggpairs(tm[, which(sapply(tm, is.numeric))]) 

# length.max and length.infinity are highly collinear, but length.infinity has 2 more NAs
tm <- select(tm, ! length.infinity)

# Save matrixes
save(tm, file = "data/intermediate/trait_matrix.RData")
save(sm, file = "data/intermediate/species_matrix.RData")

## Gower's distance with balanced contribution from each trait (De Bello et al., 2021)
set.seed(8)
gd <- gawdis::gawdis(tm,
                     w.type = "optimized",
                     opti.maxiter = 1000)  

# habitat and body.shape have unbalanced distributions
table(tm$habitat)
table(tm$body.shape)

# Check their contribution to the dissimilarity (should be similar)
attr(gd,"correls")

# PCoA of the distance matrix and variance explained by first 3 axes (84.35948 %)
pcoape <- ape::pcoa(gd)

v <- pcoape$values$Relative_eig * 100
sum(v[1:3]) 

# Check how many pcoa axes should be retained (best number has lowest rmsd - 6) 
y <- mFD::quality.fspaces(gd, deviation_weighting = "squared")
y$quality_fspaces 
plot(y$quality_fspaces[, 1], type = "b")  

# Extract PCoA scores to obtain a synthetic trait matrix
pcoa <- as.data.frame(pcoape$vectors[, 1:3])

# Check hyperspace quality (0 to 1, here 0.9893658)
BAT::hyper.quality(gd, pcoa) 

# Save the synthetic trait matrix
save(pcoa, file = "data/intermediate/synthetic_trait_matrix.RData")

rm(list = ls())

# Compute functional diversity metrics
load("tools/Fdiversity_function.RData")
load("data/intermediate/synthetic_trait_matrix.RData")
load("data/intermediate/species_matrix.RData")

set.seed(8)
Fdiversity(species_matrix = sm, trait_matrix = pcoa)

# Check results
load("data/intermediate/functional_diversity.RData")

## End 

## IF NEEDED ----

## Check the effect of different trait number 

load("data/intermediate/trait_matrix.RData")
load("data/intermediate/species_matrix.RData")

# Select a subset of communities and adjust matrixes
set.seed(8)

sm_s <- sample_n(sm, 10)
sm_sample <- sm_s[, colSums(sm_s) > 0]

tm_s <- tm[colnames(sm_sample), ]
tm_sample <- BAT::fill(tm_s, method = "similar")

# All possible trait combinations (968)
trait_list <- names(tm)
all_combos <- lapply(2:10, function(k) {
  combn(trait_list, k, simplify = FALSE)
}) |> unlist(recursive = FALSE)

# List to store results
FD_list <- list()

# Loop to compute FD metrics for all trait combinations
for (i in 1:length(all_combos)) {
  cat("Running cycle", i, "of", length(all_combos), "\n")
  
  tryCatch({  # a few combinations will not work because two species share only one assessed trait
    set.seed(8)
    
    tm <- tm_sample[, all_combos[[i]]]
    
    gd <- gawdis::gawdis(tm, w.type = "optimized", opti.maxiter = 300)
    
    pcoa <- as.data.frame(ape::pcoa(gd)$vectors[, 1:3])
    
    set.seed(8)
    hv <- BAT::kernel.build(comm = sm_sample,
                            trait = pcoa,
                            method.hv = "gaussian",
                            abund = TRUE,
                            cores = 1)
    
    kde_ric <- kernel.alpha(hv)
    kde_dis <- kernel.dispersion(hv, func = "divergence", frac = 0.1)
    kde_eve <- kernel.evenness(hv)
    
    KDE_nor <- data.frame(haul_id   = names(kde_ric),
                          kde_ric   = kde_ric,
                          kde_eve   = kde_eve,
                          kde_dis   = kde_dis,
                          n_trait   = length(all_combos[[i]]),
                          trait_comb = paste(all_combos[[i]], collapse = "_"))
    
    FD_list[[i]] <- janitor::clean_names(KDE_nor)
    
  }, error = function(e) {
    cat("❌ Error in cycle", i, ":", conditionMessage(e), "\n")
    FD_list[[i]] <- NA                                 # mark failed runs
  })
  
  rm(list = setdiff(ls(), c("tm_sample", "sm_sample", "FD_list", "all_combos")))
}

trait_sensitivity <- do.call(rbind, FD_list)

save(trait_sensitivity, file = "data/final/trait_sensitivity.RData")


## dggrid
g <- dgconstruct(spacing = 100, resround='nearest')

data <- nor %>%
  group_by(haul_id, year, latitude, longitude) %>%
  summarize(depth = median(depth)) %>%
  left_join(Feve, by = c("haul_id" = "site")) %>%
  left_join(KDE_nor_all_sp, by = "haul_id") %>%
  left_join(KDE_nor, by = "haul_id") %>%
  left_join(TD_nor, by = "haul_id") %>%
  mutate(cell = dgGEO_to_SEQNUM(g, longitude, latitude)$seqnum) %>%
  group_by(cell, year) %>%               # group spatio-temporally
  summarize(KDEeve = mean(KDEeve.y),
            KDEeve_all = mean(KDEeve.x),
            KDEric = mean(KDEric.y),
            KDEric_all = mean(KDEric.x),
            KDEdis = mean(KDEdis.y),
            KDEdis_all = mean(KDEdis.x))

# Merge the grid boundaries with the binned data
plot_data <- dgcellstogrid(g, data$cell) %>% # cell boundaries 
  merge(data, 
        by.x = "seqnum",
        by.y = "cell") %>%
  st_transform(crs = "+proj=laea +lat_0=75 +lon_0=30 +datum=WGS84 +units=km +no_defs") 


# Range
max_val <- max(plot_data$KDEeve, na.rm = TRUE)
min_val <- min(plot_data$KDEeve, na.rm = TRUE)

# Plot
ggplot() +
  geom_sf(data = plot_data, aes(fill = KDEeve), color = NA) +
  geom_sf(data = land, fill = "gray45", color = "gray45") +
  geom_sf(data = test, color = "black", size = 1) +
  scale_x_continuous(limits = xlim) +
  scale_y_continuous(limits = ylim) +
  scale_fill_gradientn(colors = c("darkblue", "white", "darkred"),
                       values = scales::rescale(c(min_val, max_val)),
                       limits = c(min_val, max_val),
                       name = "KDEeve",
                       breaks = seq(min_val, max_val, length.out = 3)) +
  labs(x = "\nLongitude", y = "Latitude\n", title = "Map") +
  guides(fill = guide_colorbar(barwidth = .5, barheight = 10, ticks.linewidth = .1)) +
  theme_custom()
