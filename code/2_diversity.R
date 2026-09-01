## File: biomass.R
## Purpose: compute taxonomic and functional diversity metrics
## Author: Filippomaria Cassarino
## Date: 29 Mar 2026

## Notes ----

## Library ---- 

# Install/load required packages
load("tools/install_load_packages_function.RData")

install_load_packages(
  
  # list of required packages
  c("dplyr",
    "tidyr",
    "ggplot2",
    "tibble",
    "purrr",
    "GGally",
    "gawdis",
    "mFD",
    "BAT",
    "janitor"
  ))

## Taxonomic diversity metrics ----

# This section computes the main taxonomic diversity indices at the 
# observation level (haul): richness, evenness (Pielou's), 
# and diversity (Shannon-Wiener)

# Load filtered community data and trait data
load("data/intermediate/community_filtered.RData")

# Compute richness, evenness (Simpson's), and diversity (Shannon-Wiener)
taxonomic_diversity <- community_filtered %>%
  group_by(haul_id) %>%
  summarize(
    tric = length(wgt_cpua), 
    
    # Simpson's evenness: log(inverse Simpson's diversity) / log(richness)
    teve = log(1 / sum((wgt_cpua / sum(wgt_cpua))^2)) / 
      log(length(wgt_cpua)), 
    
    # Pielou's evenness: log(exp(Shannon's diversity)) / log(richness)
    p_teve = -sum(wgt_cpua / sum(wgt_cpua) * log(wgt_cpua / sum(wgt_cpua))) /
      log(length(wgt_cpua))
  )

# Recalculate after progressively excluding top 10 species
load("data/intermediate/top_10_species.RData") # top_10 (species by biomass)

new_res <- paste0("remove_", 1:length(top_10))

for (i in 1:length(new_res)) {
  
  new_td <- community_filtered %>%
    dplyr::filter(!taxon %in% top_10[1:i]) %>%
    group_by(haul_id) %>%
    summarize(
      !!paste0(new_res[i], "_tric") := length(wgt_cpua), 
      !!paste0(new_res[i], "_teve") := log(1 / sum((wgt_cpua / sum(wgt_cpua))^2)) / 
        log(length(wgt_cpua))
    )
  
  taxonomic_diversity <- left_join(taxonomic_diversity, new_td, by = "haul_id")
  
  rm(new_td)
}

save(taxonomic_diversity, file = "data/intermediate/taxonomic_diversity.RData")

# How many species per haul on average?
mean(taxonomic_diversity$tric) # 11.9

mean(taxonomic_diversity$remove_10_tric, na.rm = TRUE) # 6.7

rm(list = ls())

## Prepare data for functional diversity estimation ----

# This section imports and inspects the trait data, associates it to the 
# community data, builds a trait matrix and a species matrix, models
# a multidimensional functional (trait) space, and extracts functional
# diversity metrics from it. 

# Output matrixes
species_matrixes <- list()

trait_matrixes <- list()

# Read trait data from Beukhof et al., 2019 -
# https://doi.org/10.1594/PANGAEA.900866
original_traits <- read.csv(
  file = "data/original/TraitCollectionFishNAtlanticNEPacificContShelf.csv",
  sep = ",")

# Community data and top 10 species
load("data/intermediate/community_filtered.RData")
load("data/intermediate/top_10_species.RData")

# Filter traits to keep the most accurate ones
traits <- original_traits %>%
  filter(
    LME == 20 # Barents Sea LME
  ) %>%
  
  # "Icelus" has no trait data for LME 20, but it does in 2 others
  rbind(filter(original_traits, taxon == "Icelus")[1, ]) %>% 
  mutate(
    habitat = case_when(
      habitat == "non-pelagic" ~ NA,  # convert "non-pelagic" to NA
      TRUE ~ habitat)) %>%
  janitor::clean_names()

# Merge traits with survey species and remove irrelevant columns
community_traits <- data.frame(taxon = unique(community_filtered$taxon)) %>%
  left_join(traits, by = c("taxon" = "taxon")) %>%
  dplyr::select(
    -which(grepl("reference", names(.)) | grepl("level", names(.))),      
    -c(
      family,
      genus,
      species,
      taxonomic_rank,
      fao,
      lme
    )) 

# Check how many traits are NAs
missing_traits <- community_traits[rowSums(is.na(community_traits)) > 0, ] 
colSums(is.na(missing_traits)) # Ar is missing for 28 species - remove
rowSums(is.na(missing_traits)) # 2 species miss all traits 
# (Liparis liparis and Zeugopterus norvegicus) 

# Species to exclude
exclude_s <- c("Liparis liparis",         # no trait data
               "Zeugopterus norvegicus",  # no trait data 
               "Somniosus microcephalus") # extreme traits, causes correlation

# How much biomass do they account for?
exclude_count <- community_filtered %>%
  filter(taxon %in% exclude_s) 

sum(exclude_count$wgt_cpua, na.rm = TRUE)/
  sum(community_filtered$wgt_cpua, na.rm = TRUE) * 100 # .3 % of total

# Traits to exclude
exclude_t <- c("habitat",         # its levels are messy: demersal is unclear
               "spawning_type",   # does not relate easily to climate or fishing 
               "ar")              # 28 species miss it


# Species matrix: row = haul, column = taxon
sm <- community_filtered %>% 
  filter(!taxon %in% exclude_s) %>%                # exclude unwanted species
  dplyr::select(haul_id, taxon, count = wgt_cpua) %>% # biomass, not number
  spread(taxon, count) %>%                         # spread to wide format
  mutate(across(everything(), ~ replace_na(., 0))) %>% # NAs are actually 0s
  filter(rowSums(across(
    -haul_id, ~ . > 0)) >= 3) %>% # exclude hauls with less than 3 species 
  column_to_rownames(var = "haul_id")              # set haul_id as row name 

# Traits matrix: row = taxon, column = trait 
tm <- community_traits %>%
  filter(!taxon %in% exclude_s) %>%     # exclude unwanted species
  dplyr::select(!all_of(exclude_t)) %>% # exclude unwanted traits
  column_to_rownames(var = "taxon") %>% # set taxon as row name
  mutate(across(where(is.character), as.factor)) %>%    # for gawdis::gawdis
  BAT::standard(trait = ., method = "z",                # scale (u = 0, sd = 1)
                convert = which(sapply(., is.numeric))) # continuous variables

# Check that tm's row names match sm's column names (must be TRUE)
all(sort(rownames(tm)) == sort(colnames(sm))) 

# Check collinearity among continuous traits
GGally::ggpairs(tm[, which(sapply(tm, is.numeric))]) 

tm <- select(tm, -c(
  length_infinity, # collinear with length_max
  age_max          # collinear with age_maturity
))

# Gower's distance with balanced contribution from each trait 
set.seed(8)
gd <- gawdis::gawdis(tm,
                     w.type = "optimized",
                     opti.maxiter = 1000)  

# body_shape and fin_shape have unbalanced distributions,
# but they are maintained
table(tm$body_shape)
table(tm$fin_shape)

# Check trait contribution to the dissimilarity - to balance trait effects, this
# should be similar
attr(gd,"correls")

# PCoA of the distance matrix 
pcoape <- ape::pcoa(gd)

# Check how many pcoa axes should be retained (best number has lowest deviation) 
y <- mFD::quality.fspaces(gd, deviation_weighting = "squared")
plot(y$quality_fspaces[, 1], type = "b")  
which.min(y$quality_fspaces[,1]) # 5 dimensions

n_dimensions <- 4 # similar score to 5, but greatly reduce computation 

# Variance explained by chosen axes 
sum(pcoape$values$Relative_eig[1:n_dimensions]) * 100 # 96.7

# Extract PCoA scores to obtain a synthetic trait matrix
pcoa <- as.data.frame(pcoape$vectors[, 1:n_dimensions])

# Check hyperspace quality (0 to 1, here 0.995)
BAT::hyper.quality(gd, pcoa) 

# Add matrixes in the lists
species_matrixes[["total_biomass"]] <- sm
trait_matrixes[["total_biomass"]] <- pcoa

# Repeat the same progressively excluding the 10 most abundant species
new_res <- paste0("remove_", 1:10)

checks <- data.frame(partial_community = new_res, 
                     best_n_dimensions = NA,
                     quality = NA)

for (i in 1:10) {
  
  # Species to exclude
  exclude_s <- c("Liparis liparis",         
                 "Zeugopterus norvegicus",  
                 "Somniosus microcephalus",
                 top_10[1:i])
  
  # Traits to exclude
  exclude_t <- c("habitat",         
                 "spawning_type",  
                 "ar",
                 "length_infinity", # drawing from total_biomass
                 "age_max")              
  
  
  # Species matrix: row = haul, column = taxon
  
  # To reduce computation, we will only extract metrics for hauls 
  # where there is also fishing data, which starts from 2012
  
  sm <- community_filtered %>% 
    filter(!taxon %in% exclude_s, year >= 2012) %>% # exclude hauls             
    dplyr::select(haul_id, taxon, count = wgt_cpua) %>% 
    spread(taxon, count) %>%                         
    mutate(across(everything(), ~ replace_na(., 0))) %>% 
    filter(rowSums(across(
      -haul_id, ~ . > 0)) >= 3) %>% 
    column_to_rownames(var = "haul_id")              
  
  # Traits matrix: row = taxon, column = trait 
  tm <- community_traits %>%
    filter(taxon %in% colnames(sm)) %>% # exclude hauls      
    dplyr::select(!all_of(exclude_t)) %>% 
    column_to_rownames(var = "taxon") %>% 
    mutate(across(where(is.character), as.factor)) %>%   
    BAT::standard(trait = ., method = "z",                
                  convert = which(sapply(., is.numeric))) 
  
  # Gower's distance with balanced contribution from each trait 
  set.seed(8)
  gd <- gawdis::gawdis(tm,
                       w.type = "optimized",
                       opti.maxiter = 1000)  
  
  # PCoA of the distance matrix 
  pcoape <- ape::pcoa(gd)
  
  # Check how many pcoa axes should be retained (best number has lowest deviation) 
  y <- mFD::quality.fspaces(gd, deviation_weighting = "squared")
  
  checks$best_n_dimensions[i] <- which.min(y$quality_fspaces[,1]) # to check
  
  n_dimensions <- 4 # keep 4 as before; check if the number is very different
  
  # Extract PCoA scores to obtain a synthetic trait matrix
  pcoa <- as.data.frame(pcoape$vectors[, 1:n_dimensions])
  
  # Check hyperspace quality 
  checks$quality[i] <- BAT::hyper.quality(gd, pcoa) 
  
  # Return in the lists
  species_matrixes[[new_res[i]]] <- sm
  
  trait_matrixes[[new_res[i]]] <- pcoa
  
  message(
    "\n=====================================================================\n",
    paste0("Matrixes for ", new_res[i], " assembled"),
    "\n=====================================================================\n")
}

# Checks are good: the best n_dimensions ranges from 4 to 6, while the quality
# is always > 0.99

# Save lists
save(trait_matrixes, species_matrixes, n_dimensions,
     file = "data/intermediate/functional_diversity_matrixes.RData")

rm(list = ls())

## Compute binary-hypervolume-based (convex hull) functional diversity -----

load("data/intermediate/functional_diversity_matrixes.RData")

# Convert to matrix objects and remove hauls with n species <= dimensions
t <- as.matrix(trait_matrixes[["total_biomass"]])

s <- species_matrixes[["total_biomass"]] %>% 
  filter(
    rowSums(. > 0, na.rm = TRUE) > n_dimensions ) %>% # convex hull requirement(4)
  as.matrix()

# Compute
fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord = t,    # trait matrix
  asb_sp_w = s,          # species matrix
  ind_vect = c(  
    "fric", # functional richness
    "feve" # functional evenness
  ))

# Extract the indexes from the list
convex_diversity <- fd$functional_diversity_indices %>%
  as.data.frame() %>%         
  tibble::rownames_to_column(var = "haul_id") %>%
  dplyr::select(
    haul_id,
    "ch_fric" = "fric",
    "ch_feve" = "feve"
  )

# Save
save(convex_diversity, 
     file = "data/intermediate/convex_diversity.RData")

rm(list = ls())

## Compute probabilistic-hypervolume-based (KDE) functional diversity ----

load("data/intermediate/functional_diversity_matrixes.RData")

# Computer cores
my.cores <- parallel::detectCores() - 1 

# Build probabilistic hypervolumes
set.seed(8)

hv <- BAT::kernel.build(  
  comm = species_matrixes[["total_biomass"]],       # species matrix
  trait = trait_matrixes[["total_biomass"]],        # synthetic trait matrix
  method.hv = "gaussian",   # recommended
  abund = TRUE,             # abundance data
  cores = my.cores
) 

# Extract metrics
kde_fric <- kernel.alpha(hv)    # functional richness 

kde_feve <- kernel.evenness(hv) # functional evenness 

kernel_diversity <- data.frame(
  haul_id = names(kde_fric),
  kde_fric,
  kde_feve
)

# Save
save(kernel_diversity, 
     file = "data/intermediate/kernel_diversity.RData")

rm(list = ls())

## Same for partial community matrixes (takes several hours to days)
load("data/intermediate/functional_diversity_matrixes.RData")

new_res <- paste0("remove_", 1:10)

my.cores <- parallel::detectCores() - 1 

# Output list
kd_list <- list()

for (i in 1:length(new_res)) {                                    
  
  # Keep the desired subset of hauls
  sm <- species_matrixes[[new_res[i]]]
  
  # Build probabilistic hypervolumes
  set.seed(8)
  
  hv <- BAT::kernel.build(  
    comm = sm,   # species matrix
    trait = trait_matrixes[[new_res[i]]],  # synthetic trait matrix
    method.hv = "gaussian",   # recommended
    abund = TRUE,             # abundance data
    cores = my.cores
  ) 
  
  # Extract metrics
  kde_fric <- kernel.alpha(hv)    # functional richness 
  
  kde_feve <- kernel.evenness(hv) # functional evenness 
  
  kd_list[[new_res[i]]] <- data.frame(
    haul_id = names(kde_fric),
    value1 = as.numeric(kde_fric),
    value2 = as.numeric(kde_feve)
  ) |>
    setNames(c(
      "haul_id",
      paste0(new_res[i], "_kde_fric"),
      paste0(new_res[i], "_kde_feve")
    ))
  
  rm(hv)
  
  message(
    "\n=====================================================================\n",
    paste0("Metrics have been successfully extracted for ", new_res[i]),
    "\n=====================================================================\n")
}

# Merge list elements
kernel_diversity_partial_biomass <- purrr::reduce(
  kd_list,
  left_join,
  by = "haul_id"
) 

# Save
save(kernel_diversity_partial_biomass, 
     file = "data/intermediate/kernel_diversity_partial_biomass.RData")

rm(list = ls())

## Merge functional diversity files and save ----

load("data/intermediate/convex_diversity.RData")

load("data/intermediate/kernel_diversity.RData")

load("data/intermediate/kernel_diversity_partial_biomass.RData")

# Join 
functional_diversity <- kernel_diversity %>%
  left_join(convex_diversity, by = "haul_id") %>%
  left_join(kernel_diversity_partial_biomass, by = "haul_id")

# Save
save(functional_diversity,
     file = "data/intermediate/functional_diversity.RData")

## End