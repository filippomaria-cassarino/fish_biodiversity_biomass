## File: biomass.R
## Purpose: compute taxonomic and functional diversity metrics
## Author: Filippomaria Cassarino
## Date: 29 Mar 2026

## Notes ----

# Run again to check that the exclusion of fdis has not caused any issues

## Library ---- 

# Install/load required packages
load("tools/install_load_packages_function.RData")

install_load_packages(
  
  # list of required packages
  c("dplyr",
    "tidyr",
    "ggplot2",
    "vegan",
    "tibble",
    "GGally",
    "mFD",
    "BAT",
    "modi",
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
  group_by(haul_id, haul_dur) %>%
  summarize(
    tric = length(num_cpua), 
    
    # Simpson's evenness: log(inverse simpson's diversity) / log(richness)
    teve = log(1 / sum((num_cpua / sum(num_cpua))^2)) / 
      log(length(num_cpua)),
    
    # Pielou's evenness: log(exp(Shannon's diversity)) / log(richness)
    p_teve = -sum(num_cpua / sum(num_cpua) * log(num_cpua / sum(num_cpua))) /
      log(length(num_cpua))
  )

save(taxonomic_diversity, file = "data/intermediate/taxonomic_diversity.RData")

# How many species per haul on average?
mean(taxonomic_diversity$tric) # 11.9

## Prepare data for functional diversity estimation ----

# This section imports and inspects the trait data, associates it to the 
# community data, builds a trait matrix and a species matrix, models
# a multidimensional functional (trait) space, and extracts functional
# diversity metrics from it. 

# Read trait data from Beukhof et al., 2019 -
# https://doi.org/10.1594/PANGAEA.900866
original_traits <- read.csv(
  file = "data/original/TraitCollectionFishNAtlanticNEPacificContShelf.csv",
  sep = ",")

# Community data
load("data/intermediate/community_filtered.RData")

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

# Traits to exclude
exclude_t <- c("habitat",         # its levels are messy: demersal is unclear
               "spawning_type",   # does not relate easily to climate or fishing 
               "ar")              # 28 species miss it
                    

# Species matrix: row = haul, column = taxon
sm <- community_filtered %>% 
  filter(!taxon %in% exclude_s) %>%                # exclude unwanted species
  dplyr::select(haul_id, taxon, count = num_cpua) %>% # biomass, not number
  spread(taxon, count) %>%                         # spread to wide format
  filter(rowSums(!is.na(.)) > 2) %>%  # exclude hauls with less than 3 species 
  mutate(across(everything(), ~ replace_na(., 0))) %>% # NAs are actually 0s
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
sum(pcoape$values$Relative_eig[1:n_dimensions]) * 100 # 96.9

# Extract PCoA scores to obtain a synthetic trait matrix
pcoa <- as.data.frame(pcoape$vectors[, 1:n_dimensions])

# Check hyperspace quality (0 to 1, here 0.99)
BAT::hyper.quality(gd, pcoa) 

rm(list = setdiff(ls(), c("pcoa", "sm", "n_dimensions")))

## Compute functional diversity -----

# This section computes functional diversity metrics 
# using probabilistic hypervolumes and convex hulls

# Computer cores
my.cores <- parallel::detectCores() - 1 

# Probabilistic hypervolumes
hv <- BAT::kernel.build(  
  comm = sm,                # species matrix
  trait = pcoa,             # synthetic trait matrix
  method.hv = "gaussian",   # recommended
  abund = TRUE,             # abundance data
  cores = my.cores
) 

kde_fric <- kernel.alpha(hv)    # functional richness 

kde_feve <- kernel.evenness(hv) # functional evenness 

kernel_diversity <- data.frame(
  haul_id = names(kde_fric),
  kde_fric,
  kde_feve
 
)

# Same diversity dimensions with a different methods (convex hull)
t <- as.matrix(pcoa)
s <- sm %>% filter(
  rowSums(. > 0, na.rm = TRUE) > n_dimensions ) %>% # convex hull requirement
  as.matrix()

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

# Join with main functional diversity
functional_diversity <- left_join(
  kernel_diversity,
  convex_diversity,
  by = "haul_id"
  )

## Save ----
save(functional_diversity,
     file = "data/intermediate/functional_diversity.RData")

## End