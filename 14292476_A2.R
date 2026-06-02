# =========================================================
# 0. Project setup
# =========================================================
setwd("F:/14292476_A2")
set.seed(33)

library(sf)
library(terra)
library(dplyr)
library(tidyr)

# =========================================================
# 1. Read and preprocess spatial data
# =========================================================
# Read data
birds_grid <- st_read("./data/Birds/GM_Birds_2025.shp")
birds_lcm  <- rast("./data/Birds/gm_lcm_2022.tif")

# Data preprocess
# geometry correction
birds_grid <- st_make_valid(birds_grid)
# Create grid ID
birds_grid$grid_id <- seq_len(nrow(birds_grid))

# Check CRS
crs(birds_grid) == crs(birds_lcm)

# =========================================================
# 2. Clip and extract land-cover raster values
# =========================================================
# Clip land-cover
birds_vect <- vect(birds_grid)

birds_lcm_clip <- crop(birds_lcm, birds_vect)
birds_lcm_mask <- mask(birds_lcm_clip, birds_vect)

# Extract land-cover in each grid
lc_extract <- terra::extract( birds_lcm_mask, birds_vect )
lc_extract <- lc_extract[, 1:2]
names(lc_extract) <- c("grid_id", "lc_class")
# Delete 0 = NA
lc_extract <- subset(lc_extract,!is.na(lc_class) & lc_class != 0)

# =========================================================
# 3. Reclassify land-cover codes
# =========================================================
# Land-cover class table
lc_n2c <- data.frame(
  lc_class = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  lc_name = c("Woodland","Conifer","Arable",
              "Grassland","Moor","Bog",
              "Water","Saltmarsh","Urban")
)
lc_extract <- merge(lc_extract,lc_n2c,by = "lc_class",all.x = TRUE)

# =========================================================
# 4. Calculate land-cover composition per grid cell
# =========================================================
# Count the number for each land-cover
lc_counts <- aggregate(
  lc_class ~ grid_id + lc_name,
  data = lc_extract,
  FUN = length )

names(lc_counts)[3] <- "n_pixels"

# Count total pixels by fishnet
lc_totals <- aggregate(
  n_pixels ~ grid_id, 
  data = lc_counts, 
  FUN = sum )

names(lc_totals)[2] <- "total_pixels"

# Join pixel counts and total pixels
lc_prop_long <- merge( 
  lc_counts, lc_totals,
  by = "grid_id" )

# Calculate land-cover proportions within each grid cell
lc_prop_long$prop <- lc_prop_long$n_pixels / lc_prop_long$total_pixels

# Count land-cover types by grid cell
lc_n_types <- aggregate(
  lc_name ~ grid_id,
  data = lc_prop_long,
  FUN = function(x) length(unique(x))
)

names(lc_n_types)[2] <- "lc_n_types"

# Create land-cover combination by grid cell
lc_combo <- aggregate(
  lc_name ~ grid_id,
  data = lc_prop_long,
  FUN = function(x) paste(sort(unique(x)), collapse = " + ")
)

names(lc_combo)[2] <- "lc_combo"

# Reshape from long to wide format
lc_prop_wide <- tidyr::pivot_wider(
  lc_prop_long[, c("grid_id", "lc_name", "prop")],
  names_from = lc_name,
  values_from = prop,
  values_fill = 0
)

# Combine type count, combination, and proportions
lc_grid_summary <- merge( lc_n_types, lc_combo, by = "grid_id" )
lc_grid_summary <- merge( lc_grid_summary, lc_prop_wide, by = "grid_id" )

# =========================================================
# 5. Summarise land-cover combinations
# =========================================================
# Count different land-cover combinations
n_lc_combos <- length(unique(lc_grid_summary$lc_combo))

# Count the number of fishnet cells for each combination
lc_combo_count <- as.data.frame(table(lc_grid_summary$lc_combo))
names(lc_combo_count) <- c("lc_combo", "n_fishnet")

# Sort from most common to least common
lc_combo_count <- lc_combo_count[
  order(lc_combo_count$n_fishnet, decreasing = TRUE),
]

n_lc_combos
lc_combo_count

# =========================================================
# 6. Cluster grid cells by land-cover composition
# =========================================================
# Prepare land-cover proportion data
lc_prop_cols <- lc_n2c$lc_name
missing_cols <- setdiff(lc_prop_cols,names(lc_grid_summary))
for (col in missing_cols) {lc_grid_summary[[col]] <- 0}

# land-cover proportion columns for clustering
lc_kmeans_data <- lc_grid_summary[, lc_prop_cols]

# Replace NA with zero
lc_kmeans_data[is.na(lc_kmeans_data)] <- 0

# Use grid_id as row names
rownames(lc_kmeans_data) <- lc_grid_summary$grid_id


# Elbow method for choosing k

max_k <- 15
wss <- numeric(max_k)

for (k in 1:max_k) {
  
  km <- kmeans(
    lc_kmeans_data,
    centers = k,
    nstart = 100
  )
  
  wss[k] <- km$tot.withinss
  
}

plot(
  1:max_k,
  wss,
  type = "b",
  pch = 19,
  xlab = "Number of clusters (k)",
  ylab = "Total within-cluster sum of squares",
  main = "Elbow method for land-cover clustering"
)


# Final k-means clustering

best_k <- 5

lc_kmeans <- kmeans(
  lc_kmeans_data,
  centers = best_k,
  nstart = 100
)

# Add hard cluster labels
lc_grid_summary$lc_cluster <- lc_kmeans$cluster


# Summarise land-cover composition of each cluster
lc_cluster_summary <- aggregate(
  lc_grid_summary[, lc_prop_cols],
  by = list(lc_cluster = lc_grid_summary$lc_cluster),
  FUN = mean
)

lc_cluster_summary

# Map land-cover clusters

birds_grid_map <- merge(
  birds_grid,
  lc_grid_summary[, c("grid_id", "lc_cluster")],
  by = "grid_id",
  all.x = TRUE
)

birds_grid_map$lc_cluster <- as.factor(birds_grid_map$lc_cluster)

plot(
  birds_grid_map["lc_cluster"],
  main = "Spatial distribution of land-cover clusters"
)

# =========================================================
# 7. Calculate fuzzy membership to land-cover clusters
# =========================================================
# Calculate distance to each cluster center

cluster_centres <- lc_kmeans$centers

dist_to_centres <- matrix(
  NA,
  nrow = nrow(lc_kmeans_data),
  ncol = nrow(cluster_centres)
)

for (i in 1:nrow(lc_kmeans_data)) {
  
  for (j in 1:nrow(cluster_centres)) {
    
    dist_to_centres[i, j] <- sqrt(
      sum(
        (
          as.numeric(lc_kmeans_data[i, ]) -
            as.numeric(cluster_centres[j, ])
        )^2
      )
    )
    
  }
  
}

colnames(dist_to_centres) <- paste0("dist_cluster_", 1:best_k)
rownames(dist_to_centres) <- rownames(lc_kmeans_data)

head(dist_to_centres)


# Convert distance to fuzzy membership

# Convert distance to similarity.
similarity <- 1 / (dist_to_centres + 1e-6)

cluster_membership <- similarity / rowSums(similarity)

cluster_membership <- as.data.frame(cluster_membership)

names(cluster_membership) <- paste0("mem_cluster_", 1:best_k)

cluster_membership$grid_id <- lc_grid_summary$grid_id

head(cluster_membership)

# =========================================================
# 8. Join bird data with land-cover cluster membership
# =========================================================
# Join fuzzy membership back to land-cover grid summary
lc_grid_summary <- merge(
  lc_grid_summary,
  cluster_membership,
  by = "grid_id",
  all.x = TRUE
)

names(lc_grid_summary)

# Join bird response variables with cluster membership
bird_fuzzy_data <- merge(
  birds_grid,
  lc_grid_summary[, c("grid_id", paste0("mem_cluster_", 1:best_k), "lc_cluster")],
  by = "grid_id",
  all.x = TRUE
)

# =========================================================
# 9. Prepare bird response variables
# =========================================================
bird_fuzzy_data$SR <- as.numeric(bird_fuzzy_data$SR)
bird_fuzzy_data$Abs <- as.numeric(bird_fuzzy_data$Abs)

# Standardise SR and Abs
bird_fuzzy_data$SR_z <- as.numeric(scale(bird_fuzzy_data$SR))

bird_fuzzy_data$Abs_log <- log1p(bird_fuzzy_data$Abs)

bird_fuzzy_data$Abs_log_z <- as.numeric(scale(bird_fuzzy_data$Abs_log))

# =========================================================
# 10. Estimate cluster contribution from bird responses
# =========================================================
mem_cols <- paste0("mem_cluster_", 1:best_k)

sr_formula <- as.formula(
  paste(
    "SR_z ~ 0 +",
    paste(mem_cols, collapse = " + ")
  )
)

sr_cluster_model <- lm(
  sr_formula,
  data = bird_fuzzy_data
)

summary(sr_cluster_model)


# Model cluster contribution to Abs

abs_formula <- as.formula(
  paste(
    "Abs_log_z ~ 0 +",
    paste(mem_cols, collapse = " + ")
  )
)

abs_cluster_model <- lm(
  abs_formula,
  data = bird_fuzzy_data
)

summary(abs_cluster_model)

# Summarise cluster contribution

sr_contribution <- data.frame(
  lc_cluster = 1:best_k,
  SR_contribution = coef(sr_cluster_model)
)

abs_contribution <- data.frame(
  lc_cluster = 1:best_k,
  Abs_contribution = coef(abs_cluster_model)
)

cluster_contribution <- merge(
  sr_contribution,
  abs_contribution,
  by = "lc_cluster"
)

cluster_contribution

mem_cols <- paste0("mem_cluster_", 1:best_k)
mem_cols %in% names(bird_fuzzy_data)

# =========================================================
# 11. Calculate habitat suitability score
# =========================================================
# Extract membership columns and drop geometry
membership_df <- st_drop_geometry(
  bird_fuzzy_data[, mem_cols]
)
# Convert all membership columns to numeric
membership_df[] <- lapply(
  membership_df,
  function(x) as.numeric(as.character(x))
)

# Convert to numeric matrix
membership_mat <- as.matrix(membership_df)

# Make sure contribution values are ordered correctly
cluster_contribution <- cluster_contribution[
  order(cluster_contribution$lc_cluster),
]

sr_coef <- cluster_contribution$SR_contribution
abs_coef <- cluster_contribution$Abs_contribution

# =========================================================
# Calculate SR-based habitat score
bird_fuzzy_data$habitat_SR <- as.numeric(
  membership_mat %*% sr_coef
)
# Calculate Abs-based habitat score
bird_fuzzy_data$habitat_Abs <- as.numeric(
  membership_mat %*% abs_coef
)

# Standardise habitat scores to 0-1
bird_fuzzy_data$habitat_SR_01 <- (
  bird_fuzzy_data$habitat_SR -
    min(bird_fuzzy_data$habitat_SR, na.rm = TRUE)
) / (
  max(bird_fuzzy_data$habitat_SR, na.rm = TRUE) -
    min(bird_fuzzy_data$habitat_SR, na.rm = TRUE)
)

bird_fuzzy_data$habitat_Abs_01 <- (
  bird_fuzzy_data$habitat_Abs -
    min(bird_fuzzy_data$habitat_Abs, na.rm = TRUE)
) / (
  max(bird_fuzzy_data$habitat_Abs, na.rm = TRUE) -
    min(bird_fuzzy_data$habitat_Abs, na.rm = TRUE)
)


# Calculate combined habitat score
bird_fuzzy_data$habitat_score <- (
  bird_fuzzy_data$habitat_SR_01 *
    bird_fuzzy_data$habitat_Abs_01
) 


# Plot habitat scores

par(mfrow = c(1, 3))

plot(
  bird_fuzzy_data["habitat_SR_01"],
  main = "SR-based habitat score"
)

plot(
  bird_fuzzy_data["habitat_Abs_01"],
  main = "Abs-based habitat score"
)

plot(
  bird_fuzzy_data["habitat_score"],
  main = "Combined habitat score"
)

par(mfrow = c(1, 1))

# Check habitat score summary
summary(
  bird_fuzzy_data[, c(
    "habitat_SR",
    "habitat_Abs",
    "habitat_SR_01",
    "habitat_Abs_01",
    "habitat_score"
  )]
)

# =========================================================
# 12. Calculate neighbourhood effect
# =========================================================
hab_var <- "habitat_score"

touch_list <- st_touches(bird_fuzzy_data)

# Count number of neighbours for each cell
bird_fuzzy_data$n_neighbours <- lengths(touch_list)

# Calculate mean habitat score of neighbouring cells

bird_fuzzy_data$neighbour_habitat_mean <- NA

for (i in seq_len(nrow(bird_fuzzy_data))) {
  
  neigh_ids <- touch_list[[i]]
  
  if (length(neigh_ids) > 0) {
    
    bird_fuzzy_data$neighbour_habitat_mean[i] <- mean(
      bird_fuzzy_data[[hab_var]][neigh_ids],
      na.rm = TRUE
    )
    
  }
}

# Calculate whether neighbourhood effect is positive or negative

global_habitat_mean <- mean(
  bird_fuzzy_data[[hab_var]],
  na.rm = TRUE
)

bird_fuzzy_data$neighbour_effect <- 
  bird_fuzzy_data$neighbour_habitat_mean - global_habitat_mean

bird_fuzzy_data$neighbour_effect_type <- ifelse(
  bird_fuzzy_data$neighbour_effect >= 0,
  "Positive",
  "Negative"
)

bird_fuzzy_data$neighbour_effect_type <- as.factor(
  bird_fuzzy_data$neighbour_effect_type
)

# Strengthen or weaken habitat score

alpha <- 0.3

bird_fuzzy_data$habitat_score_context <- 
  bird_fuzzy_data[[hab_var]] +
  alpha * bird_fuzzy_data$neighbour_effect

# Keep final score between 0 and 1
bird_fuzzy_data$habitat_score_context <- pmax(
  0,
  pmin(1, bird_fuzzy_data$habitat_score_context)
)


# Map neighbourhood effect

par(mfrow = c(1, 3))

plot(
  bird_fuzzy_data["neighbour_habitat_mean"],
  main = "Mean habitat score of 8 neighbours"
)

plot(
  bird_fuzzy_data["neighbour_effect"],
  main = "Neighbourhood effect"
)

plot(
  bird_fuzzy_data["habitat_score_context"],
  main = "Context-adjusted habitat score"
)

par(mfrow = c(1, 1))

# =========================================================
# 13. Define potential habitat cells
# =========================================================
hab_threshold <- quantile(
  bird_fuzzy_data$habitat_score_context,
  probs = 0.70,
  na.rm = TRUE
)

hab_threshold

bird_fuzzy_data$potential_habitat <- ifelse(
  bird_fuzzy_data$habitat_score_context >= hab_threshold,
  1,0
)

bird_fuzzy_data$potential_habitat <- as.factor(
  bird_fuzzy_data$potential_habitat
)

plot(
  bird_fuzzy_data["potential_habitat"],
  main = "Potential habitat cells"
)
