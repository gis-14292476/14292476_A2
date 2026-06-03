# =========================================================
# 0. Project setup
# =========================================================
setwd("F:/14292476_A2")
set.seed(33)

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
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

# Spatial block validation and stable cluster contribution

# Create six spatial blocks

cell_xy <- st_coordinates(
  st_centroid(bird_fuzzy_data)
)

bird_fuzzy_data$cell_x <- cell_xy[, 1]
bird_fuzzy_data$cell_y <- cell_xy[, 2]

bird_fuzzy_data$x_block <- cut(
  bird_fuzzy_data$cell_x,
  breaks = quantile(
    bird_fuzzy_data$cell_x,
    probs = seq(0, 1, length.out = 4),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = c("X1", "X2", "X3")
)

bird_fuzzy_data$y_block <- cut(
  bird_fuzzy_data$cell_y,
  breaks = quantile(
    bird_fuzzy_data$cell_y,
    probs = seq(0, 1, length.out = 3),
    na.rm = TRUE
  ),
  include.lowest = TRUE,
  labels = c("Y1", "Y2")
)

bird_fuzzy_data$spatial_block <- interaction(
  bird_fuzzy_data$x_block,
  bird_fuzzy_data$y_block,
  drop = TRUE
)

bird_fuzzy_data$spatial_block <- as.factor(
  bird_fuzzy_data$spatial_block
)

table(bird_fuzzy_data$spatial_block)

plot(
  bird_fuzzy_data["spatial_block"],
  main = "Six spatial blocks for validation"
)

# Prepare model variables

bird_fuzzy_data$SR <- as.numeric(bird_fuzzy_data$SR)
bird_fuzzy_data$Abs <- as.numeric(bird_fuzzy_data$Abs)
bird_fuzzy_data$Abs_log <- log1p(bird_fuzzy_data$Abs)

mem_cols <- paste0("mem_cluster_", 1:best_k)

for (col in mem_cols) {
  bird_fuzzy_data[[col]] <- as.numeric(
    as.character(bird_fuzzy_data[[col]])
  )
}

all_blocks <- levels(bird_fuzzy_data$spatial_block)

train_block_combinations <- combn(
  all_blocks,
  4,
  simplify = FALSE
)

length(train_block_combinations)  # should be 15

# Helper functions

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

r2_pred <- function(obs, pred) {
  1 - sum((obs - pred)^2, na.rm = TRUE) /
    sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
}

standardise_by_train <- function(x, train_index) {
  (x - mean(x[train_index], na.rm = TRUE)) /
    sd(x[train_index], na.rm = TRUE)
}

# Exhaustive spatial block validation

coef_results <- data.frame()
validation_results <- data.frame()

for (i in seq_along(train_block_combinations)) {
  
  train_blocks_i <- train_block_combinations[[i]]
  test_blocks_i  <- setdiff(all_blocks, train_blocks_i)
  
  data_i <- bird_fuzzy_data
  
  data_i$split_i <- ifelse(
    data_i$spatial_block %in% train_blocks_i,
    "Train",
    "Test"
  )
  
  train_index <- data_i$split_i == "Train"
  test_index  <- data_i$split_i == "Test"
  
  # Standardise responses using training blocks only
  data_i$SR_z_i <- standardise_by_train(
    data_i$SR,
    train_index
  )
  
  data_i$Abs_log_z_i <- standardise_by_train(
    data_i$Abs_log,
    train_index
  )
  
  train_i <- data_i[train_index, ]
  test_i  <- data_i[test_index, ]
  
  # Build models
  sr_formula_i <- as.formula(
    paste(
      "SR_z_i ~ 0 +",
      paste(mem_cols, collapse = " + ")
    )
  )
  
  abs_formula_i <- as.formula(
    paste(
      "Abs_log_z_i ~ 0 +",
      paste(mem_cols, collapse = " + ")
    )
  )
  
  sr_model_i <- lm(
    sr_formula_i,
    data = train_i
  )
  
  abs_model_i <- lm(
    abs_formula_i,
    data = train_i
  )
  
  # Predict
  train_i$SR_pred_i  <- predict(sr_model_i, newdata = train_i)
  test_i$SR_pred_i   <- predict(sr_model_i, newdata = test_i)
  train_i$Abs_pred_i <- predict(abs_model_i, newdata = train_i)
  test_i$Abs_pred_i  <- predict(abs_model_i, newdata = test_i)
  
  # Store validation results
  validation_results <- rbind(
    validation_results,
    data.frame(
      split_id = i,
      response = c("SR", "SR", "Abs", "Abs"),
      dataset = c("Train", "Test", "Train", "Test"),
      RMSE = c(
        rmse(train_i$SR_z_i, train_i$SR_pred_i),
        rmse(test_i$SR_z_i, test_i$SR_pred_i),
        rmse(train_i$Abs_log_z_i, train_i$Abs_pred_i),
        rmse(test_i$Abs_log_z_i, test_i$Abs_pred_i)
      ),
      pred_R2 = c(
        r2_pred(train_i$SR_z_i, train_i$SR_pred_i),
        r2_pred(test_i$SR_z_i, test_i$SR_pred_i),
        r2_pred(train_i$Abs_log_z_i, train_i$Abs_pred_i),
        r2_pred(test_i$Abs_log_z_i, test_i$Abs_pred_i)
      ),
      train_blocks = paste(train_blocks_i, collapse = " + "),
      test_blocks = paste(test_blocks_i, collapse = " + ")
    )
  )
  
  # Store contribution coefficients
  coef_results <- rbind(
    coef_results,
    data.frame(
      split_id = i,
      response = "SR",
      cluster = names(coef(sr_model_i)),
      coefficient = as.numeric(coef(sr_model_i)),
      train_blocks = paste(train_blocks_i, collapse = " + "),
      test_blocks = paste(test_blocks_i, collapse = " + ")
    ),
    data.frame(
      split_id = i,
      response = "Abs",
      cluster = names(coef(abs_model_i)),
      coefficient = as.numeric(coef(abs_model_i)),
      train_blocks = paste(train_blocks_i, collapse = " + "),
      test_blocks = paste(test_blocks_i, collapse = " + ")
    )
  )
}

# Summarise coefficient stability

coef_stability <- coef_results %>%
  group_by(response, cluster) %>%
  summarise(
    mean_coef = mean(coefficient, na.rm = TRUE),
    sd_coef = sd(coefficient, na.rm = TRUE),
    min_coef = min(coefficient, na.rm = TRUE),
    max_coef = max(coefficient, na.rm = TRUE),
    positive_rate = mean(coefficient > 0, na.rm = TRUE),
    negative_rate = mean(coefficient < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    response,
    desc(mean_coef)
  )

coef_stability

# Summarise validation performance

validation_summary <- validation_results %>%
  group_by(response, dataset) %>%
  summarise(
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    sd_RMSE = sd(RMSE, na.rm = TRUE),
    mean_pred_R2 = mean(pred_R2, na.rm = TRUE),
    sd_pred_R2 = sd(pred_R2, na.rm = TRUE),
    min_pred_R2 = min(pred_R2, na.rm = TRUE),
    max_pred_R2 = max(pred_R2, na.rm = TRUE),
    .groups = "drop"
  )

validation_summary


# Plot coefficient stability

ggplot(
  coef_results,
  aes(
    x = cluster,
    y = coefficient
  )
) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(~ response) +
  theme_bw() +
  labs(
    title = "Stability of cluster contribution coefficients",
    x = "Fuzzy land-cover cluster",
    y = "Coefficient"
  )

# Plot validation performance

ggplot(
  validation_results,
  aes(
    x = dataset,
    y = pred_R2
  )
) +
  geom_boxplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(~ response) +
  theme_bw() +
  labs(
    title = "Spatial block validation performance",
    x = "Dataset",
    y = "Predictive R-squared"
  )

# Calculate habitat score using stable SR coefficients

stable_SR_coef <- coef_stability %>%
  filter(response == "SR") %>%
  mutate(
    cluster_number = as.numeric(
      gsub("mem_cluster_", "", cluster)
    )
  ) %>%
  arrange(cluster_number)

stable_SR_coef

sr_coef <- stable_SR_coef$mean_coef

membership_df <- st_drop_geometry(
  bird_fuzzy_data[, mem_cols]
)

membership_df[] <- lapply(
  membership_df,
  function(x) as.numeric(as.character(x))
)

membership_mat <- as.matrix(membership_df)

bird_fuzzy_data$habitat_SR_raw <- as.numeric(
  membership_mat %*% sr_coef
)

bird_fuzzy_data$habitat_score <- (
  bird_fuzzy_data$habitat_SR_raw -
    min(bird_fuzzy_data$habitat_SR_raw, na.rm = TRUE)
) / (
  max(bird_fuzzy_data$habitat_SR_raw, na.rm = TRUE) -
    min(bird_fuzzy_data$habitat_SR_raw, na.rm = TRUE)
)

# Plot stable SR contribution and habitat score
ggplot(
  stable_SR_coef,
  aes(
    x = reorder(cluster, mean_coef),
    y = mean_coef
  )
) +
  geom_col() +
  geom_errorbar(
    aes(
      ymin = mean_coef - sd_coef,
      ymax = mean_coef + sd_coef
    ),
    width = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Stable SR contribution of fuzzy land-cover clusters",
    x = "Fuzzy land-cover cluster",
    y = "Mean SR contribution coefficient"
  )

plot(
  bird_fuzzy_data["habitat_score"],
  main = "Habitat score based on stable SR contribution"
)

summary(
  bird_fuzzy_data[, c(
    "habitat_SR_raw",
    "habitat_score"
  )]
)

# =========================================================
# 11.Calculate neighbourhood effect
# =========================================================

hab_var <- "habitat_score"

# Find neighbouring cells
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

# Calculate relative neighbourhood effect
# Positive value:
# neighbouring cells are better than the focal cell

# Negative value:
# neighbouring cells are worse than the focal cell

bird_fuzzy_data$neighbour_effect <- 
  bird_fuzzy_data$neighbour_habitat_mean -
  bird_fuzzy_data[[hab_var]]


bird_fuzzy_data$neighbour_effect_type <- ifelse(
  bird_fuzzy_data$neighbour_effect > 0,
  "Positive",
  ifelse(
    bird_fuzzy_data$neighbour_effect < 0,
    "Negative",
    "Neutral"
  )
)

bird_fuzzy_data$neighbour_effect_type <- as.factor(
  bird_fuzzy_data$neighbour_effect_type
)

# Adjust habitat score using neighbourhood effect
alpha <- 0.3

bird_fuzzy_data$habitat_score_context <- 
  bird_fuzzy_data[[hab_var]] +
  alpha * bird_fuzzy_data$neighbour_effect


# Keep final score between 0 and 1
bird_fuzzy_data$habitat_score_context <- pmax(
  0,
  pmin(
    1,
    bird_fuzzy_data$habitat_score_context
  )
)

# Compare original and context-adjusted habitat scores

summary(
  bird_fuzzy_data[, c(
    "habitat_score",
    "neighbour_habitat_mean",
    "neighbour_effect",
    "habitat_score_context"
  )]
)

# Map neighbourhood effect

par(mfrow = c(1, 4))

plot(
  bird_fuzzy_data["habitat_score"],
  main = "Original habitat score"
)

plot(
  bird_fuzzy_data["neighbour_habitat_mean"],
  main = "Neighbour mean score"
)

plot(
  bird_fuzzy_data["neighbour_effect"],
  main = "Neighbour effect"
)

plot(
  bird_fuzzy_data["habitat_score_context"],
  main = "Context-adjusted score"
)

par(mfrow = c(1, 1))

# Map neighbourhood effect type
plot(
  bird_fuzzy_data["neighbour_effect_type"],
  main = "Neighbourhood effect type"
)

# =========================================================
# 12 Define potential habitat cells
# =========================================================
# Threshold sensitivity test

threshold_probs <- c(
  0.50,
  0.60,
  0.70,
  0.80,
  0.90
)

threshold_summary <- data.frame()

for (p in threshold_probs) {
  
  threshold_i <- quantile(
    bird_fuzzy_data$habitat_score_context,
    probs = p,
    na.rm = TRUE
  )
  
  habitat_i <- ifelse(
    bird_fuzzy_data$habitat_score_context >= threshold_i,
    1,
    0
  )
  
  threshold_summary <- rbind(
    threshold_summary,
    data.frame(
      threshold_prob = p,
      threshold_value = as.numeric(threshold_i),
      n_habitat_cells = sum(habitat_i == 1, na.rm = TRUE),
      prop_habitat_cells = mean(habitat_i == 1, na.rm = TRUE)
    )
  )
}

threshold_summary

# Choose main habitat threshold

# Use the top 30% of context-adjusted habitat scores as potential habitat.
# This corresponds to the 0.70 quantile threshold.

hab_threshold_prob <- 0.80

hab_threshold <- quantile(
  bird_fuzzy_data$habitat_score_context,
  probs = hab_threshold_prob,
  na.rm = TRUE
)

hab_threshold

# Define potential habitat cells

bird_fuzzy_data$potential_habitat <- ifelse(
  bird_fuzzy_data$habitat_score_context >= hab_threshold,
  1,
  0
)

bird_fuzzy_data$potential_habitat <- as.factor(
  bird_fuzzy_data$potential_habitat
)

# Plot potential habitat cells

plot(
  bird_fuzzy_data["potential_habitat"],
  main = "Potential habitat cells"
)

# =========================================================
# 13.Select potential habitat cells
# =========================================================

hab_cells <- bird_fuzzy_data %>%
  filter(potential_habitat == 1)

if (nrow(hab_cells) == 0) {
  
  stop("No habitat cells found. Try a lower threshold.")
  
}

# Identify connected habitat patches

if (nrow(hab_cells) == 1) {
  
  hab_cells$patch_id <- 1
  
} else {
  
  # Find touching habitat cells
  touch_list_hab <- st_touches(hab_cells)
  
  # Build graph from touching relationships
  hab_graph <- graph_from_adj_list(
    touch_list_hab,
    mode = "all"
  )
  
  # Find connected components
  hab_comp <- components(hab_graph)
  
  # Add patch ID
  hab_cells$patch_id <- hab_comp$membership
}

# Dissolve cells into habitat patches

hab_patches <- hab_cells %>%
  group_by(patch_id) %>%
  summarise(
    n_cells = n(),
    mean_habitat_score = mean(
      habitat_score_context,
      na.rm = TRUE
    ),
    max_habitat_score = max(
      habitat_score_context,
      na.rm = TRUE
    ),
    mean_neighbour_effect = mean(
      neighbour_effect,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Calculate patch-level metrics

hab_patches$patch_area_m2 <- as.numeric(
  st_area(hab_patches)
)

hab_patches$patch_area_km2 <- hab_patches$patch_area_m2 / 1e6

hab_patches$patch_perimeter_m <- as.numeric(
  st_length(
    st_boundary(hab_patches)
  )
)

hab_patches$patch_perimeter_km <- hab_patches$patch_perimeter_m / 1000

# Shape complexity:

hab_patches$shape_complexity <- 
  hab_patches$patch_perimeter_km /
  (2 * sqrt(pi * hab_patches$patch_area_km2))

# Summarise habitat patch structure

patch_summary <- data.frame(
  n_habitat_cells = nrow(hab_cells),
  n_patches = nrow(hab_patches),
  total_habitat_area_km2 = sum(
    hab_patches$patch_area_km2,
    na.rm = TRUE
  ),
  mean_patch_area_km2 = mean(
    hab_patches$patch_area_km2,
    na.rm = TRUE
  ),
  max_patch_area_km2 = max(
    hab_patches$patch_area_km2,
    na.rm = TRUE
  ),
  mean_patch_score = mean(
    hab_patches$mean_habitat_score,
    na.rm = TRUE
  )
)

patch_summary


# Plot habitat patches

plot(
  hab_patches["patch_id"],
  main = "Potential habitat patches"
)

plot(
  hab_patches["mean_habitat_score"],
  main = "Mean habitat score by patch"
)

# Check largest patches

largest_patches <- hab_patches %>%
  st_drop_geometry() %>%
  arrange(desc(patch_area_km2)) %>%
  select(
    patch_id,
    n_cells,
    patch_area_km2,
    mean_habitat_score,
    max_habitat_score,
    shape_complexity
  )

largest_patches


# =========================================================
# 15. Validate whether potential habitat captures Abs
# =========================================================

# Calculate Abs captured by potential habitat cells

abs_capture_summary <- bird_fuzzy_data %>%
  st_drop_geometry() %>%
  group_by(potential_habitat) %>%
  summarise(
    n_cells = n(),
    total_Abs = sum(Abs, na.rm = TRUE),
    mean_Abs = mean(Abs, na.rm = TRUE),
    median_Abs = median(Abs, na.rm = TRUE),
    .groups = "drop"
  )

abs_capture_summary$total_cells <- sum(abs_capture_summary$n_cells)

abs_capture_summary$total_Abs_all <- sum(
  abs_capture_summary$total_Abs,
  na.rm = TRUE
)

abs_capture_summary$cell_share <- 
  abs_capture_summary$n_cells / abs_capture_summary$total_cells

abs_capture_summary$Abs_share <- 
  abs_capture_summary$total_Abs / abs_capture_summary$total_Abs_all

abs_capture_summary$Abs_enrichment <- 
  abs_capture_summary$Abs_share / abs_capture_summary$cell_share

abs_capture_summary

# Plot cell share vs Abs share

abs_capture_plot <- abs_capture_summary %>%
  filter(potential_habitat == 1)

capture_compare <- data.frame(
  metric = c(
    "Cell share",
    "Abs share"
  ),
  value = c(
    abs_capture_plot$cell_share,
    abs_capture_plot$Abs_share
  )
)

ggplot(
  capture_compare,
  aes(
    x = metric,
    y = value
  )
) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Abs captured by potential habitat cells",
    x = "",
    y = "Proportion"
  )