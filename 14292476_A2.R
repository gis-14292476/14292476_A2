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

birds_grid <- st_make_valid(birds_grid)
birds_grid$grid_id <- seq_len(nrow(birds_grid))

# Check CRS
crs(birds_grid) == crs(birds_lcm)

# =========================================================
# 2. Clip and extract land-cover raster values
# =========================================================
birds_vect <- vect(birds_grid)

birds_lcm_clip <- crop(birds_lcm, birds_vect)
birds_lcm_mask <- mask(birds_lcm_clip, birds_vect)

lc_extract <- terra::extract( birds_lcm_mask, birds_vect )
lc_extract <- lc_extract[, 1:2]
names(lc_extract) <- c("grid_id", "lc_class")
lc_extract <- subset(lc_extract,!is.na(lc_class) & lc_class != 0)

# =========================================================
# 3. Reclassify land-cover codes
# =========================================================
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
lc_counts <- aggregate(
  lc_class ~ grid_id + lc_name,
  data = lc_extract,
  FUN = length )

names(lc_counts)[3] <- "n_pixels"

lc_totals <- aggregate(
  n_pixels ~ grid_id,
  data = lc_counts,
  FUN = sum )

names(lc_totals)[2] <- "total_pixels"

lc_prop_long <- merge(
  lc_counts, lc_totals,
  by = "grid_id" )

lc_prop_long$prop <- lc_prop_long$n_pixels / lc_prop_long$total_pixels

lc_n_types <- aggregate(
  lc_name ~ grid_id,
  data = lc_prop_long,
  FUN = function(x) length(unique(x))
)

names(lc_n_types)[2] <- "lc_n_types"

lc_combo <- aggregate(
  lc_name ~ grid_id,
  data = lc_prop_long,
  FUN = function(x) paste(sort(unique(x)), collapse = " + ")
)

names(lc_combo)[2] <- "lc_combo"

lc_prop_wide <- tidyr::pivot_wider(
  lc_prop_long[, c("grid_id", "lc_name", "prop")],
  names_from = lc_name,
  values_from = prop,
  values_fill = 0
)

lc_grid_summary <- merge( lc_n_types, lc_combo, by = "grid_id" )
lc_grid_summary <- merge( lc_grid_summary, lc_prop_wide, by = "grid_id" )

# =========================================================
# 5. Summarise land-cover combinations
# =========================================================
n_lc_combos <- length(unique(lc_grid_summary$lc_combo))

lc_combo_count <- as.data.frame(table(lc_grid_summary$lc_combo))
names(lc_combo_count) <- c("lc_combo", "n_fishnet")

lc_combo_count <- lc_combo_count[
  order(lc_combo_count$n_fishnet, decreasing = TRUE),
]

n_lc_combos
lc_combo_count

# =========================================================
# 6. Cluster grid cells by land-cover composition
# =========================================================
lc_prop_cols <- lc_n2c$lc_name
missing_cols <- setdiff(lc_prop_cols,names(lc_grid_summary))
for (col in missing_cols) {lc_grid_summary[[col]] <- 0}

lc_kmeans_data <- lc_grid_summary[, lc_prop_cols]

lc_kmeans_data[is.na(lc_kmeans_data)] <- 0

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

lc_grid_summary$lc_cluster <- lc_kmeans$cluster


lc_cluster_summary <- aggregate(
  lc_grid_summary[, lc_prop_cols],
  by = list(lc_cluster = lc_grid_summary$lc_cluster),
  FUN = mean
)

lc_cluster_summary


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


similarity <- 1 / (dist_to_centres + 1e-6)

cluster_membership <- similarity / rowSums(similarity)

cluster_membership <- as.data.frame(cluster_membership)

names(cluster_membership) <- paste0("mem_cluster_", 1:best_k)

cluster_membership$grid_id <- lc_grid_summary$grid_id

head(cluster_membership)

# =========================================================
# 8. Join bird data with land-cover cluster membership
# =========================================================
lc_grid_summary <- merge(
  lc_grid_summary,
  cluster_membership,
  by = "grid_id",
  all.x = TRUE
)

names(lc_grid_summary)

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

bird_fuzzy_data$SR_z <- as.numeric(scale(bird_fuzzy_data$SR))

bird_fuzzy_data$Abs_log <- log1p(bird_fuzzy_data$Abs)

bird_fuzzy_data$Abs_log_z <- as.numeric(scale(bird_fuzzy_data$Abs_log))

# Spatial block validation and stable cluster contribution


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

length(train_block_combinations)

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
# 11. Cluster-patch-level boundary and high-value adjustment
# =========================================================

# 1. Cells belonging to the same lc_cluster and touching each other
# 2. If a cluster patch contains high habitat scores, the whole patch
# 3. Cells on the boundary of a cluster patch are slightly weakened
# 4. The direction of adjustment is controlled by the SR contribution

hab_var <- "habitat_score"


# ---------------------------------------------------------
# 11.1 Prepare coordinate keys
# ---------------------------------------------------------

if (!all(c("cell_x", "cell_y") %in% names(bird_fuzzy_data))) {
  
  cell_xy <- st_coordinates(
    st_centroid(bird_fuzzy_data)
  )
  
  bird_fuzzy_data$cell_x <- cell_xy[, 1]
  bird_fuzzy_data$cell_y <- cell_xy[, 2]
}

bird_fuzzy_data$x_key <- round(
  bird_fuzzy_data$cell_x,
  6
)

bird_fuzzy_data$y_key <- round(
  bird_fuzzy_data$cell_y,
  6
)

bird_fuzzy_data$xy_key <- paste(
  bird_fuzzy_data$x_key,
  bird_fuzzy_data$y_key,
  sep = "_"
)

x_vals <- sort(
  unique(
    bird_fuzzy_data$x_key
  )
)

y_vals <- sort(
  unique(
    bird_fuzzy_data$y_key
  )
)

grid_dx <- median(
  diff(x_vals),
  na.rm = TRUE
)

grid_dy <- median(
  diff(y_vals),
  na.rm = TRUE
)

cell_lookup <- data.frame(
  row_id = seq_len(
    nrow(bird_fuzzy_data)
  ),
  xy_key = bird_fuzzy_data$xy_key
)


# ---------------------------------------------------------
# 11.2 Helper function: get same-cluster ratio
# ---------------------------------------------------------

get_same_cluster_ratio <- function(x_shift, y_shift) {
  
  ratio_out <- rep(
    NA_real_,
    nrow(bird_fuzzy_data)
  )
  
  n_out <- rep(
    NA_integer_,
    nrow(bird_fuzzy_data)
  )
  
  for (i in seq_len(nrow(bird_fuzzy_data))) {
    
    target_key <- paste(
      round(
        bird_fuzzy_data$x_key[i] + x_shift,
        6
      ),
      round(
        bird_fuzzy_data$y_key[i] + y_shift,
        6
      ),
      sep = "_"
    )
    
    neigh_id <- cell_lookup$row_id[
      match(
        target_key,
        cell_lookup$xy_key
      )
    ]
    
    neigh_id <- neigh_id[
      !is.na(neigh_id)
    ]
    
    if (length(neigh_id) > 0) {
      
      same_cluster <-
        bird_fuzzy_data$lc_cluster[neigh_id] ==
        bird_fuzzy_data$lc_cluster[i]
      
      n_out[i] <- length(neigh_id)
      
      ratio_out[i] <- mean(
        same_cluster,
        na.rm = TRUE
      )
    }
  }
  
  data.frame(
    n = n_out,
    ratio = ratio_out
  )
}


# ---------------------------------------------------------
# 11.3 Calculate rook same-cluster ratio
# ---------------------------------------------------------

rook_N <- get_same_cluster_ratio(0, grid_dy)
rook_S <- get_same_cluster_ratio(0, -grid_dy)
rook_W <- get_same_cluster_ratio(-grid_dx, 0)
rook_E <- get_same_cluster_ratio(grid_dx, 0)

rook_ratio_mat <- cbind(
  rook_N$ratio,
  rook_S$ratio,
  rook_W$ratio,
  rook_E$ratio
)

bird_fuzzy_data$n_rook_neighbours <- rowSums(
  !is.na(rook_ratio_mat)
)

bird_fuzzy_data$rook_same_cluster_ratio <- rowMeans(
  rook_ratio_mat,
  na.rm = TRUE
)

bird_fuzzy_data$rook_same_cluster_ratio[
  is.na(bird_fuzzy_data$rook_same_cluster_ratio)
] <- 0.5


bird_fuzzy_data$strict_rook_interior <-
  bird_fuzzy_data$n_rook_neighbours == 4 &
  bird_fuzzy_data$rook_same_cluster_ratio == 1

bird_fuzzy_data$edge_boundary_type <- ifelse(
  bird_fuzzy_data$strict_rook_interior,
  "Interior",
  "Boundary"
)

bird_fuzzy_data$edge_boundary_type <- as.factor(
  bird_fuzzy_data$edge_boundary_type
)

table(
  bird_fuzzy_data$edge_boundary_type
)


# ---------------------------------------------------------
# 11.4 Identify connected cluster patches
# ---------------------------------------------------------


touch_list_all <- st_touches(
  bird_fuzzy_data
)

cluster_patch_graph_edges <- data.frame(
  from = integer(),
  to = integer()
)

for (i in seq_len(nrow(bird_fuzzy_data))) {
  
  neigh_i <- touch_list_all[[i]]
  
  if (length(neigh_i) > 0) {
    
    same_cluster_neigh <- neigh_i[
      bird_fuzzy_data$lc_cluster[neigh_i] ==
        bird_fuzzy_data$lc_cluster[i]
    ]
    
    if (length(same_cluster_neigh) > 0) {
      
      cluster_patch_graph_edges <- rbind(
        cluster_patch_graph_edges,
        data.frame(
          from = i,
          to = same_cluster_neigh
        )
      )
    }
  }
}


if (nrow(cluster_patch_graph_edges) == 0) {
  
  bird_fuzzy_data$cluster_patch_id <- seq_len(
    nrow(bird_fuzzy_data)
  )
  
} else {
  
  cluster_patch_graph <- graph_from_data_frame(
    cluster_patch_graph_edges,
    directed = FALSE,
    vertices = data.frame(
      name = seq_len(
        nrow(bird_fuzzy_data)
      )
    )
  )
  
  cluster_patch_comp <- components(
    cluster_patch_graph
  )
  
  bird_fuzzy_data$cluster_patch_id <- cluster_patch_comp$membership[
    as.character(
      seq_len(nrow(bird_fuzzy_data))
    )
  ]
}

bird_fuzzy_data$cluster_patch_id <- as.integer(
  bird_fuzzy_data$cluster_patch_id
)

length(
  unique(
    bird_fuzzy_data$cluster_patch_id
  )
)


# ---------------------------------------------------------
# 11.5 Calculate cluster-patch-level habitat signal
# ---------------------------------------------------------

cluster_patch_summary <- bird_fuzzy_data %>%
  st_drop_geometry() %>%
  group_by(
    cluster_patch_id
  ) %>%
  summarise(
    lc_cluster = first(
      lc_cluster
    ),
    patch_n_cells = n(),
    patch_mean_habitat_score = mean(
      habitat_score,
      na.rm = TRUE
    ),
    patch_max_habitat_score = max(
      habitat_score,
      na.rm = TRUE
    ),
    patch_median_habitat_score = median(
      habitat_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

cluster_patch_summary


cluster_patch_high_threshold <- quantile(
  bird_fuzzy_data$habitat_score,
  probs = 0.80,
  na.rm = TRUE
)

cluster_patch_summary$patch_high_value <- ifelse(
  cluster_patch_summary$patch_max_habitat_score >=
    cluster_patch_high_threshold,
  1,
  0
)

table(
  cluster_patch_summary$patch_high_value
)


bird_fuzzy_data <- bird_fuzzy_data %>%
  select(
    -any_of(
      c(
        "patch_n_cells",
        "patch_mean_habitat_score",
        "patch_max_habitat_score",
        "patch_median_habitat_score",
        "patch_high_value"
      )
    )
  ) %>%
  left_join(
    cluster_patch_summary %>%
      select(
        cluster_patch_id,
        patch_n_cells,
        patch_mean_habitat_score,
        patch_max_habitat_score,
        patch_median_habitat_score,
        patch_high_value
      ),
    by = "cluster_patch_id"
  )


# ---------------------------------------------------------
# 11.6 Prepare signed cluster weights
# ---------------------------------------------------------

cluster_weight_table <- stable_SR_coef %>%
  mutate(
    lc_cluster = as.numeric(
      gsub(
        "mem_cluster_",
        "",
        cluster
      )
    ),
    cluster_weight_raw = mean_coef,
    cluster_weight = cluster_weight_raw /
      max(
        abs(cluster_weight_raw),
        na.rm = TRUE
      )
  ) %>%
  select(
    lc_cluster,
    cluster_weight_raw,
    cluster_weight
  )

cluster_weight_table

bird_fuzzy_data <- bird_fuzzy_data %>%
  select(
    -any_of(
      c(
        "cluster_weight_raw",
        "cluster_weight"
      )
    )
  ) %>%
  left_join(
    cluster_weight_table,
    by = "lc_cluster"
  )

bird_fuzzy_data$cluster_weight[
  is.na(bird_fuzzy_data$cluster_weight)
] <- 0


# ---------------------------------------------------------
# 11.7 Sensitivity test for patch and boundary adjustment
# ---------------------------------------------------------


patch_alpha_values <- c(
  0,
  0.02,
  0.05,
  0.08,
  0.10
)

boundary_alpha_values <- c(
  0,
  0.02,
  0.05,
  0.08,
  0.10
)

patch_boundary_tuning_summary <- data.frame()

for (patch_alpha_i in patch_alpha_values) {
  
  for (boundary_alpha_i in boundary_alpha_values) {
    
    
    patch_effect_i <- ifelse(
      bird_fuzzy_data$patch_high_value == 1,
      0.5,
      0
    )
    
    
    boundary_effect_i <- ifelse(
      bird_fuzzy_data$strict_rook_interior,
      0.5,
      -0.5
    )
    
    
    signed_patch_effect_i <-
      bird_fuzzy_data$cluster_weight *
      patch_effect_i
    
    signed_boundary_effect_i <-
      bird_fuzzy_data$cluster_weight *
      boundary_effect_i
    
    score_i <- bird_fuzzy_data[[hab_var]] *
      (
        1 +
          patch_alpha_i * signed_patch_effect_i +
          boundary_alpha_i * signed_boundary_effect_i
      )
    
    score_i <- pmax(
      0,
      pmin(
        1,
        score_i
      )
    )
    
    threshold_i <- quantile(
      score_i,
      probs = 0.80,
      na.rm = TRUE
    )
    
    habitat_i <- ifelse(
      score_i >= threshold_i,
      1,
      0
    )
    
    cell_share_i <- mean(
      habitat_i == 1,
      na.rm = TRUE
    )
    
    abs_share_i <- sum(
      bird_fuzzy_data$Abs[
        habitat_i == 1
      ],
      na.rm = TRUE
    ) /
      sum(
        bird_fuzzy_data$Abs,
        na.rm = TRUE
      )
    
    patch_boundary_tuning_summary <- rbind(
      patch_boundary_tuning_summary,
      data.frame(
        patch_alpha = patch_alpha_i,
        boundary_alpha = boundary_alpha_i,
        threshold_prob = 0.80,
        n_habitat_cells = sum(
          habitat_i == 1,
          na.rm = TRUE
        ),
        cell_share = cell_share_i,
        Abs_share = abs_share_i,
        Abs_enrichment = abs_share_i / cell_share_i
      )
    )
  }
}

patch_boundary_tuning_summary <- patch_boundary_tuning_summary %>%
  arrange(
    desc(
      Abs_enrichment
    )
  )

patch_boundary_tuning_summary


ggplot(
  patch_boundary_tuning_summary,
  aes(
    x = patch_alpha,
    y = Abs_enrichment,
    colour = as.factor(boundary_alpha),
    group = boundary_alpha
  )
) +
  geom_line() +
  geom_point(
    size = 2
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    title = "Tuning cluster-patch and boundary adjustment",
    x = "Patch alpha",
    y = "Abs enrichment",
    colour = "Boundary alpha"
  )


# ---------------------------------------------------------
# 11.8 Apply manually selected patch-boundary adjustment
# ---------------------------------------------------------


final_patch_alpha <- 0.1
final_boundary_alpha <- 0.1


selected_patch_boundary_setting <- patch_boundary_tuning_summary %>%
  filter(
    patch_alpha == final_patch_alpha,
    boundary_alpha == final_boundary_alpha
  )

selected_patch_boundary_setting


bird_fuzzy_data$patch_high_value_effect <- ifelse(
  bird_fuzzy_data$patch_high_value == 1,
  0.5,
  0
)


bird_fuzzy_data$edge_boundary_effect <- ifelse(
  bird_fuzzy_data$strict_rook_interior,
  0.5,
  -0.5
)


bird_fuzzy_data$signed_patch_effect <-
  bird_fuzzy_data$cluster_weight *
  bird_fuzzy_data$patch_high_value_effect

bird_fuzzy_data$signed_edge_effect <-
  bird_fuzzy_data$cluster_weight *
  bird_fuzzy_data$edge_boundary_effect

bird_fuzzy_data$signed_patch_effect[
  is.na(bird_fuzzy_data$signed_patch_effect)
] <- 0

bird_fuzzy_data$signed_edge_effect[
  is.na(bird_fuzzy_data$signed_edge_effect)
] <- 0


bird_fuzzy_data$combined_adjustment <-
  final_patch_alpha * bird_fuzzy_data$signed_patch_effect +
  final_boundary_alpha * bird_fuzzy_data$signed_edge_effect

bird_fuzzy_data$habitat_score_context <-
  bird_fuzzy_data[[hab_var]] *
  (
    1 + bird_fuzzy_data$combined_adjustment
  )

bird_fuzzy_data$habitat_score_context <- pmax(
  0,
  pmin(
    1,
    bird_fuzzy_data$habitat_score_context
  )
)


# ---------------------------------------------------------
# 11.9 Classify and check final adjustment
# ---------------------------------------------------------

bird_fuzzy_data$edge_boundary_type <- ifelse(
  bird_fuzzy_data$strict_rook_interior,
  "Interior",
  "Boundary"
)

bird_fuzzy_data$edge_boundary_type <- as.factor(
  bird_fuzzy_data$edge_boundary_type
)

bird_fuzzy_data$patch_value_type <- ifelse(
  bird_fuzzy_data$patch_high_value == 1,
  "High-value cluster patch",
  "Other cluster patch"
)

bird_fuzzy_data$patch_value_type <- as.factor(
  bird_fuzzy_data$patch_value_type
)

bird_fuzzy_data$combined_adjustment_type <- ifelse(
  bird_fuzzy_data$combined_adjustment > 0,
  "Strengthened",
  ifelse(
    bird_fuzzy_data$combined_adjustment < 0,
    "Weakened",
    "Neutral"
  )
)

bird_fuzzy_data$combined_adjustment_type <- as.factor(
  bird_fuzzy_data$combined_adjustment_type
)

summary(
  bird_fuzzy_data[, c(
    "habitat_score",
    "habitat_score_context",
    "lc_cluster",
    "cluster_patch_id",
    "patch_n_cells",
    "patch_max_habitat_score",
    "patch_high_value",
    "cluster_weight",
    "strict_rook_interior",
    "signed_patch_effect",
    "signed_edge_effect",
    "combined_adjustment"
  )]
)

table(
  bird_fuzzy_data$edge_boundary_type
)

table(
  bird_fuzzy_data$patch_value_type
)

table(
  bird_fuzzy_data$combined_adjustment_type
)


ggplot(
  st_drop_geometry(bird_fuzzy_data),
  aes(
    x = habitat_score,
    y = habitat_score_context,
    colour = combined_adjustment_type
  )
) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    title = "Cluster-patch-level habitat adjustment",
    x = "Original habitat score",
    y = "Patch-adjusted habitat score",
    colour = "Adjustment"
  )


# =========================================================
# 12. Define potential habitat cells
# =========================================================


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
      n_habitat_cells = sum(
        habitat_i == 1,
        na.rm = TRUE
      ),
      prop_habitat_cells = mean(
        habitat_i == 1,
        na.rm = TRUE
      )
    )
  )
}

threshold_summary


hab_threshold_prob <- 0.80

hab_threshold <- quantile(
  bird_fuzzy_data$habitat_score_context,
  probs = hab_threshold_prob,
  na.rm = TRUE
)

hab_threshold


bird_fuzzy_data$potential_habitat <- ifelse(
  bird_fuzzy_data$habitat_score_context >= hab_threshold,
  1,
  0
)

bird_fuzzy_data$potential_habitat <- as.factor(
  bird_fuzzy_data$potential_habitat
)


plot(
  bird_fuzzy_data["potential_habitat"],
  main = "Potential habitat cells"
)

# =========================================================
# 13. Select potential habitat cells and create patches
# =========================================================

hab_cells <- bird_fuzzy_data %>%
  filter(
    potential_habitat == 1
  )

if (nrow(hab_cells) == 0) {
  
  stop("No habitat cells found. Try a lower threshold.")
  
}


if (nrow(hab_cells) == 1) {
  
  hab_cells$patch_id <- 1
  
} else {
  
  touch_list_hab <- st_touches(
    hab_cells
  )
  
  hab_graph <- graph_from_adj_list(
    touch_list_hab,
    mode = "all"
  )
  
  hab_comp <- components(
    hab_graph
  )
  
  hab_cells$patch_id <- hab_comp$membership
}


hab_patches <- hab_cells %>%
  group_by(
    patch_id
  ) %>%
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
    mean_edge_boundary_effect = mean(
      edge_boundary_effect,
      na.rm = TRUE
    ),
    mean_signed_edge_effect = mean(
      signed_edge_effect,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


hab_patches$patch_area_m2 <- as.numeric(
  st_area(
    hab_patches
  )
)

hab_patches$patch_area_km2 <-
  hab_patches$patch_area_m2 / 1e6

hab_patches$patch_perimeter_m <- as.numeric(
  st_length(
    st_boundary(
      hab_patches
    )
  )
)

hab_patches$patch_perimeter_km <-
  hab_patches$patch_perimeter_m / 1000


hab_patches$shape_complexity <-
  hab_patches$patch_perimeter_km /
  (
    2 * sqrt(
      pi * hab_patches$patch_area_km2
    )
  )


patch_summary <- data.frame(
  n_habitat_cells = nrow(
    hab_cells
  ),
  n_patches = nrow(
    hab_patches
  ),
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
  ),
  mean_patch_edge_effect = mean(
    hab_patches$mean_edge_boundary_effect,
    na.rm = TRUE
  ),
  mean_patch_signed_effect = mean(
    hab_patches$mean_signed_edge_effect,
    na.rm = TRUE
  )
)

patch_summary


plot(
  hab_patches["patch_id"],
  main = "Potential habitat patches"
)

plot(
  hab_patches["mean_habitat_score"],
  main = "Mean habitat score by patch"
)


largest_patches <- hab_patches %>%
  st_drop_geometry() %>%
  arrange(
    desc(
      patch_area_km2
    )
  ) %>%
  select(
    patch_id,
    n_cells,
    patch_area_km2,
    mean_habitat_score,
    max_habitat_score,
    mean_edge_boundary_effect,
    mean_signed_edge_effect,
    shape_complexity
  )

largest_patches


# =========================================================
# 15. Validate whether potential habitat captures Abs
# =========================================================

abs_capture_summary <- bird_fuzzy_data %>%
  st_drop_geometry() %>%
  group_by(
    potential_habitat
  ) %>%
  summarise(
    n_cells = n(),
    total_Abs = sum(
      Abs,
      na.rm = TRUE
    ),
    mean_Abs = mean(
      Abs,
      na.rm = TRUE
    ),
    median_Abs = median(
      Abs,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

abs_capture_summary$total_cells <- sum(
  abs_capture_summary$n_cells
)

abs_capture_summary$total_Abs_all <- sum(
  abs_capture_summary$total_Abs,
  na.rm = TRUE
)

abs_capture_summary$cell_share <-
  abs_capture_summary$n_cells /
  abs_capture_summary$total_cells

abs_capture_summary$Abs_share <-
  abs_capture_summary$total_Abs /
  abs_capture_summary$total_Abs_all

abs_capture_summary$Abs_enrichment <-
  abs_capture_summary$Abs_share /
  abs_capture_summary$cell_share

abs_capture_summary


abs_capture_plot <- abs_capture_summary %>%
  filter(
    potential_habitat == 1
  )

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
  aes( x = metric, y = value )) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Abs captured by potential habitat cells",
    x = "",
    y = "Proportion"
  )

boundary_response_summary <- bird_fuzzy_data %>%
  st_drop_geometry() %>%
  group_by(edge_boundary_type) %>%
  summarise(
    n_cells = n(),
    mean_SR = mean(SR, na.rm = TRUE),
    mean_Abs = mean(Abs, na.rm = TRUE),
    median_Abs = median(Abs, na.rm = TRUE),
    mean_habitat_score = mean(habitat_score, na.rm = TRUE),
    mean_habitat_score_context = mean(habitat_score_context, na.rm = TRUE),
    .groups = "drop"
  )

boundary_response_summary
