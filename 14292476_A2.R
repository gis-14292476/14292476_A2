
setwd("F:/14292476_A2")
set.seed(33)


library(sf)
library(terra)
library(dplyr)
library(tidyr)





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

# Land-cover class table
lc_n2c <- data.frame(
  lc_class = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  lc_name = c("Woodland","Conifer","Arable",
              "Grassland","Moor","Bog",
              "Water","Saltmarsh","Urban")
  )
lc_extract <- merge(lc_extract,lc_n2c,by = "lc_class",all.x = TRUE)






# Calculate land-cover proportions by grid cell

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






# Cluster fishnet cells by land-cover proportions
lc_prop_cols <- lc_n2c$lc_name
missing_cols <- setdiff(lc_prop_cols, names(lc_grid_summary))

for (col in missing_cols) {lc_grid_summary[[col]] <- 0}

# Prepare k-means data
lc_kmeans_data <- lc_grid_summary[, lc_prop_cols]
lc_kmeans_data[is.na(lc_kmeans_data)] <- 0
rownames(lc_kmeans_data) <- lc_grid_summary$grid_id


# Elbow method
max_k <- 20
wss <- numeric(max_k)

for (k in 1:max_k) {
  km <- kmeans(lc_kmeans_data, centers = k, nstart = 100)
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

best_k <- 7

lc_kmeans <- kmeans(
  lc_kmeans_data,
  centers = best_k,
  nstart = 100
)

lc_grid_summary$lc_cluster <- lc_kmeans$cluster


# Summarise clusters

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