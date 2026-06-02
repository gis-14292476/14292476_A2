
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

lc_grid_summary


