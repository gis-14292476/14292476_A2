
setwd("F:/14292476_A2")
set.seed(33)


library(sf)
library(terra)
library(dplyr)
library(tidyr)

# read data
birds_grid <- st_read("./data/Birds/GM_Birds_2025.shp")
birds_lcm  <- rast("./data/Birds/gm_lcm_2022.tif")

# data preprocess
# geometry correction
birds_grid <- st_make_valid(birds_grid)
# create grid ID
birds_grid$grid_id <- seq_len(nrow(birds_grid))

# check CRS
crs(birds_grid) == crs(birds_lcm)

# clip land-cover
birds_vect <- vect(birds_grid)

birds_lcm_clip <- crop(birds_lcm, birds_vect)
birds_lcm_mask <- mask(birds_lcm_clip, birds_vect)

# extract land-cover in each grid
lc_extract <- terra::extract( birds_lcm_mask, birds_vect )
lc_extract <- lc_extract[, 1:2]
names(lc_extract) <- c("grid_id", "lc_class")
# delete 0 = NA
lc_extract <- subset(lc_extract,!is.na(lc_class) & lc_class != 0)

# land-cover class table
lc_n2c <- data.frame(
  lc_class = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  lc_name = c("Woodland","Conifer","Arable",
              "Grassland","Moor","Bog",
              "Water","Saltmarsh","Urban")
  )
lc_extract <- merge(lc_extract,lc_n2c,by = "lc_class",all.x = TRUE)

