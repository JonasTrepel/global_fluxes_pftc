library(sf)
library(tidyverse)

d1 <- read_csv("data/spatial_data/pftc_coordinates.csv")
p <- st_read("data/spatial_data/coordinates.gpkg")
p %>% print(n = 222)

# Step 1: Create a convex hull for all points within each area
hulls <- p %>%
  group_by(area) %>%
  summarise(geometry = st_convex_hull(st_union(geom))) %>%
  st_as_sf()

# Step 2: Create a 5km buffer around each convex hull
buffered_hulls <- st_buffer(hulls, dist = 5000)

# Step 3: Create a square (bounding box) for each buffered hull
square_polygons <- buffered_hulls %>%
  group_by(area) %>%
  mutate(geometry = st_as_sfc(st_bbox(geometry)))

plot(square_polygons)

points_sf <- st_as_sf(d1, coords = c("longitude_e", "latitude_n"), crs = st_crs(square_polygons))

# Check which points are within any square_polygons
points_within <- points_sf %>%
  mutate(
    inside_any = st_intersects(geometry, square_polygons, sparse = FALSE) %>%
      apply(1, any)
  )

all(points_within$inside_any)
table(points_within$inside_any)

st_write(square_polygons, "data/spatial_data/area_polygons.gpkg")
