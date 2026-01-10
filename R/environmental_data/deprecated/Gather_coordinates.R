library(tidyverse)
library(sf)

# China

p_china <- tibble(area = "China",
                  site = c("High Elevation","Alpine","Middle","Lowland"),
                  elev = c(4100, 3850, 3500, 3000),
                  lat = c(29.90742, 29.88911, 29.86192, 29.84347),
                  lon = c(102.0118, 102.0173, 102.0360, 102.0343)) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326)

# USA, Colorado

p_usa <- read_csv("data/raw_data/colorado/rmbl_site_info.csv") %>% 
  select(-lat_long) %>% 
  rename(lon = long) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>% 
  mutate(area = "USA") %>% 
  relocate(area)

# Peru

p_peru <- read_csv("data/raw_data/peru/PU.10_PFTC3.10_2020_Peru_Coordinates.csv") %>% 
  select(-Comment, -Burn_year) %>% 
  rename(lon = Longitude, lat = Latitude, elev = Elevation, 
         site = Site, treatment = Treatment, plot = PlotID) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 4326) %>% 
  mutate(area = "Peru") %>% 
  relocate(area)

# South Africa

p_sa <- st_read("data/raw_data/south_africa/PFTC7_plot_coordinates.gpkg") %>% 
  select(SiteID, aspect, PlotID, elevation_true) %>% 
  rename(elev = elevation_true,
         site = SiteID, plot = PlotID,
         treatment = aspect,
         geometry = geom) %>% 
  mutate(area = "south_africa") %>% 
  relocate(area)

# Svalbard

p_svalbard <- bind_rows(read_csv("data/raw_data/svalbard/PFTC4_Svalbard_Coordinates_ITEX.csv") %>% 
                          select(Treatment:Longitude_E) %>% 
                          drop_na() %>% 
                          rename(treatment = Treatment, site = Site, 
                                 elev = Elevation_m, lat = Latitude_N, lon = Longitude_E) %>% 
                          st_as_sf(coords = c("lon","lat"), crs = 4326) %>% 
                          mutate(area = "Svalbard") %>% 
                          relocate(area),
                        read_csv("data/raw_data/svalbard/PFTC4_Svalbard_Coordinates_Gradient.csv") %>% 
                          rename(treatment = Gradient, site = Site, plot = PlotID,
                                 elev = Elevation_m, lat = Latitude_N, lon = Longitude_E) %>% 
                          st_as_sf(coords = c("lon","lat"), crs = 4326) %>% 
                          mutate(area = "Svalbard") %>% 
                          relocate(area) %>% 
                          mutate(site = as.character(site)))
  



# Norway
d <- read_csv("data/processed_data/preliminary_data/prelim_fluxes.csv") %>% 
  filter(tier == "Norway_2022") %>% 
  group_by(site, plot_id) %>% 
  summarise(elevation = mean(elevation),
            latitude = mean(latitude),
            longitude = mean(longitude)) %>% 
  ungroup

p_norway <- d %>% 
  st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>% 
  rename(plot = plot_id,
         elev = elevation) %>% 
  mutate(area = "Norway")

# Combine

p_all <- bind_rows(p_china, 
                   p_usa,
                   p_peru %>% mutate(plot = as.character(plot)),
                   p_sa,
                   p_svalbard,
                   p_norway)

st_write(p_all, "data/coordinates.gpkg", append = FALSE)

