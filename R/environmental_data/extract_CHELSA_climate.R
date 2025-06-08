library(terra)
library(sf)
library(tidyverse)

p <- st_read("data/coordinates.gpkg")

# MEAN TEMPERATURE 

dl_dir <- "/scratch/project_2007415/temp/"

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          1,"_1981-2010_V.2.1.tif")

download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, basename(tif_to_download)))

pt <- p %>% 
  st_transform(crs(r, proj = T)) %>% 
  st_buffer(10)

res <- terra::extract(r, pt, fun = mean, weights=TRUE, na.rm = TRUE)
names(res)[2] <- "annual_mean_temperature"

d1 <- bind_cols(pt %>% st_drop_geometry(), res %>% select(-ID)) %>% 
  mutate(annual_mean_temperature = round(annual_mean_temperature,2))

d1 %>% filter(is.na(annual_mean_temperature))

# DIURNAL TEMPERATURE RANGE

dl_dir <- "/scratch/project_2007415/temp/"

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          2,"_1981-2010_V.2.1.tif")

download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, basename(tif_to_download)))

pt <- p %>% 
  st_transform(crs(r, proj = T)) %>% 
  st_buffer(10)

res <- terra::extract(r, pt, fun = mean, weights=TRUE, na.rm = TRUE)
names(res)[2] <- "diurnal_temperature_range"

d2 <- bind_cols(pt %>% st_drop_geometry(), res %>% select(-ID)) %>% 
  mutate(diurnal_temperature_range = round(diurnal_temperature_range,2))

d2 %>% filter(is.na(diurnal_temperature_range))

# Temperature seasonality

dl_dir <- "/scratch/project_2007415/temp/"

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          4,"_1981-2010_V.2.1.tif")

download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, basename(tif_to_download)))

pt <- p %>% 
  st_transform(crs(r, proj = T)) %>% 
  st_buffer(10)

res <- terra::extract(r, pt, fun = mean, weights=TRUE, na.rm = TRUE)
names(res)[2] <- "temperature_seasonality"

d3 <- bind_cols(pt %>% st_drop_geometry(), res %>% select(-ID)) %>% 
  mutate(temperature_seasonality = round(temperature_seasonality/100,2))

d3 %>% filter(is.na(temperature_seasonality))

# Annual precipitation

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          12,"_1981-2010_V.2.1.tif")

download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, basename(tif_to_download)))

pt <- p %>% 
  st_transform(crs(r, proj = T)) %>% 
  st_buffer(10)

res <- terra::extract(r, pt, fun = mean, weights=TRUE, na.rm = TRUE)
names(res)[2] <- "annual_precipitation"

d4 <- bind_cols(pt %>% st_drop_geometry(), res %>% select(-ID)) %>% 
  mutate(annual_precipitation = round(annual_precipitation,1))

d4 %>% filter(is.na(annual_precipitation))

# precipitation seasonality

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",
                          15,"_1981-2010_V.2.1.tif")

download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, basename(tif_to_download)))

pt <- p %>% 
  st_transform(crs(r, proj = T)) %>% 
  st_buffer(10)

res <- terra::extract(r, pt, fun = mean, weights=TRUE, na.rm = TRUE)
names(res)[2] <- "precipitation_seasonality"

d5 <- bind_cols(pt %>% st_drop_geometry(), res %>% select(-ID)) %>% 
  mutate(precipitation_seasonality = round(precipitation_seasonality,1))

d5 %>% filter(is.na(precipitation_seasonality))

# Köppen Geiger climate classification

tif_to_download <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_kg2_1981-2010_V.2.1.tif")

download.file(url = tif_to_download, destfile = paste0(dl_dir, basename(tif_to_download)))

r <- rast(paste0(dl_dir, basename(tif_to_download)))

pt <- p %>% 
  st_transform(crs(r, proj = T))

res <- terra::extract(r, pt)
names(res)[2] <- "KG_climate_classification"

d6 <- bind_cols(pt %>% st_drop_geometry(), res %>% select(-ID))

d6 %>% filter(is.na(KG_climate_classification))

d <- full_join(d1,d2) %>% 
  full_join(., d3) %>% 
  full_join(., d4) %>% 
  full_join(., d5) %>% 
  full_join(., d6)

write_csv(d, "data/environmental_data/climate.csv")

