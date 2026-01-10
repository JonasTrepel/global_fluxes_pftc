library(tidyverse)
library(terra)
library(sf)

aoi <- st_read("data/processed_data/preliminary_data/prelim_flux_loc.gpkg") %>% 
  mutate(tier2 = gsub("South_Africa","SouthAfrica",tier)) %>% 
  mutate(area = unlist(lapply(tier2, function(x){ str_split_i(x, "_", 1)}))) %>%
  st_transform(crs = 4326) %>% 
  mutate(name = as.character(area))

dl_dir <- "/scratch/project_2003061/ERA5/PFTC"

all <- tibble()
for(i in unique(aoi$area)){
  print(i)
  ii <- i
  if(i == "Colorado"){
    ii <- "USA"
  }
  if(i == "SouthAfrica"){
    ii <- "south_africa"
  }
  
  # Temperature
  fl <- list.files(paste0(dl_dir,"/",ii), pattern = "_T2m.tif$", full.names = TRUE)
  
  pt <- aoi %>% 
    filter(area == i) %>% 
    st_transform(crs = st_crs(rast(fl[[1]]))) %>% 
    select(plot_id) %>% 
    distinct()
  
  temps <- terra::extract(rast(fl), pt)
  era_times <- names(temps)[-1]
  era_times <- ifelse(grepl(" ", era_times), era_times, paste0(era_times, " 00:00:00"))
  
  temps <- t(temps)
  temps <- temps[-1,]
  colnames(temps) <- pt$plot_id
  
  temps <- temps %>% 
    as.data.frame() %>% 
    mutate(datetime = ymd_hms(era_times)) %>% 
    pivot_longer(cols = c(-datetime), names_to = "plot_id", values_to = "T2m")
  
  
  # Relative humidity
  fl <- list.files(paste0(dl_dir,"/",ii), pattern = "_relhum.tif$", full.names = TRUE)
  
  relhum <- terra::extract(rast(fl), pt)
  era_times <- names(relhum)[-1]
  era_times <- ifelse(grepl(" ", era_times), era_times, paste0(era_times, " 00:00:00"))
  
  relhum <- t(relhum)
  relhum <- relhum[-1,]
  colnames(relhum) <- pt$plot_id
  
  relhum <- relhum %>% 
    as.data.frame() %>% 
    mutate(datetime = ymd_hms(era_times)) %>% 
    pivot_longer(cols = c(-datetime), names_to = "plot_id", values_to = "relhum")
  
  # Relative humidity
  fl <- list.files(paste0(dl_dir,"/",ii), pattern = "_windspeed.tif$", full.names = TRUE)
  
  windspeed <- terra::extract(rast(fl), pt)
  era_times <- names(windspeed)[-1]
  era_times <- ifelse(grepl(" ", era_times), era_times, paste0(era_times, " 00:00:00"))
  
  windspeed <- t(windspeed)
  windspeed <- windspeed[-1,]
  colnames(windspeed) <- pt$plot_id
  
  windspeed <- windspeed %>% 
    as.data.frame() %>% 
    mutate(datetime = ymd_hms(era_times)) %>% 
    pivot_longer(cols = c(-datetime), names_to = "plot_id", values_to = "windspeed")
  
  # Join
  
  all <- bind_rows(all,
                   full_join(temps, relhum) %>% 
                     full_join(., windspeed) %>% 
                     mutate(area = i,
                            T2m = T2m/10,
                            relhum = relhum/10,
                            windspeed = windspeed/100) %>% 
                     distinct() %>% 
                     relocate(area, plot_id))
  
}

dd <- all %>% 
  mutate(date = as_date(datetime)) %>% 
  group_by(area, plot_id, date) %>% 
  summarise(across(T2m:windspeed, mean), .groups = "drop")

# Calculate VPD
dd <- dd %>%
  mutate(
    # Calculate saturation vapor pressure (es) using Tetens equation
    es = 0.6108 * exp((17.27 * T2m) / (T2m + 237.3)),
    
    # Calculate actual vapor pressure (ea)
    ea = (es * relhum) / 100,
    
    # Calculate Vapor Pressure Deficit (VPD)
    VPD = es - ea
  ) %>% 
  select(-es, -ea)

dm <- dd %>% 
  mutate(m = month(date)) %>% 
  group_by(area, m) %>% 
  summarise(meanT = mean(T2m, na.rm = TRUE))

dm <- dm %>% 
  filter(!(area == "Peru" & m %in% c(1, 5:12)))

dm <- dm %>% 
  slice_max(meanT, n = 3) %>% as.data.frame()

d3m <- right_join(dd %>% mutate(m = month(date)), 
          dm) %>% 
  group_by(area, plot_id) %>% 
  summarise(across(T2m:VPD, mean)) %>% 
  ungroup

summary(d3m)

d3m %>% 
  group_by(area) %>% 
  summarise(across(T2m:VPD, max))

cor(d3m$VPD, d3m$relhum)

write_csv(d3m, "output/downscaled_climate.csv")

aoi %>% 
  group_by(tier) %>% 
  summarise(start_date = min(date),
            end_date = max(date))

