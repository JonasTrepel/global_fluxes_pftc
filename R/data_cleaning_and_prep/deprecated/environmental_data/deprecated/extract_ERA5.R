library(tidyverse)
library(reticulate)
library(rgee, lib.loc = "/projappl/project_2003061/Rpackages/")
library(googledrive)
library(parallel)
library(sf)

ee_check()
# ee_install_upgrade()

ee_Initialize("pekkaniittynen", drive = T)

ee$Authenticate(auth_mode='notebook')
ee$Initialize(project='radiant-planet-412111')
ee$String('Hello from the Earth Engine servers!')$getInfo()

d <- st_read("data/coordinates.gpkg")

d$Lat <- st_coordinates(d)[,"Y"]
d$Lon <- st_coordinates(d)[,"X"]

co <- d %>% 
  st_drop_geometry() %>% 
  mutate(area_id = paste0(area, "_", rownames(.)))

sdate <- "2014-01-01"
edate <- "2023-12-31"

df <- lapply(co$area_id, function(i){
  # i <- "Svalbard_177"
  print(i)
  aoii <- co %>% 
    filter(area_id == i) %>% 
    st_as_sf(coords = c("Lon","Lat"), crs = 4326)
  
  aoi_ee <- aoii %>% st_transform(crs = 4326) %>% 
    st_geometry() %>% 
    sf_as_ee()
  
  dataset <- ee$ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')$filterDate(sdate, as.character(ymd(edate) + days(1)))
  dataset <- dataset$filterBounds(geometry = aoi_ee)
  dataset <- dataset$select(c("temperature_2m"))
  dataset <- dataset$toBands()
  # di <- dataset$getInfo()
  # lapply(di$bands, function(x) x$id) %>% unlist
  
  valuesList <- dataset$reduceRegion(
    reducer= ee$Reducer$toList() ,
    geometry=  aoi_ee
  )$values()
  
  temps <- valuesList$getInfo() %>% unlist
  
  if(!is.null(temps)){
    df <- tibble(area_id = i,
                 date = seq.Date(ymd(sdate), ymd(edate), by = "day"),
                 temperature_2m = temps - 273.15)
  } else {
    df <- tibble(area_id = i,
                 date = seq.Date(ymd(sdate), ymd(edate), by = "day"),
                 temperature_2m = NA)
  }
  
  
  # df %>% 
  #   ggplot(aes(y = temperature_2m, x = date)) +
  #   geom_line()
  
  return(df)
}) %>% bind_rows()

df %>%
  drop_na() %>% 
  ggplot(aes(y = temperature_2m, x = date, group = area_id, color = area_id)) +
  geom_line() +
  theme(legend.position = "none")

df %>% write_csv("data/environmental_data/ERA5_temps.csv")
