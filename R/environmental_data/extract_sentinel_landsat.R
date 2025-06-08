library(rgee)
library(tidyverse)
library(sf)

ee_Initialize(drive = "T")

ee_check()

p <- st_read("data/coordinates.gpkg")

for(i in unique(p$area)){
  
}
