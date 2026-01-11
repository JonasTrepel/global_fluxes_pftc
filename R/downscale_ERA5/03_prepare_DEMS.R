library(sf)
library(tidyverse)
library(terra)
library(rstac)

p <- st_read("data/spatial_data/area_polygons.gpkg")
load(file = "data/spatial_data/utm_zones.rda")

dl_dir <- "/scratch/project_2003061/ERA5/PFTC"

for(area_id in p$area){
  print(area_id)
  
  pt <- p %>% filter(area == area_id)
  aoi_mid <- pt %>%
    st_centroid() %>%
    st_transform(crs = 4326)
  
  utmall <- utm_zones
  
  # WGS84 UTM zones to set the correct projection
  utm <- utmall[aoi_mid,] # Which zone the study points falls in
  
  if (nrow(utm) == 0) {
    utm <- utmall[st_nearest_feature(aoi_mid, utmall),]
  }
  
  lat <- st_coordinates(aoi_mid)[,"Y"]
  utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0", utm$ZONE), utm$ZONE)
  epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
  
  aoi <- pt %>% 
    st_transform(epsg)
  
  area_dl_dir <- paste0(dl_dir, "/", area_id, "/")
  
  if(!dir.exists(area_dl_dir)){
    dir.create(area_dl_dir)
  }
  
  ############ Fine-scale DEM
  # ALOS DEM from STAC
  s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")
  it_obj <- s_obj %>%
    stac_search(collections = "alos-dem",
                bbox = st_bbox(pt %>% st_transform(4326)),
                limit = 1000) %>%
    get_request()
  
  # Internal Function to Create a Base URL with Microsoft Planetary Computer
  make_vsicurl_url_dem <- function(base_url) {
    paste0(
      "/vsicurl", 
      "?pc_url_signing=yes",
      "&pc_collection=alos-dem",
      "&url=",
      base_url
    )
  }
  
  juuh <- lapply(it_obj$features, function(ft) {
    full_url <- make_vsicurl_url_dem(assets_url(ft) %>% sort)
    full_url <- full_url[endsWith(full_url, "_DSM.tif")]
    file_names <- gsub("TIF$", "tif", basename(full_url))
    
    juuh <- lapply(seq_len(length(full_url)), function(nr) {
      e <- try({
        gdal_utils(
          "warp",
          source = full_url[[nr]],
          destination = paste0(area_dl_dir, "/", file_names[[nr]]),
          options = c(
            "-t_srs", st_crs(aoi)$wkt,
            "-te", st_bbox(aoi),
            "-tr", c(30, 30)
          )
        )
      }, silent = TRUE)
      if (class(e)[[1]] == "try-error") {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
  })
  
  dems <- list.files(area_dl_dir, pattern = "_DSM.tif", full.names = TRUE)
  dems <- lapply(dems, function(x) {
    dem <- rast(x)
    dem[dem == 0] <- NA
    return(dem)
  })
  
  dem <- sprc(dems)
  dem <- mosaic(dem)
  dem[dem < 0] <- 0
  dem[is.na(dem)] <- 0
  plot(dem, main = area_id)
  
  writeRaster(dem, paste0(area_dl_dir, "/ALOSDEM.tif"), overwrite = TRUE)
  unlink(list.files(area_dl_dir, pattern = "_DSM.tif", full.names = TRUE))
  
  # ERA5-Land DEM
  
  dem <- rast("/scratch/project_2003061/ERA5/ERA5_LAND_DEM.nc") %>% project(rast(list.files(area_dl_dir, pattern = "ERA5_COMB", full.names = TRUE)[[1]]), method = "bilinear")
  lsm <- rast("/scratch/project_2003061/ERA5/ERA5_LAND_water.nc") %>% project(rast(list.files(area_dl_dir, pattern = "ERA5_COMB", full.names = TRUE)[[1]]), method = "bilinear")
  
  writeRaster(dem, paste0(area_dl_dir, "/ERA5-Land_DEM.tif"), overwrite = TRUE)
  writeRaster(lsm, paste0(area_dl_dir, "/ERA5-Land_water.tif"), overwrite = TRUE)
  
}
