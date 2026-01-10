library(sf)
library(tidyverse)
library(mcera5)
library(ecmwfr)
library(ncdf4)
library(keyring)
library(terra)

p <- st_read("data/spatial_data/area_polygons.gpkg") %>% 
  left_join(.,
            tibble(area = c("China","Norway","Peru","Svalbard","USA","south_africa"),
                   year = c("2015","2022","2018","2018","2017","2023")))

temp_dir <- "/scratch/project_2003061/temp/"
dl_dir <- "/scratch/project_2003061/ERA5/PFTC"

# Set here your ecmwf credentials
uid <- ""
cds_access_token <- ""

ecmwfr::wf_set_key(user = uid,
                   key = cds_access_token)

for(area_id in p$area){
  # area_id <- "Svalbard"
  print(area_id)
  pt <- p %>% filter(area == area_id)
  
  area_dl_dir <- paste0(dl_dir, "/", area_id, "/")
  
  if(!dir.exists(area_dl_dir)){
    dir.create(area_dl_dir)
  }
  
  # bounding coordinates (in WGS84 / EPSG:4326) 
  xmn <- st_bbox(pt)$xmin
  xmx <- st_bbox(pt)$xmax
  ymn <- st_bbox(pt)$ymin
  ymx <- st_bbox(pt)$ymax
  
  first_date <- paste0(pt$year, "-01-01")
  last_date <- paste0(pt$year, "-12-31")
  
  dates <- tibble(date = seq.Date(ymd(first_date), ymd(last_date), by = "days")) %>% 
    mutate(y = year(date),
           m = month(date)) %>% 
    group_by(y, m) %>% 
    summarise(st_date = min(date),
              en_date = max(date)) %>% 
    ungroup()
  
  for(i in seq_len(nrow(dates))){
    print(dates %>% slice(i) %>% pull(st_date))
    
    st_time <- as.POSIXlt(paste0(dates %>% slice(i) %>% pull(st_date), " 00:00"), tz = "UTC")
    en_time <- as.POSIXlt(paste0(dates %>% slice(i) %>% pull(en_date), " 00:00"), tz = "UTC")
    
    if(!file.exists(paste0(area_dl_dir, "/ERA5_SST_", year(st_time), "_", month(st_time), ".nc"))){
      # temporal extent
      
      reql <- build_era5_land_request(xmin = xmn, xmax = xmx,
                                      ymin = ymn, ymax = ymx,
                                      start_time = st_time,
                                      end_time = en_time,
                                      outfile_name = "ERA5_LAND")
      reql[[1]]$variable <- c("2m_temperature",
                              "2m_dewpoint_temperature",
                              "surface_pressure",
                              "10m_u_component_of_wind",
                              "10m_v_component_of_wind",
                              "total_precipitation",
                              "surface_solar_radiation_downwards")
      suppressMessages({
        request_era5(request = reql, uid = uid, out_path = temp_dir, overwrite = TRUE)
      })
      
      req <- build_era5_request(xmin = xmn, xmax = xmx,
                                ymin = ymn, ymax = ymx,
                                start_time = st_time,
                                end_time = en_time,
                                outfile_name = "ERA5_ORI")
      req[[1]]$variable <- c("sea_surface_temperature",
                             "total_cloud_cover",
                             "mean_surface_net_long_wave_radiation_flux",
                             "mean_surface_downward_long_wave_radiation_flux",
                             "total_sky_direct_solar_radiation_at_surface")
      
      suppressMessages({
        request_era5(request = req, uid = uid, out_path = temp_dir, overwrite = TRUE)
      })
      ########################################################
      # Preprocess and combine ERA5 Land & ERA5
      
      r1 <- rast(paste0(temp_dir, "/ERA5_LAND_", year(st_time), "_", month(st_time), ".nc"))
      r2 <- rast(paste0(temp_dir, "/ERA5_ORI_", year(st_time), "_", month(st_time), ".nc"))
      sst <- r2[[(grepl("^sst", names(r2)))]]
      r2 <- r2[[!(grepl("^sst", names(r2)))]]
      
      r1 <- focal(r1, w=3, fun=mean, na.policy="only", na.rm=T)
      r2 <- focal(r2, w=3, fun=mean, na.policy="only", na.rm=T)
      
      r3 <- terra::project(r2, r1)
      sst <- terra::project(sst, r1)
      r1 <- c(r1, r3)
      
      terra::time(r1) <- as.numeric(str_split_i(names(r1), "_valid_time=", 2))
      names(r1) <- str_split_i(names(r1), "_valid_time=", 1)
      terra::time(sst) <- as.numeric(str_split_i(names(sst), "_valid_time=", 2))
      names(sst) <- str_split_i(names(sst), "_valid_time=", 1)
      varnames(sst) <- "sst"
      
      all_vars <- unique(names(r1))
      
      raster_list <- lapply(all_vars, function(var){
        layers_r1 <- r1[[names(r1) == var]]
        return(layers_r1)
      })
      
      # Create a SpatRasterDataset from the list of SpatRaster objects
      sds_out <- sds(raster_list)
      
      writeCDF(sds_out, paste0(area_dl_dir, "/ERA5_COMB_", year(st_time), "_", month(st_time), ".nc"))
      writeCDF(sst, paste0(area_dl_dir, "/ERA5_SST_", year(st_time), "_", month(st_time), ".nc"))
      
      unlink(paste0(temp_dir, "/ERA5_LAND_", year(st_time), "_", month(st_time), ".nc"))
      unlink(paste0(temp_dir, "/ERA5_ORI_", year(st_time), "_", month(st_time), ".nc"))
      
    } else {
      print("File exists. Skipping dates...")
    }
  }
}
