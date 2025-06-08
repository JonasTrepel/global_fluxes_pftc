# Needed libraries
library(terra)
library(rgee)
library(sf)
library(tidyverse)

# ee_install()
ee_Initialize(drive = T)
# reticulate::py_available()
# ee_check()
ee$Authenticate(auth_mode='notebook')
ee$Initialize(project='radiant-planet-412111')
ee$String('Hello from the Earth Engine servers!')$getInfo()

# Direction where imagery will be stores (in area-specific sub-directories)
base_dem_dir <- "C://MyTemp/DEMS"
if(!dir.exists(base_dem_dir)){
  dir.create(base_dem_dir)
}

utmall <- st_read("data/utm_zones.gpkg") %>% 
  filter(ZONE != 0)

# Read in the study area points
aoi <- st_read("data/coordinates.gpkg") %>%
  st_transform(crs = 4326) %>% 
  mutate(name = as.character(area))

for(i in unique(aoi$name)){
  # i <- aoi$name[1]
  # WGS84 UTM zones to set the correct projection
  
  aoii <- aoi %>% filter(name == i) %>% 
    st_buffer(1000) %>% 
    st_bbox() %>% 
    st_as_sfc() %>% 
    st_as_sf() %>% 
    mutate(name = i)
  
  if(!file.exists(paste0(base_dem_dir, "/",aoii$name,"/ESALC.tif"))){
    print(i)
    
    if(!dir.exists(paste0(base_dem_dir, "/",aoii$name))){
      dir.create(paste0(base_dem_dir, "/",aoii$name))
    }
    
    
    # WGS84 UTM zones to set the correct projection
    utm <- utmall[aoii,] # Which zone the study points falls in
    
    if(nrow(utm) == 0){
      utm <- utmall[st_nearest_feature(aoii, utmall),]
    }
    
    lat <- st_coordinates(aoii %>% st_centroid())[,"Y"]
    lon <- st_coordinates(aoii %>% st_centroid())[,"X"]
    utm$ZONE <- ifelse(nchar(utm$ZONE) == 1, paste0("0",utm$ZONE), utm$ZONE)
    epsg <- as.numeric(ifelse(lat > 0, paste0(326, utm$ZONE), paste0(327, utm$ZONE)))
    
    # From point to polygon (20km x 20km)
    # aoii <- aoii %>% st_transform(crs = epsg) %>% 
    #   st_buffer(1000) %>% st_bbox() %>% 
    #   st_as_sfc() %>% st_as_sf() %>% 
    #   mutate(name = aoii$name)
    
    # AOI polygon to GEE format
    aoi_ee <- aoii %>% st_transform(crs = 4326) %>% 
      st_geometry() %>% 
      sf_as_ee()
    
    # ALOS DEM
    dataset <- ee$ImageCollection('JAXA/ALOS/AW3D30/V3_2')$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("DSM")
    
    e <- try({ei <- ee_print(dataset, quiet = T)}, silent = T)
    
    if(class(e) != "try-error"){
      dataset <- dataset$max()
      dataset <- dataset$reproject(crs = paste0("EPSG:",epsg), scale = 30)
      # ee_print(dataset)
      task_img <- ee_image_to_drive(
        image = dataset,
        fileFormat = "GEO_TIFF",folder = paste0("dems_",i),
        region = aoi_ee,
        scale = 30,
        fileNamePrefix = paste0(aoii$name,"_ALOSDEM")
      )
      
      task_img$start()
      ee_monitoring(task_img, max_attempts = 150)
      
      ee_drive_to_local(task = task_img, 
                        dsn = paste0(base_dem_dir, "/",aoii$name,"/ALOSDEM"))
      
    }
    
    # ESA LAND COVER
    dataset <- ee$ImageCollection('ESA/WorldCover/v100')$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("Map")
    
    e <- try({ei <- ee_print(dataset, quiet = T)}, silent = T)
    
    if(class(e) != "try-error"){
      dataset <- dataset$max()
      dataset <- dataset$reproject(crs = paste0("EPSG:",epsg), scale = 10)
      # ee_print(dataset)
      task_img <- ee_image_to_drive(
        image = dataset,
        fileFormat = "GEO_TIFF",folder = paste0("dems_",i),
        region = aoi_ee,
        scale = 10,
        fileNamePrefix = paste0(aoii$name,"_ESALC")
      )
      
      task_img$start()
      ee_monitoring(task_img, max_attempts = 150)
      
      ee_drive_to_local(task = task_img, 
                        dsn = paste0(base_dem_dir, "/",aoii$name,"/ESALC"))
      
    }
    
    ee_clean_container(name = paste0("dems_",i), type = "drive", quiet = FALSE)
    
  }
}

getmode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

all <- lapply(unique(aoi$name), function(i){
  print(i)
  # i <- "China"
  # ESA WorldCover
  if(file.exists(paste0(base_dem_dir,"/",i,"/ESALC.tif"))){
    r <- rast(paste0(base_dem_dir,"/",i,"/ESALC.tif"))
    r[r %in% c(70,80)] <- NA
    
    pt <- aoi %>% 
      filter(name == i) %>% 
      st_transform(crs(r, proj = T)) %>% 
      st_buffer(10)
    
    res <- terra::extract(r, pt, fun = getmode)
    names(res)[2] <- "esa_worldcover_10"
    
    dt <- pt %>% 
      st_drop_geometry() %>% 
      mutate(esa_worldcover_10 = res$esa_worldcover_10)
    
    pt <- aoi %>% 
      filter(name == i) %>% 
      st_transform(crs(r, proj = T)) %>% 
      st_buffer(50)
    
    res <- terra::extract(r, pt, fun = getmode)
    names(res)[2] <- "esa_worldcover_50"
    
    dt <- dt %>% 
      mutate(esa_worldcover_50 = res$esa_worldcover_50)
    
  } else {
    pt <- aoi %>% 
      filter(name == i) %>% 
      st_drop_geometry() %>% 
      mutate(esa_worldcover_50 = NA,
             esa_worldcover_500 = NA)
  }
  
  # ELevation
  if(file.exists(paste0(base_dem_dir,"/",i,"/ALOSDEM.tif"))){
    r <- rast(paste0(base_dem_dir,"/",i,"/ALOSDEM.tif"))
    
    pt <- aoi %>% 
      filter(name == i) %>% 
      st_transform(crs(r, proj = T)) %>% 
      st_buffer(10)
    
    res <- terra::extract(r, pt, fun = mean, na.rm = T)
    names(res)[2] <- "elevation"
    
    dt <- dt %>% 
      mutate(elevation = round(res$elevation))
    
  } else {
    dt <- dt %>% 
      mutate(elevation = NA)
  }
  
  return(dt)
})

all <- bind_rows(all) %>% 
  select(-name)

write_csv(all, "data/environmental_data/worldcover.csv")
