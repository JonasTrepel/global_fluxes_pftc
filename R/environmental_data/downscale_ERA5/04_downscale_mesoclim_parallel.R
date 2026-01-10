# library(terra)
# library(tidync, lib.loc = "/projappl/project_2003061/Rpackages")
# library(microclima, lib.loc = "/projappl/project_2003061/Rpackages")
# library(microclimf, lib.loc = "/projappl/project_2003061/Rpackages")
# library(micropoint, lib.loc = "/projappl/project_2003061/Rpackages")
# library(NicheMapR, lib.loc = "/projappl/project_2003061/Rpackages")
# library(mesoclim, lib.loc = "/projappl/project_2003061/Rpackages")
# library(microclimdata, lib.loc = "/projappl/project_2003061/Rpackages")
# library(tidyverse)
# library(ncdf4)
# library(curl)
# library(keyring)
# library(abind)

library(terra)
library(mesoclim)
library(tidyverse)
library(sf)

source("scr/functions/extract_clima_own.R")
source("scr/functions/extract_clima_own2.R")

p <- st_read("data/spatial_data/area_polygons.gpkg")

p %>% st_area()/1e+06

dl_dir <- "/scratch/project_2003061/ERA5/PFTC"

for(area_id in p$area){
  # area_id <- "Svalbard"
  print(area_id)
  
  pt <- p %>% filter(area == area_id)
  
  # if(as.numeric(st_area(pt)/1e6) > 400){
  #   
  #   pt <- st_make_grid(pt, n = 2) %>% 
  #     st_buffer(500, endCapStyle = "SQUARE") %>% 
  #     st_as_sf() %>% 
  #     mutate(area = area_id)
  #   
  # }
  pt <- pt %>% 
    mutate(segment = row_number())
  
  area_dl_dir <- paste0(dl_dir, "/", area_id, "/")
  
  for(sid in seq_along(pt$segment)){
    # sid <- 1
    pt2 <- pt %>% filter(segment == sid)
    
    dem <- rast(paste0(area_dl_dir, "/ERA5-Land_DEM.tif"))
    lsm <- rast(paste0(area_dl_dir, "/ERA5-Land_water.tif"))
    
    fdem <- rast(paste0(area_dl_dir, "/ALOSDEM.tif")) %>% 
      terra::aggregate(2, na.rm = TRUE)
    
    fdem <- fdem %>% 
      crop(st_transform(pt2, crs = st_crs(fdem)), snap = "out")
    
    dem <- crop(dem, project(fdem, crs(dem)), snap = "out")
    lsm <- crop(lsm, project(fdem, crs(lsm)), snap = "out")
    
    f <- list.files(area_dl_dir, pattern = "ERA5_COMB", full.names = TRUE)
    
    for(im in f){
      # im <- f[[1]]
      if(!file.exists(gsub("ERA5_COMB_","mesoclimate_",gsub(".nc$","_windspeed.tif",im)))){
        sst <- rast(gsub("ERA5_COMB_","ERA5_SST_",im)) - 273.15
        sst <- crop(sst, dem, snap = "out")
        
        era5input <- era5toclimarray_own(ncfile = im, toArrays = FALSE,
                                         dtmc = dem, lsm = lsm, aoi = dem, dtr_cor_fac = 1.285)
        terra::time(sst) <- terra::time(era5input$temp)
        
        # fdem <- terra::aggregate(fdem, 4, na.rm = TRUE)
        dtmm <- terra::aggregate(fdem, 4, na.rm = TRUE)
        # plot(fdem)
        system.time({
          mesoclimate <- spatialdownscale_own_limited(climdata=era5input, sst=sst, dtmf = fdem, dtmm = dtmm, 
                                              basins = NA,  wca=NA, skyview=NA, horizon=FALSE, cad=FALSE, 
                                              coastal = TRUE, thgto =2, whgto=2, include_tmean=TRUE, 
                                              rhmin = 20, pksealevel = FALSE, patchsim = FALSE,  
                                              terrainshade = FALSE, precipmethod = "Tps", fast = TRUE, 
                                              noraincut = 0.01, toArrays=FALSE)
        })
        varnames(mesoclimate$temp) <- "T2m"
        varnames(mesoclimate$relhum) <- "relhum"
        varnames(mesoclimate$windspeed) <- "windspeed"
        mesoclimate$temp <- round(mesoclimate$temp*10)
        mesoclimate$relhum <- round(mesoclimate$relhum*10)
        mesoclimate$windspeed <- round(mesoclimate$windspeed*100)
        writeRaster(mesoclimate$temp, gsub("ERA5_COMB_","mesoclimate_",gsub(".nc$","_T2m.tif",im)), overwrite = TRUE, datatype = "INT2S")
        writeRaster(mesoclimate$relhum, gsub("ERA5_COMB_","mesoclimate_",gsub(".nc$","_relhum.tif",im)), overwrite = TRUE, datatype = "INT2S")
        writeRaster(mesoclimate$windspeed, gsub("ERA5_COMB_","mesoclimate_",gsub(".nc$","_windspeed.tif",im)), overwrite = TRUE, datatype = "INT2S")
        
      }
    }
  }
}
