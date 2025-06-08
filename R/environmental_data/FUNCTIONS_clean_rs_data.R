dn_to_ref <- function(r){
  r <- r*2.75e-05-0.2 # Scale and offset from: https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2
  return(r)
}

# Function NDVI
ndvi_fun <- function(red, nir){
  ndvi <- (nir-red)/(nir+red)
  return(ndvi)
}

# Function to extract cloud pixels from QA layer
mask_snow <- function(x) {as.numeric(intToBits(x)[6])}
mask_clear <- function(x) {as.numeric(intToBits(x)[7])}

clean_landsat <- function(ld, sample_id_from){
  
  # parse coords
  coords <- str_extract(string = ld$.geo, pattern = "(?<=\\[).*(?=\\])")
  coords <- matrix(unlist(strsplit(coords, ',')), ncol = 2, byrow = T)
  ld$Y <- as.numeric(coords[,2])
  ld$X <- as.numeric(coords[,1])
  
  ld <- ld %>% 
    filter(max_extent == 0) %>% 
    select(all_of(sample_id_from), X, Y, DATE_ACQUIRED, SR_B1:SR_B7, QA_PIXEL, QA_RADSAT,
           SPACECRAFT_ID, LANDSAT_PRODUCT_ID, CLOUD_COVER, GEOMETRIC_RMSE_MODEL, SUN_ELEVATION,
           PROCESSING_LEVEL, COLLECTION_NUMBER) %>% 
  filter(QA_PIXEL > 0, (!is.na(QA_PIXEL)))
  
  # rename bands seperatly lsat 5/7 and lsat 8
  ld57 <- ld %>% 
    filter(SPACECRAFT_ID %in% c("LANDSAT_5", "LANDSAT_7")) %>% 
    rename(blue = SR_B1,
           green = SR_B2,
           red = SR_B3,
           nir = SR_B4,
           swir1 = SR_B5,
           swir2 = SR_B7) %>% 
    select(-SR_B6)
  
  ld89 <- ld %>% 
    filter(SPACECRAFT_ID %in% c("LANDSAT_8", "LANDSAT_9")) %>% 
    rename(ublue = SR_B1,
           blue = SR_B2,
           green = SR_B3,
           red = SR_B4,
           nir = SR_B5,
           swir1 = SR_B6,
           swir2 = SR_B7)
  
  ld <- bind_rows(ld89, ld57) %>% 
    rename(date = DATE_ACQUIRED,
           satellite = SPACECRAFT_ID)
  
  ld <- ld %>% 
    filter(unlist(lapply(QA_PIXEL, mask_clear)) == 1 |
             unlist(lapply(QA_PIXEL, mask_snow)) == 1) %>% 
    mutate(across(ublue:swir2, ~(.x*0.0000275-0.2))) %>% 
    filter(if_all(blue:nir, ~(.x > 0.005))) %>%
    filter(if_all(blue:nir, ~(.x < 1)))
  
  ld <- ld %>% 
    mutate(NDVI = ndvi_fun(red, nir)) %>% 
    filter(NDVI > (-1), NDVI < 1) %>% 
    relocate(NDVI, .after = swir2)
  
  return(ld)
  
}


clean_sentinel <- function(s2d, sample_id_from = "SITE_ID"){
  
  # parse coords
  coords <- str_extract(string = s2d$.geo, pattern = "(?<=\\[).*(?=\\])")
  coords <- matrix(unlist(strsplit(coords, ',')), ncol = 2, byrow = T)
  s2d$Y <- as.numeric(coords[,2])
  s2d$X <- as.numeric(coords[,1])
  
  s2d <- s2d %>% 
    filter(max_extent == 0) %>% 
    mutate(DATE_ACQUIRED = ymd(unlist(lapply(PRODUCT_ID, function(x) as.character(substr(str_split(x,"_")[[1]][3], 1, 8)) )))) %>% 
    select(all_of(sample_id_from), X, Y, DATE_ACQUIRED, B1:B9, SCL,
           SPACECRAFT_NAME, PRODUCT_ID, CLOUDY_PIXEL_PERCENTAGE, MEAN_SOLAR_ZENITH_ANGLE) %>% 
    filter(SCL %in% c(4:5,11), (!is.na(SCL)))
  
  s2d <- s2d %>% 
    rename(ublue = B1,
           blue = B2,
           green = B3,
           red = B4,
           rededge1 = B5,
           rededge2 = B6,
           rededge3 = B7,
           nir = B8,
           rededge4 = B8A,
           watervapor = B9,
           swir1 = B11,
           swir2 = B12) %>% 
    relocate(swir1:swir2, .after = watervapor)
    
  
  s2d <- s2d %>% 
    rename(date = DATE_ACQUIRED,
           satellite = SPACECRAFT_NAME)
  
  s2d <- s2d %>% 
    mutate(across(ublue:swir2, ~(.x*0.0001))) %>% 
    filter(if_all(blue:nir, ~(.x > 0.005))) %>% 
    filter(if_all(blue:nir, ~(.x < 1)))
  
  s2d <- s2d %>% 
    mutate(NDVI = ndvi_fun(red, nir)) %>% 
    filter(NDVI > (-1), NDVI < 1) %>% 
    relocate(NDVI, .after = swir2)
  
  return(s2d)
  
}
