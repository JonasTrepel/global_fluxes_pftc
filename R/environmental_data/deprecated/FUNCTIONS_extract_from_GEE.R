Landsat_export_points <- function(pixel_coords_sf,
                               sample_id_from = "sample_id",
                               chunks_from = NULL,
                               max_chunk_size = 300,
                               drive_export_dir = "landsat_export",
                               file_prefix = "landsat_export",
                               start_doy = 0,
                               end_doy = 366,
                               start_date = "2015-01-01",
                               end_date = "2023-12-31",
                               scale = 30,
                               mask_value = 0,
                               cloud_max = 80,
                               geom_rmse_max = 15,
                               sun_elev_min = 0
){
  
  # confirm rgee is initialized
  tryCatch(rgee::ee_user_info(quiet = T),
           error = function(e) {
             stop("rgee not initialized!\nPlease intialize
                rgee. See: https://r-spatial.github.io/rgee/index.html")
           })
  
  # Check whether end_date was supplied and if not set to today's.
  if(end_date == "today") end_date <- as.character(Sys.Date())
  
  # Check whether pixel_coords_sf is an sf object with an sfc of points
  if(!("sfc_POINT" %in% class(sf::st_geometry(pixel_coords_sf)))) {
    stop("Invalid argument supplied for pixel_coords_sf!\n",
         "Please supply an object with 'sfc_POINT' geometries.")
  }
  
  # Check wether number of point coordinates exceeds 100 000
  if(length(sf::st_geometry(pixel_coords_sf)) > 100000){
    cat(crayon::red("Warning: Extraction requested for more than 100000 point locations!\n"))
    warning("Extraction requested for more than 100000 point locations!")
    answer <- readline("Would you like to continue nonetheless? (not recommended!) [y/n]: ")
    if(answer == "n"){
      cat("Okay, stopping extraction.\n")
      return(NULL)
    } else if(answer == "y"){
      cat("Okay, continuing...\n")
    } else {
      cat("Invalid answer, stopping extraction.\n")
      return(NULL)
    }
  }
  
  # Check whether columns exists if column selectors were supplied
  if((sample_id_from != "sample_id") & !(sample_id_from %in% names(pixel_coords_sf))) {
    stop("Invalid columns specificed for 'sample_id_from': ", sample_id_from)
  }
  
  # Prep Landsat Time series
  bands <- list("SR_B1",
                "SR_B2",
                "SR_B3",
                "SR_B4",
                "SR_B5",
                "SR_B6",
                "SR_B7",
                "QA_PIXEL",
                "QA_RADSAT")
  BAND_LIST <- rgee::ee$List(bands)
  
  # addon assets and bands
  ADDON <- rgee::ee$Image('JRC/GSW1_0/GlobalSurfaceWater')$
    float()$unmask(mask_value)
  ADDON_BANDLIST <- rgee::ee$List(list("max_extent"));
  
  # Blank image for "SR_B6" to replace in collections earlier than LS8
  ZERO_IMAGE <- rgee::ee$Image(0)$select(list("constant"),
                                         list("SR_B6"))$selfMask()
  
  # Set image properties to export
  PROPERTIES <- list("CLOUD_COVER",
                     "COLLECTION_NUMBER",
                     "DATE_ACQUIRED",
                     "GEOMETRIC_RMSE_MODEL",
                     "LANDSAT_PRODUCT_ID",
                     'LANDSAT_SCENE_ID',
                     "PROCESSING_LEVEL",
                     "SPACECRAFT_ID",
                     "SUN_ELEVATION")
  
  # Landsat Surface Reflectance collections
  ls5_1 <- rgee::ee$ImageCollection("LANDSAT/LT05/C02/T1_L2")
  ls5_2 <- rgee::ee$ImageCollection("LANDSAT/LT05/C02/T2_L2")
  ls7_1 <- rgee::ee$ImageCollection("LANDSAT/LE07/C02/T1_L2")
  ls7_2 <- rgee::ee$ImageCollection("LANDSAT/LE07/C02/T2_L2")
  ls8_1 <- rgee::ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")
  ls8_2 <- rgee::ee$ImageCollection("LANDSAT/LC08/C02/T2_L2")
  ls9_1 <- rgee::ee$ImageCollection("LANDSAT/LC09/C02/T1_L2")
  ls9_2 <- rgee::ee$ImageCollection("LANDSAT/LC09/C02/T2_L2")
  
  ALL_BANDS <- BAND_LIST$cat(ADDON_BANDLIST)
  
  # merge all collections into one
  LS_COLL <- ls5_1$
    merge(ls7_1$
            merge(ls8_1$
                    merge(ls9_1$
                            merge(ls5_2$
                                    merge(ls7_2$
                                            merge(ls8_2$
                                                    merge(ls9_2)))))))$
    filterDate(start_date, end_date)$
    filter(rgee::ee$Filter$calendarRange(start_doy,
                                         end_doy,
                                         "day_of_year"))$
    map(function(image){
      image = rgee::ee$Algorithms$If(image$bandNames()$size()$
                                       eq(rgee::ee$Number(10)),
                                     image,
                                     image$addBands(ZERO_IMAGE))
      return(image)})$
    map(function(image) {return(image$addBands(ADDON, ADDON_BANDLIST))})$
    select(ALL_BANDS)$
    map(function(image){ return(image$float())} )
  
  # FILter imagery by metadata
  LS_COLL <- LS_COLL$filterMetadata("CLOUD_COVER", "less_than", cloud_max)
  LS_COLL <- LS_COLL$filterMetadata("GEOMETRIC_RMSE_MODEL", "less_than", geom_rmse_max)
  LS_COLL <- LS_COLL$filterMetadata("SUN_ELEVATION", "greater_than", sun_elev_min)
  
  # Check if chunks_from was specified, if not determine chunks
  if(!is.null(chunks_from)){
    if(!(chunks_from %in% colnames(pixel_coords_sf))) {
      stop("Invalid colum name specified for chunks_from: ", chunks_from)
    }
  } else {
    n_chunks <- floor(nrow(pixel_coords_sf) / max_chunk_size) + 1
    pixel_coords_sf$chunk_id <-
      paste0("chunk_",
             sort(rep(1:n_chunks, max_chunk_size)))[1:nrow(pixel_coords_sf)]
    chunks_from <- "chunk_id"
  }
  
  # Status:
  cat(paste0("Exporting time-series for ",
             nrow(pixel_coords_sf),
             " pixels",
             " in ",
             length(unique(sf::st_drop_geometry(pixel_coords_sf) %>%
                             as.data.frame() %>%
                             .[,chunks_from])),
             " chunks.\n"))
  
  # Retrieve time-series by chunk
  task_list <- pixel_coords_sf %>%
    split(., sf::st_drop_geometry(.)[,chunks_from]) %>%
    purrr::map(function(chunk){
      # Status
      cat(paste0("Submitting task to EE for chunk_id: ",
                 sf::st_drop_geometry(chunk) %>%
                   as.data.frame() %>%
                   .[1, chunks_from],
                 ".\n"))
      # Upload chunk to sf to reduce size keep only necessary columns
      ee_chunk <- rgee::sf_as_ee(chunk[,c("geometry",
                                          sample_id_from,
                                          chunks_from)])
      # Retrieve Landsat time-series
      ee_chunk_export <- ee_chunk$map(function(feature){
        return(
          # Create FC containing a single empty image
          # This will ensure all bands are present in the export
          rgee::ee$ImageCollection$fromImages(
            list(rgee::ee$Image(list(0,0,0,0,0,0,0,0,0,0))$
                   select(list(0,1,2,3,4,5,6,7,8,9), ALL_BANDS)$
                   copyProperties(ls8_1$first())))$
            # Merge with extraction of time-series form Landsat collection
            merge(LS_COLL$filterBounds(feature$geometry()))$
            # For each image in the collection ....
            map(function(image){
              # Create a feature
              return(rgee::ee$Feature(feature$geometry(),
                                      # fill it with the point value extracted with
                                      # reduceRegion with first() reducer at scale
                                      image$reduceRegion(rgee::ee$Reducer$first(),
                                                         feature$geometry(),
                                                         scale))$
                       # copy the image properties to the feature
                       # (incl. date and image metadata) as specified above
                       copyProperties(image, PROPERTIES)$
                       # assign a pixel and chunk id columns for identification
                       set(sample_id_from, feature$get(sample_id_from))$
                       set(chunks_from, feature$get(chunks_from)))
            }))
      })$flatten()
      # Prepare export task
      chunk_task <- rgee::ee_table_to_drive(
        collection = ee_chunk_export,
        description = paste0("lsatTS_export_",
                             sf::st_drop_geometry(chunk) %>%
                               as.data.frame() %>%
                               .[1, chunks_from]),
        folder = drive_export_dir,
        fileNamePrefix = paste0(file_prefix,
                                "_",
                                sf::st_drop_geometry(chunk) %>%
                                  as.data.frame() %>%
                                  .[1, chunks_from]),
        timePrefix = F,
        fileFormat = "csv")
      # Submit export task
      chunk_task$start()
      
      # Return nothing
      return(chunk_task)
    })
  # Status update
  cat(crayon::green("Done!\n"))
  cat("You can monitor the progress of the task(s)",
      "using rgee's ee_monitoring() or the GEE WebAPI.\n")
  return(task_list)
}

Sentinel_export_points <- function(pixel_coords_sf,
                                  sample_id_from = "sample_id",
                                  chunks_from = NULL,
                                  max_chunk_size = 300,
                                  drive_export_dir = "sentinel_export",
                                  file_prefix = "sentinel_export",
                                  start_doy = 0,
                                  end_doy = 366,
                                  start_date = "2015-01-01",
                                  end_date = "2023-12-31",
                                  scale = 30,
                                  mask_value = 0,
                                  cloud_max = 80,
                                  sun_elev_min = 0
){
  
  # confirm rgee is initialized
  tryCatch(rgee::ee_user_info(quiet = T),
           error = function(e) {
             stop("rgee not initialized!\nPlease intialize
                rgee. See: https://r-spatial.github.io/rgee/index.html")
           })
  
  # Check whether end_date was supplied and if not set to today's.
  if(end_date == "today") end_date <- as.character(Sys.Date())
  
  # Check whether pixel_coords_sf is an sf object with an sfc of points
  if(!("sfc_POINT" %in% class(sf::st_geometry(pixel_coords_sf)))) {
    stop("Invalid argument supplied for pixel_coords_sf!\n",
         "Please supply an object with 'sfc_POINT' geometries.")
  }
  
  # Check wether number of point coordinates exceeds 100 000
  if(length(sf::st_geometry(pixel_coords_sf)) > 100000){
    cat(crayon::red("Warning: Extraction requested for more than 100000 point locations!\n"))
    warning("Extraction requested for more than 100000 point locations!")
    answer <- readline("Would you like to continue nonetheless? (not recommended!) [y/n]: ")
    if(answer == "n"){
      cat("Okay, stopping extraction.\n")
      return(NULL)
    } else if(answer == "y"){
      cat("Okay, continuing...\n")
    } else {
      cat("Invalid answer, stopping extraction.\n")
      return(NULL)
    }
  }
  
  # Check whether columns exists if column selectors were supplied
  if((sample_id_from != "sample_id") & !(sample_id_from %in% names(pixel_coords_sf))) {
    stop("Invalid columns specificed for 'sample_id_from': ", sample_id_from)
  }
  
  # Prep Landsat Time series
  bands <- list("B1",
                "B2",
                "B3",
                "B4",
                "B5",
                "B6",
                "B7",
                "B8",
                "B8A",
                "B9",
                "B11",
                "B12",
                "SCL",
                "MSK_CLDPRB",
                "MSK_SNWPRB")
  BAND_LIST <- rgee::ee$List(bands)
  
  # addon assets and bands
  ADDON <- rgee::ee$Image('JRC/GSW1_0/GlobalSurfaceWater')$
    float()$unmask(mask_value)
  ADDON_BANDLIST <- rgee::ee$List(list("max_extent"))
  
  # Set image properties to export
  PROPERTIES <- list("CLOUDY_PIXEL_PERCENTAGE",
                     "DATASTRIP_ID",
                     "DATATAKE_IDENTIFIER",
                     "DATATAKE_TYPE",
                     "GENERAL_QUALITY",
                     "HIGH_PROBA_CLOUDS_PERCENTAGE",
                     "MEAN_SOLAR_ZENITH_ANGLE",
                     "PRODUCT_ID",
                     "SPACECRAFT_NAME",
                     "GENERATION_TIME",
                     "GEOMETRIC_QUALITY")
  
  # Landsat Surface Reflectance collections
  s2 <- rgee::ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  
  ALL_BANDS <- BAND_LIST$cat(ADDON_BANDLIST)
  
  # merge all collections into one
  LS_COLL <- s2$
    filterDate(start_date, end_date)$
    filter(rgee::ee$Filter$calendarRange(start_doy,
                                         end_doy,
                                         "day_of_year"))$
    map(function(image) {return(image$addBands(ADDON, ADDON_BANDLIST))})$
    select(ALL_BANDS)$
    map(function(image){ return(image$float())} )
  
  # FILter imagery by metadata
  LS_COLL <- LS_COLL$filterMetadata("CLOUDY_PIXEL_PERCENTAGE", "less_than", cloud_max)
  # LS_COLL <- LS_COLL$filterMetadata("GEOMETRIC_RMSE_MODEL", "less_than", geom_rmse_max)
  LS_COLL <- LS_COLL$filterMetadata("MEAN_SOLAR_ZENITH_ANGLE", "less_than", (90 - sun_elev_min))
  
  # Check if chunks_from was specified, if not determine chunks
  if(!is.null(chunks_from)){
    if(!(chunks_from %in% colnames(pixel_coords_sf))) {
      stop("Invalid colum name specified for chunks_from: ", chunks_from)
    }
  } else {
    n_chunks <- floor(nrow(pixel_coords_sf) / max_chunk_size) + 1
    pixel_coords_sf$chunk_id <-
      paste0("chunk_",
             sort(rep(1:n_chunks, max_chunk_size)))[1:nrow(pixel_coords_sf)]
    chunks_from <- "chunk_id"
  }
  
  # Status:
  cat(paste0("Exporting time-series for ",
             nrow(pixel_coords_sf),
             " pixels",
             " in ",
             length(unique(sf::st_drop_geometry(pixel_coords_sf) %>%
                             as.data.frame() %>%
                             .[,chunks_from])),
             " chunks.\n"))
  
  # Retrieve time-series by chunk
  task_list <- pixel_coords_sf %>%
    split(., sf::st_drop_geometry(.)[,chunks_from]) %>%
    purrr::map(function(chunk){
      # Status
      cat(paste0("Submitting task to EE for chunk_id: ",
                 sf::st_drop_geometry(chunk) %>%
                   as.data.frame() %>%
                   .[1, chunks_from],
                 ".\n"))
      # Upload chunk to sf to reduce size keep only necessary columns
      ee_chunk <- rgee::sf_as_ee(chunk[,c("geometry",
                                          sample_id_from,
                                          chunks_from)])
      # Retrieve Landsat time-series
      ee_chunk_export <- ee_chunk$map(function(feature){
        return(
          # Create FC containing a single empty image
          # This will ensure all bands are present in the export
          LS_COLL$filterBounds(feature$geometry())$
            # For each image in the collection ....
            map(function(image){
              # Create a feature
              return(rgee::ee$Feature(feature$geometry(),
                                      # fill it with the point value extracted with
                                      # reduceRegion with first() reducer at scale
                                      image$reduceRegion(rgee::ee$Reducer$first(),
                                                         feature$geometry(),
                                                         scale))$
                       # copy the image properties to the feature
                       # (incl. date and image metadata) as specified above
                       copyProperties(image, PROPERTIES)$
                       # assign a pixel and chunk id columns for identification
                       set(sample_id_from, feature$get(sample_id_from))$
                       set(chunks_from, feature$get(chunks_from)))
            }))
      })$flatten()
      # Prepare export task
      chunk_task <- rgee::ee_table_to_drive(
        collection = ee_chunk_export,
        description = paste0("ssatTS_export_",
                             sf::st_drop_geometry(chunk) %>%
                               as.data.frame() %>%
                               .[1, chunks_from]),
        folder = drive_export_dir,
        fileNamePrefix = paste0(file_prefix,
                                "_",
                                sf::st_drop_geometry(chunk) %>%
                                  as.data.frame() %>%
                                  .[1, chunks_from]),
        timePrefix = F,
        fileFormat = "csv")
      # Submit export task
      chunk_task$start()
      
      # Return nothing
      return(chunk_task)
    })
  # Status update
  cat(crayon::green("Done!\n"))
  cat("You can monitor the progress of the task(s)",
      "using rgee's ee_monitoring() or the GEE WebAPI.\n")
  return(task_list)
}
