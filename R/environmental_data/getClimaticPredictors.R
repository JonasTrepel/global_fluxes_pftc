# extract environmental Covs 

library(exactextractr)
library(tidyverse)
library(data.table)
library(sf)
library(terra)

dtCoords <- read_sf("data/processed_data/preliminary_data/prelim_flux_loc.gpkg")


#### EXTRACT CONTINUOUS COVS #### ----------------------------------

colNames <- c(
  #### Environmental ####
  "MMP", ## MAP
  "MAT", ## MAT
  "MaxTemp", ## MaxTemp
  "MinTemp", ##MinTemp
  "BurnedAreaMean"
)

covPaths <- c(
  #### Environmental ####
  "../GlobalChangeEcosystemFunctioningGitR/data/spatialData/climateData/currentMeanMonthlyPrec20092019.tif", ## MMP
  "../GlobalChangeEcosystemFunctioningGitR/data/spatialData/climateData/currentMeanTemp20092019.tif", ##MAT
  "../GlobalChangeEcosystemFunctioningGitR/data/spatialData/climateData/currentMaxTemp20092019.tif", ##Max Temp
  "../GlobalChangeEcosystemFunctioningGitR/data/spatialData/climateData/currentMinTemp20092019.tif", ## MinTemp
  "../GlobalChangeEcosystemFunctioningGitR/data/spatialData/otherCovariates/BurnedAreaMean20012023.tif"
)

covs <- data.table(
  colName = colNames, 
  covPath = covPaths
) %>% filter(!is.na(covPaths))


pasRawCovs <- dtCoords %>% as.data.table() %>% mutate(geom = NULL)

for(i in 1:nrow(covs)){
  
  covR <- rast(covs[i, ]$covPath)
  
  worldGridTrans <- st_transform(dtCoords, crs = st_crs(covR)) %>% st_buffer(100)
  
  extr <- exactextractr::exact_extract(covR, 
                                       worldGridTrans, 
                                       fun = "mean")
  extrDT <- data.table(
    extrCol = extr
  )
  setnames(extrDT, "extrCol", covs[i, ]$colName)
  
  pasRawCovs <- cbind(pasRawCovs, extrDT)
  
  print(paste0(i, "/", nrow(covs)))
  
}

pasCovsDT <- pasRawCovs %>% 
  as.data.table() %>% 
  mutate(x = NULL) %>% 
  dplyr::select(plot_id, tier, MMP, MAT, MaxTemp, MinTemp, BurnedAreaMean) %>% 
  unique()

fwrite(pasCovsDT, "data/processed_data/data_fragments/climaticPredictors.csv")
