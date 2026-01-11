spatialdownscale_own_limited <- function(climdata, sst, dtmf, dtmm = NA, basins = NA, wca = NA, 
                                 skyview = NA, horizon = NA, cad = TRUE, coastal = TRUE, thgto = 2, 
                                 whgto = 2, include_tmean = FALSE, rhmin = 20, pksealevel = TRUE, 
                                 patchsim = FALSE, terrainshade = TRUE, precipmethod = "Elev", 
                                 fast = TRUE, noraincut = 0, toArrays = FALSE) 
{
  input_class <- lapply(lapply(climdata, class), "[", 1)
  if (any(input_class == "PackedSpatRaster")){
    climdata[which(input_class == "PackedSpatRaster")] <- lapply(climdata[which(input_class == 
                                                                                  "PackedSpatRaster")], unwrap)
  }
  if (any(input_class == "array")){
    climdata[which(input_class == "array")] <- lapply(climdata[which(input_class == 
                                                                       "array")], .rast, tem = climdata$dtm)
  }
  if (inherits(sst, "PackedSpatRaster")){
    sst <- unwrap(sst)
  }
  
  tme <- as.POSIXlt(climdata$tme, tz = "UTC")
  tint <- as.numeric(tme[2]) - as.numeric(tme[1])
  hourly <- TRUE
  if (abs(tint - 86400) < 5){
    hourly = FALSE
  }
  mmry <- mem_info(dtmf, n = length(tme) * 24, print = FALSE)
  if (mmry["needed"] > (0.5 * mmry["available"]) & mmry["needed"] < 
      mmry["available"]){
    warning("High free memory use predicted - consider using spatialdownscale_tile!!!")
  }
  
  if (mmry["needed"] > mmry["available"]){
    warning("Memory demand predicted to exceed available memory - use spatialdownscale_tile option or reduce time period!!!")
  }
  
  if (hourly) {
    tc <- climdata$temp
  } else {
    tmin <- climdata$tmin
    tmax <- climdata$tmax
    tmean <- .hourtoday(temp_dailytohourly(tmin, tmax, climdata$tme), 
                        mean)
  }
  dtmc <- climdata$dtm
  rh <- climdata$relhum
  pk <- climdata$pres
  wspeed <- climdata$windspeed
  wdir <- climdata$winddir
  swrad <- climdata$swrad
  lwrad <- climdata$lwrad
  prec <- climdata$prec
  whgti <- climdata$windheight_m
  thgti <- climdata$tempheight_m
  message("Downscaling wind...")
  uzf <- winddownscale(wspeed, wdir, dtmf, dtmm, dtmc, wca, 
                       whgti, whgto)
  gc()
  message("Downscaling temperature...")
  tcf <- temphrly_downscale(climdata, sst, dtmf, dtmm, 
                            basins, uzf, cad, coastal, thgto, whgto)
  gc()
  message("Downscaling relative humidity")
  rhf <- relhumdownscale(rh, tc, tcf, dtmc, rhmin)
  gc()
  
  message("Formatting output...")
  if (hourly) {
    terra::time(tcf) <- climdata$tme
    names(tcf) <- climdata$tme
  } else {
    terra::time(tminf) <- climdata$tme
    terra::time(tmaxf) <- climdata$tme
    terra::time(tmeanf) <- climdata$tme
    names(tminf) <- climdata$tme
    names(tmaxf) <- climdata$tme
    names(tmeanf) <- climdata$tme
  }
  terra::time(rhf) <- climdata$tme
  terra::time(uzf) <- climdata$tme
  names(rhf) <- climdata$tme
  names(uzf) <- climdata$tme
  out <- list(tme = climdata$tme, windheight_m = whgto, 
              tempheight_m = thgto)
  if (hourly) {
    climout <- list(temp = tcf, relhum = rhf,
                    windspeed = uzf)
    if (toArrays) 
      climout <- lapply(climout, as.array)
    out <- c(out, climout)
  }
  else {
    if (!include_tmean) 
      climout <- list(tmin = tminf, tmax = tmaxf, relhum = rhf, 
                      windspeed = uzf)
    if (include_tmean) 
      climout <- list(tmean = tmeanf, tmin = tminf, tmax = tmaxf, 
                      relhum = rhf, windspeed = uzf)
    if (toArrays) 
      climout <- lapply(climout, as.array)
    out <- c(out, climout)
  }
  return(out)
}

spatialdownscale_own <- function(climdata, sst, dtmf, dtmm = NA, basins = NA, wca = NA, 
          skyview = NA, horizon = NA, cad = TRUE, coastal = TRUE, thgto = 2, 
          whgto = 2, include_tmean = FALSE, rhmin = 20, pksealevel = TRUE, 
          patchsim = FALSE, terrainshade = TRUE, precipmethod = "Elev", 
          fast = TRUE, noraincut = 0, toArrays = FALSE) 
{
  input_class <- lapply(lapply(climdata, class), "[", 1)
  if (any(input_class == "PackedSpatRaster")){
    climdata[which(input_class == "PackedSpatRaster")] <- lapply(climdata[which(input_class == 
                                                                                  "PackedSpatRaster")], unwrap)
  }
  if (any(input_class == "array")){
    climdata[which(input_class == "array")] <- lapply(climdata[which(input_class == 
                                                                       "array")], .rast, tem = climdata$dtm)
  }
  if (inherits(sst, "PackedSpatRaster")){
    sst <- unwrap(sst)
  }
  
  tme <- as.POSIXlt(climdata$tme, tz = "UTC")
  tint <- as.numeric(tme[2]) - as.numeric(tme[1])
  hourly <- TRUE
  if (abs(tint - 86400) < 5){
    hourly = FALSE
  }
  mmry <- mem_info(dtmf, n = length(tme) * 24, print = FALSE)
  if (mmry["needed"] > (0.5 * mmry["available"]) & mmry["needed"] < 
      mmry["available"]){
    warning("High free memory use predicted - consider using spatialdownscale_tile!!!")
  }
    
  if (mmry["needed"] > mmry["available"]){
    warning("Memory demand predicted to exceed available memory - use spatialdownscale_tile option or reduce time period!!!")
  }
  
  if (hourly) {
    tc <- climdata$temp
  } else {
    tmin <- climdata$tmin
    tmax <- climdata$tmax
    tmean <- .hourtoday(temp_dailytohourly(tmin, tmax, climdata$tme), 
                        mean)
  }
  dtmc <- climdata$dtm
  rh <- climdata$relhum
  pk <- climdata$pres
  wspeed <- climdata$windspeed
  wdir <- climdata$winddir
  swrad <- climdata$swrad
  lwrad <- climdata$lwrad
  prec <- climdata$prec
  whgti <- climdata$windheight_m
  thgti <- climdata$tempheight_m
  message("Downscaling wind...")
  uzf <- winddownscale(wspeed, wdir, dtmf, dtmm, dtmc, wca, 
                       whgti, whgto)
  uu <- wspeed * cos(wdir * pi/180)
  vv <- wspeed * sin(wdir * pi/180)
  if (crs(uu) != crs(dtmf)) 
    uu <- project(uu, crs(dtmf))
  if (crs(vv) != crs(dtmf)) 
    vv <- project(vv, crs(dtmf))
  uu <- .resample(uu, dtmf, method = "cubicspline")
  vv <- .resample(vv, dtmf, method = "cubicspline")
  wdf <- atan2(vv, uu) * 180/pi
  wdf <- .rast(.is(wdf)%%360, dtmf)
  wdf <- mask(wdf, dtmf)
  message("Downscaling temperature...")
  tcf <- temphrly_downscale(climdata, sst, dtmf, dtmm, 
                            basins, uzf, cad, coastal, thgto, whgto)
  
  message("Downscaling relative humidity")
  rhf <- relhumdownscale(rh, tc, tcf, dtmc, rhmin)
  
  message("Downscaling pressure...")
  pkf <- presdownscale(pk, dtmf, dtmc, sealevel = pksealevel)
  
  message("Downscaling SW radiation...")
  if (patchsim) {
    mws <- mean(as.vector(uzf), na.rm = TRUE)
    tstep <- as.numeric(tme[2]) - as.numeric(tme[1])
    dtr <- tstep * mws
    sze <- max(dim(dtmf)[1] * res(dtmf)[1], dim(dtmf)[2] * 
                 res(dtmf)[2])
    nsim <- round((dtr/sze) * dim(swrad)[3]) + 1
    if (nsim > dim(swrad)[3]){
      nsim <- dim(swrad)[3]
    }
  } else {
    nsim <- dim(swrad)[3]
  } 
  swf <- swdownscale(swrad, tme, dtmf, dtmc, patchsim, nsim, 
                     horizon, terrainshade)
  totswrad <- swf$swf
  if (terrainshade) {
    difrad <- swf$drf
  } else {
    difrad = NA
  }
  message("Downscaling LW radiation...")
  if (hourly){
    lwf <- lwdownscale(lwrad, tc, tcf, tme, dtmf, dtmc, skyview, 
                       terrainshade)
  } else {
    lwf <- lwdownscale(lwrad, tmean, tmeanf, tme, dtmf, 
                       dtmc, skyview, terrainshade)
  } 
  message("Downscaling precipitation...")
  precf <- precipdownscale(prec, dtmf, dtmc, precipmethod, 
                           fast = FALSE, noraincut, patchsim, nsim)
  plot(dtmf)
  message("Formatting output...")
  if (hourly) {
    terra::time(tcf) <- climdata$tme
    names(tcf) <- climdata$tme
  }
  else {
    terra::time(tminf) <- climdata$tme
    terra::time(tmaxf) <- climdata$tme
    terra::time(tmeanf) <- climdata$tme
    names(tminf) <- climdata$tme
    names(tmaxf) <- climdata$tme
    names(tmeanf) <- climdata$tme
  }
  terra::time(rhf) <- climdata$tme
  terra::time(pkf) <- climdata$tme
  terra::time(totswrad) <- climdata$tme
  terra::time(lwf) <- climdata$tme
  terra::time(uzf) <- climdata$tme
  terra::time(wdf) <- climdata$tme
  terra::time(precf) <- climdata$tme
  names(rhf) <- climdata$tme
  names(pkf) <- climdata$tme
  names(totswrad) <- climdata$tme
  names(lwf) <- climdata$tme
  names(uzf) <- climdata$tme
  names(wdf) <- climdata$tme
  names(precf) <- climdata$tme
  out <- list(dtm = dtmf, tme = climdata$tme, windheight_m = whgto, 
              tempheight_m = thgto)
  if (hourly) {
    climout <- list(temp = tcf, relhum = rhf, pres = pkf, 
                    swrad = totswrad, lwrad = lwf, windspeed = uzf, winddir = wdf, 
                    prec = precf)
    if (toArrays) 
      climout <- lapply(climout, as.array)
    out <- c(out, climout)
  }
  else {
    if (!include_tmean) 
      climout <- list(tmin = tminf, tmax = tmaxf, relhum = rhf, 
                      pres = pkf, swrad = totswrad, lwrad = lwf, windspeed = uzf, 
                      winddir = wdf, prec = precf)
    if (include_tmean) 
      climout <- list(tmean = tmeanf, tmin = tminf, tmax = tmaxf, 
                      relhum = rhf, pres = pkf, swrad = totswrad, lwrad = lwf, 
                      windspeed = uzf, winddir = wdf, prec = precf)
    if (toArrays) 
      climout <- lapply(climout, as.array)
    out <- c(out, climout)
  }
  if (terrainshade) {
    terra::time(difrad) <- climdata$tme
    names(difrad) <- climdata$tme
    out$difrad <- difrad
  }
  return(out)
}


era5toclimarray_own <- function (ncfile, dtmc, lsm, aoi = NA, dtr_cor_fac = 1.285, toArrays = TRUE, 
          zo = 10) 
{
  if (class(aoi)[1] != "logical") {
    if (!class(aoi)[1] %in% c("SpatRaster", "SpatVector", 
                              "sf")) 
      stop("Parameter aoi NOT of suitable spatial class ")
    if (class(aoi)[1] == "sf") 
      aoi <- vect(aoi)
    if (ext(dtmc) < ext(project(aoi, crs(dtmc)))) 
      stop("dtmc smaller than aoi!!!")
    if (ext(lsm) < ext(project(aoi, crs(dtmc)))) 
      stop("dtmc smaller than aoi!!!")
  }
  if (class(aoi)[1] == "logical"){
    aoi <- dtmc
  }
  units(dtmc) <- "m"
  names(dtmc) <- "Elevation"
  
  era5vars <- c("t2m","d2m","sp",
                "u10","v10","tp",
                "ssrd","fdir","avg_snlwrf",
                "avg_sdlwrf","tcc")
  
  varname_list <- c("t2m","d2m","sp",
                    "u10","v10","tp",
                    "ssrd","fdir","avg_snlwrf",
                    "avg_sdlwrf","tcc")
  
  nc_dat = ncdf4::nc_open(ncfile)
  timedim <- extract_timedim(nc_dat)
  timedim$units <- "seconds since 1970-01-01"
  base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), 
                              tz = "UTC")
  nc_datetimes <- c(timedim$vals)
  nc_datetimes <- nc_datetimes * ifelse(grepl("hours", timedim$units), 
                                        3600, 1)
  
  tme <- as.POSIXlt(nc_datetimes, tz = "UTC")
  
  # v <- "sp"
  var_list <- lapply(varname_list, function(v) {
    v2 <- paste0("var_", which(varname_list == v))
    
    r <- terra::rast(ncfile, subds = v2)
    r <- terra::crop(r, aoi)
    varnames(r) <- v
    names(r) <- rep(v, times = nlyr(r))
    
    return(r)
  })
  
  names(var_list) <- varname_list
  
  t2m <- var_list$t2m
  d2m <- var_list$d2m
  sp <- var_list$sp
  tp <- var_list$tp
  u10 <- var_list$u10
  v10 <- var_list$v10
  tcc <- var_list$tcc
  msnlwrf <- var_list$avg_snlwrf
  msdwlwrf <- var_list$avg_sdlwrf
  fdir <- var_list$fdir/3600
  ssrd <- var_list$ssrd/3600
  pres <- var_list$sp / 1000
  
  t2m <- t2m - 273.15
  lsm_e <- crop(lsm, dtmc)
  tmn <- .ehr(.hourtoday(as.array(t2m), min))
  a <- as.array(rep(lsm_e, dim(tmn)[3]))
  mu <- (1 - a) * dtr_cor_fac + 1
  tc <- .rast(((as.array(t2m)) - tmn) * mu + tmn, t2m)
  tc <- terra::project(tc, crs(aoi))
  d2m <- terra::project(d2m, crs(aoi))
  pres <- terra::project(pres, crs(aoi))
  u10 <- terra::project(u10, crs(aoi))
  v10 <- terra::project(v10, crs(aoi))
  tp <- terra::project(tp, crs(aoi))
  msdwlwrf <- terra::project(msdwlwrf, crs(aoi))
  fdir <- terra::project(fdir, crs(aoi))
  ssrd <- terra::project(ssrd, crs(aoi))
  dtmc <- project(dtmc, crs(aoi))
  lsm <- project(lsm, crs(aoi))
  ea <- .rast(.satvap(as.array(d2m) - 273.15), tc)
  temp <- as.array(tc)
  pres <- as.array(pres)
  relhum <- (as.array(ea)/.satvap(temp)) * 100
  swrad <- as.array(ssrd)
  difrad <- swrad - as.array(fdir)
  lwrad <- as.array(msdwlwrf)
  windspeed <- sqrt(as.array(u10)^2 + as.array(v10)^2) * log(67.8 * 
                                                               zo - 5.42)/log(67.8 * 10 - 5.42)
  winddir <- as.array((terra::atan2(u10, v10) * 180/pi + 180)%%360)
  prec <- as.array(tp) * 1000
  out <- list(dtm = dtmc, tme = tme, windheight_m = zo, tempheight_m = 2)
  climout <- list(temp = temp, relhum = relhum, pres = pres, 
                  swrad = swrad, difrad = difrad, lwrad = lwrad, windspeed = windspeed, 
                  winddir = winddir, prec = prec)
  climunits <- c("degC", "%", "kPa", "watt/m^2", "watt/m^2", 
                 "watt/m^2", "m/s", "deg", "mm")
  climround <- c(2, 1, 1, 1, 1, 1, 3, 1, 4)
  for (n in 1:length(climout)) {
    v <- names(climout)[n]
    climout[[v]] <- round(climout[[v]], climround[n])
  }
  if (toArrays == FALSE) {
    climout <- lapply(climout, .rast, tem = dtmc)
    for (n in 1:length(climout)) {
      v <- names(climout)[n]
      u <- climunits[n]
      terra::time(climout[[v]]) <- tme
      names(climout[[v]]) <- tme
      units(climout[[v]]) <- u
    }
  }
  out <- c(out, climout)
  return(out)
}

# nc <- "/scratch/project_2003061/temp/ERA5_COMB_2024_1.nc"
# lsm_path <- "/scratch/project_2003061/ERA5/ERA5_LAND_water.nc"
extract_clima_own <- function (nc, lsm_path, long_min, long_max, lat_min, lat_max, start_time, 
          end_time, dtr_cor = TRUE, dtr_cor_fac = 1.285, format = "microclimf") 
{
  nc_dat = ncdf4::nc_open(nc)
  timedim <- extract_timedim(nc_dat)
  timedim$units <- "seconds since 1970-01-01"
  base_datetime <- as.POSIXct(gsub(".*since ", "", timedim$units), 
                              tz = "UTC")
  nc_datetimes <- c(timedim$vals)
  nc_datetimes <- nc_datetimes * ifelse(grepl("hours", timedim$units), 
                                        3600, 1)
  first_timestep <- nc_datetimes[1]
  last_timestep <- utils::tail(nc_datetimes, n = 1)
  if (any(!class(start_time) %in% c("Date", "POSIXct", "POSIXt", 
                                    "POSIXlt")) | any(!class(end_time) %in% c("Date", "POSIXct", 
                                                                              "POSIXt", "POSIXlt"))) {
    stop("`start_time` and `end_time` must be provided as date-time objects.")
  }
  if (any(class(start_time) != class(end_time))) {
    stop("`start_time` and `end_time` must be of the same date-time class.")
  }
  start <- base_datetime + first_timestep
  if (start_time < start) {
    stop("Requested start time is before the beginning of time series of the ERA5 netCDF.")
  }
  end <- base_datetime + last_timestep
  if (end_time > end) {
    stop("Requested end time is after the end of time series of the ERA5 netCDF.")
  }
  if (long_max <= long_min) {
    stop("Maximum longitude must be greater than minimum longitude.")
  }
  if (lat_max <= lat_min) {
    stop("Maximum longitude must be greater than minimum longitude.")
  }
  if (abs(long_min) > 180 | abs(long_max) > 180 | abs(lat_min) > 
      90 | abs(lat_max) > 90) {
    stop("Coordinates must be provided in decimal degrees (longitude between -180 and 180, latitude between -90 and 90).")
  }
  if (long_min < (min(nc_dat$dim$lon$vals) - 0.125) | 
      long_max > (max(nc_dat$dim$lon$vals) + 0.125)) {
    long_out <- TRUE
  } else {
    long_out <- FALSE
  }
  if (lat_min < (min(nc_dat$dim$lat$vals) - 0.125) | lat_min > 
      (max(nc_dat$dim$lat$vals) + 0.125)) {
    lat_out <- TRUE
  } else {
    lat_out <- FALSE
  }
  ncdf4::nc_close(nc_dat)
  if (long_out & lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (both longitude and latitude out of range).")
  }
  if (long_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (longitude out of range).")
  }
  if (lat_out) {
    stop("Requested coordinates are not represented in the ERA5 netCDF (latitude out of range).")
  }
  if (lubridate::tz(start_time) != lubridate::tz(end_time)) {
    stop("start_time and end_time are not in the same timezone.")
  }
  if (lubridate::tz(start_time) != "UTC" | lubridate::tz(end_time) != 
      "UTC") {
    warning("provided times (start_time and end_time) are not in timezone UTC (default timezone of ERA5 data). Output will be provided in timezone UTC however.")
  }
  if (lubridate::hour(end_time) == 0) {
    end_time <- as.POSIXlt(paste0(lubridate::year(end_time), 
                                  "-", lubridate::month(end_time), "-", lubridate::day(end_time), 
                                  " 23:00"), tz = lubridate::tz(end_time))
  }
  if (!format %in% c("NicheMapR", "microclima", "microclimc", 
                     "micropoint", "microclimf")) {
    stop("Argument `format` must be one of the following values: `NicheMapR`, `microclima`, `microclimc`, `micropoint`, `microclimf`")
  }
  if (dtr_cor == TRUE & !is.numeric(dtr_cor_fac)) {
    stop("Invalid diurnal temperature range correction value provided.")
  }
  tme <- as.POSIXct(seq(start_time, end_time, by = 3600), tz = lubridate::tz(end_time))
  nc_dat <- ncdf4::nc_open(nc)
  varnames <- names(nc_dat$var)
  varname_list <- c("t2m","d2m","sp",
                    "u10","v10","tp",
                    "ssrd","fdir","avg_snlwrf",
                    "avg_sdlwrf","tcc")
  # v <- "u10"
  var_list <- lapply(varname_list, function(v) {
    v2 <- paste0("var_", which(varname_list == v))
    
    r <- terra::rast(nc, subds = v2)
    varnames(r) <- v
    r <- r[[as.POSIXct(nc_datetimes, tz = "UTC") %in% 
              tme]]
    names(r) <- tme
    
    r <- terra::crop(r, terra::ext(long_min, long_max, lat_min, 
                                   lat_max))
    return(r)
  })
  
  names(var_list) <- varname_list
  if ("msnlwrf" %in% names(var_list)) {
    var_list$avg_snlwrf <- var_list$msnlwrf
    var_list$avg_sdlwrf <- var_list$msdwlwrf
  }
  t2m <- var_list$t2m
  d2m <- var_list$d2m
  sp <- var_list$sp
  u10 <- var_list$u10
  v10 <- var_list$v10
  tcc <- var_list$tcc
  avg_snlwrf <- var_list$avg_snlwrf
  avg_sdlwrf <- var_list$avg_sdlwrf
  fdir <- var_list$fdir
  ssrd <- var_list$ssrd
  prec <- var_list$tp * 1000
  
  lsm <- rast(lsm_path)
  lsm <- terra::project(lsm, t2m)
  
  temperature <- t2m - 273.15
  if (any(terra::values(lsm) < 1)) {
    ind <- rep(1:(dim(temperature)[3]/24), each = 24)
    tmean <- terra::tapp(temperature, ind, fun = mean, na.rm = T)
    tmean <- rep(tmean, 24)
    tmean <- tmean[[paste0("X", sort(rep(seq(1:(dim(temperature)[3]/24)), 
                                         24)))]]
    m <- (1 - lsm) * dtr_cor_fac + 1
    tdif <- (temperature - tmean) * m
    temperature <- tmean + tdif
  }
  humidity <- humfromdew(d2m - 273.15, temperature, sp)
  windspeed = sqrt(u10^2 + v10^2)
  windspeed = windheight(windspeed, 10, 2)
  winddir = (terra::atan2(u10, v10) * 180/pi + 180)%%360
  cloudcover = tcc * 100
  netlong = abs(avg_snlwrf) * 0.0036
  downlong = avg_sdlwrf * 0.0036
  uplong = netlong + downlong
  emissivity = downlong/uplong
  jd = julday(lubridate::year(tme), lubridate::month(tme), 
              lubridate::day(tme))
  rad_dni = fdir * 1e-06
  rad_glbl = ssrd * 1e-06
  si <- t2m
  foo <- t2m[[1]]
  coords <- as.data.frame(terra::crds(t2m[[1]]))
  out <- array(NA, dim = c(nrow(coords), length(tme)))
  for (i in 1:nrow(coords)) {
    out[i, ] <- siflat(lubridate::hour(tme), lat = coords$y[i], 
                       long = coords$x[i], jd)
  }
  si <- terra::setValues(si, out)
  rad_dif = rad_glbl - rad_dni * si
  szenith <- t2m
  out <- array(NA, dim = c(nrow(coords), length(tme)))
  for (i in 1:nrow(coords)) {
    out[i, ] <- solalt(lubridate::hour(tme), lat = coords$y[i], 
                       long = coords$x[i], julian = jd)
  }
  szenith <- terra::setValues(szenith, out)
  names(temperature) <- names(t2m)
  terra::time(temperature) <- tme
  if (format %in% c("microclimc", "micropoint", "microclimf")) {
    pres <- sp/1000
    relhum <- humidity
    terra::values(relhum) <- converthumidity(h = terra::as.array(humidity), 
                                             intype = "specific", tc = terra::as.array(temperature), 
                                             pk = terra::as.array(pres))
    relhum[relhum > 100] <- 100
    message("note: extract_clima() does not currently convert accumulated measurements across the hour to a mean on the hour, as is performed in extract_clim(). This yields small but non-marginable differences in values.")
    raddr <- (rad_dni * si)/0.0036
    difrad <- rad_dif/0.0036
    swrad <- raddr + difrad
    downlong <- downlong/0.0036
  }
  if (format %in% c("micropoint", "microclimf")) {
    return(list(obs_time = tme, temp = terra::wrap(temperature), 
                relhum = terra::wrap(relhum), pres = terra::wrap(pres), 
                swdown = terra::wrap(swrad), difrad = terra::wrap(difrad), 
                lwdown = terra::wrap(downlong), windspeed = terra::wrap(windspeed), 
                winddir = terra::wrap(winddir), precip = terra::wrap(prec)))
  }
  if (format == "microclimc") {
    return(list(obs_time = tme, temp = terra::wrap(temperature), 
                relhum = terra::wrap(relhum), pres = terra::wrap(pres), 
                swrad = terra::wrap(swrad), difrad = terra::wrap(difrad), 
                skyem = terra::wrap(emissivity), windspeed = terra::wrap(windspeed), 
                winddir = terra::wrap(winddir), precip = terra::wrap(prec)))
  }
  if (format %in% c("microclima", "NicheMapR")) {
    return(list(obs_time = tme, temperature = terra::wrap(temperature), 
                humidity = terra::wrap(humidity), pressure = terra::wrap(sp), 
                windspeed = terra::wrap(windspeed), winddir = terra::wrap(winddir), 
                emissivity = terra::wrap(emissivity), netlong = terra::wrap(netlong), 
                uplong = terra::wrap(uplong), downlong = terra::wrap(downlong), 
                rad_dni = terra::wrap(rad_dni), rad_dif = terra::wrap(rad_dif), 
                szenith = terra::wrap(szenith), cloudcover = terra::wrap(tcc)))
  }
}