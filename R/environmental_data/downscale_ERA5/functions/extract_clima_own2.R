# climhrly <- climdata
temphrly_downscale <- function (climhrly, sst, dtmf, dtmm = NA, basins = NA, uzf = NA, 
          cad = TRUE, coastal = TRUE, thgto = 2, whgto = 2, tempvar = "temp") 
{
  input_class <- lapply(lapply(climhrly, class), "[", 1)
  if (any(input_class == "PackedSpatRaster")){
    climhrly[which(input_class == "PackedSpatRaster")] <- lapply(climhrly[which(input_class == 
                                                                                  "PackedSpatRaster")], unwrap)
  }
  if (any(input_class == "array")){
    climhrly[which(input_class == "array")] <- lapply(climhrly[which(input_class == 
                                                                       "array")], .rast, tem = climhrly$dtm)
  }
  if (inherits(sst, "PackedSpatRaster")){
    sst <- unwrap(sst)
  }
    
  rh <- climhrly$relhum
  pk <- climhrly$pres
  tc <- climhrly[[tempvar]]
  thgti <- climhrly$tempheight_m
  whgti <- climhrly$windheight_m
  tme <- climhrly$tme
  u2 <- .windhgt(climhrly$windspeed, whgti, thgto)
  swrad <- climhrly$swrad
  lwrad <- climhrly$lwrad
  dtmc <- climhrly$dtm
  if (inherits(sst, "SpatRaster")){
    if (all(is.na(values(sst[[1]])))){
      sst <- NA
    }
  }  
  if (is.logical(sst)) {
    if (coastal == TRUE){
      message("Ignoring coastal effects as no valid sea surface temperature data provided - assuming inland area!!")
    }
    coastal <- FALSE
  }
  if (is.numeric(sst)) {
    if (length(sst) != length(tme)) {
      warning("SST provided as vector is NOT the same length as timesteps/layers in temperature data - ignoring coastal effects!!!")
      coastal <- FALSE
    } else {
      sstvals <- rep(sst, each = ncell(dtmc))
      sst <- rast(crs = crs(dtmc), extent = ext(dtmc), 
                  resolution = res(dtmc), nlyr = length(tme), vals = sstvals)
      terra::time(sst) <- tme
    }
  }
  lrc <- lapserate(.is(tc), .is(rh), .is(pk))
  lrc <- .rast(lrc, dtmc)
  if (crs(lrc) != crs(dtmf)) {
    lrcp <- project(lrc, crs(dtmf))
    lrf <- .resample(lrcp, dtmf)
  } else {
    lrf <- .resample(lrc, dtmf)
  }
  tcf <- .tempelev(tc, lrf, lrc, dtmf + thgto, ifel(is.na(dtmc), 
                                                    0, dtmc + thgti))
  if (cad) {
    mu <- .cad_multiplier(dtmf, basins, refhgt = thgti)
    st <- .cad_conditions(tc, u2, swrad = 0, lwrad, dtmc, 
                          dtmf, refhgt = thgti)
    tcad <- .apply_cad(lrf, mu, st)
    tcf <- tcf + tcad
  }
  if (coastal) {
    if (any(global(sst, anyNA))){
      sstinterp <- .spatinterp(sst)
    } else {
      sstinterp <- sst
    } 
    sstinterp <- .tmeinterp(sstinterp, NA, tme)
    if (crs(sst) != crs(dtmf)){
      sstinterp <- project(sstinterp, crs(dtmf))
    }
    sstf <- .resample(sstinterp, dtmf, method = "cubic")
    if (class(uzf)[1] == "logical") 
      uzf <- winddownscale(climhrly$windspeed, climhrly$winddir, 
                           dtmf, dtmm, dtmc, whgti, thgto)
    tcf <- .tempcoastal(tc = tcf, sst = sstf, u2 = uzf, wdir = climhrly$winddir, 
                        dtmf, dtmm, dtmc)
  }
  terra::time(tcf) <- tme
  return(tcf)
}

swdownscale <- function (swrad, tme = NA, dtmf, dtmc, patchsim = FALSE, nsim = dim(swrad)[3], 
          hor = NA, terrainshade = FALSE) 
{
  if (inherits(swrad, "SpatRaster") & inherits(tme, "logical")){
    tme <- as.POSIXlt(terra::time(swrad))
  }
  if (inherits(tme, "logical")){
    stop("swdownscale requires tme parameter or spatraster with time!!!")
  }
    
  ti <- round(as.numeric(tme[2]) - as.numeric(tme[1]))
  if (ti <= 3600){
    hourly <- TRUE
  } else {
    hourly <- FALSE
  } 
  dtmc <- ifel(is.na(dtmc), 0, dtmc)
  jd <- juldayR(tme$year + 1900, tme$mon + 1, tme$mday)
  lt <- tme$hour + tme$min/60 + tme$sec/3600
  ll <- .latslonsfromr(dtmc)
  lats <- .rast(ll$lats, dtmc)
  lons <- .rast(ll$lons, dtmc)
  lats <- .is(mask(lats, dtmc))
  lons <- .is(mask(lons, dtmc))
  csr <- array(clearskyradmR(jd, lt, as.vector(lats), as.vector(lons), 
                               hourly), dim = dim(swrad))
  csfc <- .is(swrad)/csr
  csfc[is.na(csfc)] <- 0.5
  csfc[csfc > 1] <- 1
  csfc[csfc < 0] <- 0
  csfc <- .rast(csfc, dtmc)
  csf <- .cfcelev(csfc, dtmf, dtmc)
  if (patchsim) {
    if (crs(dtmc) != crs(dtmf)) {
      dtmp <- project(dtmc, crs(dtmf))
    } else {
      dtmp <- dtmc
    } 
    af <- res(dtmp)[1]/res(dtmf)[1]
    csf <- .simpatch(csf, af, n = nsim, varn = "csfrac")
  }
  ll <- .latslonsfromr(dtmf)
  lats <- .rast(ll$lats, dtmf)
  lons <- .rast(ll$lons, dtmf)
  lats <- .is(mask(lats, dtmf))
  lons <- .is(mask(lons, dtmf))
  if (hourly) {
    lats <- .rast(ll$lats, dtmf)
    lons <- .rast(ll$lons, dtmf)
    lats <- .is(mask(lats, dtmf))
    lons <- .is(mask(lons, dtmf))
    csrf <- array(clearskyradmR(jd, lt, as.vector(lats), 
                                  as.vector(lons)), dim = c(dim(dtmf)[1:2], dim(swrad)[3]))
  } else {
    lats <- as.vector(ll$lats[, 1])
    csrf <- clearskyradmR(jd, rep(12, length(jd)), lats, 
                            rep(0, length(lats)), hourly = FALSE)
    csrf <- array(apply(csrf, 2, rep, times = dim(dtmf)[2]), 
                  dim = c(dim(dtmf)[1:2], dim(swrad)[3]))
  }
  swradf <- .rast(csrf, dtmf) * csf
  if (terrainshade) {
    if (!hourly){
      swradfh <- .ehr(.is(swradf))
    } else {
      swradfh <- .is(swradf)
    } 
    if (hourly) {
      swf <- matrix(.is(swradf), ncol = dim(swradf)[3])
      dp <- difpropmR(swf, jd, lt, as.vector(lats), as.vector(lons))
      dp <- array(dp, dim = dim(swradf))
      ze <- array(solzenmR(jd, lt, as.vector(lats), as.vector(lons)), 
                  dim = c(dim(dtmf)[1:2], dim(swrad)[3]))
      azi <- .solazi(lt, mean(lats, na.rm = TRUE), mean(lons, 
                                                        na.rm = TRUE), jd)
    } else {
      csrfh <- .ehr(csrf)
      lth <- rep(c(0:23), dim(csrf)[3])
      jdh <- rep(jd, each = 24)
      ll <- .latslonsfromr(dtmf)
      lats <- .rast(ll$lats, dtmf)
      lons <- .rast(ll$lons, dtmf)
      lats <- .is(mask(lats, dtmf))
      lons <- .is(mask(lons, dtmf))
      cs <- array(clearskyradmR(jdh, lth, as.vector(lats), 
                                  as.vector(lons)), dim = c(dim(dtmf)[1:2], length(jdh)))
      rad <- cs * csrfh
      swf <- matrix(rad, ncol = length(jdh))
      dp <- difpropmR(swf, jdh, lth, as.vector(lats), 
                        as.vector(lons))
      dp <- array(dp, dim = dim(csrfh))
      ze <- array(solzenmR(jdh, lth, as.vector(lats), 
                             as.vector(lons)), dim = c(dim(dtmf)[1:2], dim(csrfh)[3]))
      azi <- .solazi(lth, mean(lats, na.rm = TRUE), mean(lons, 
                                                         na.rm = TRUE), jdh)
    }
    if (inherits(hor, "logical")) {
      hor <- array(NA, dim = c(dim(dtmf)[1:2], 24))
      for (i in 1:24) hor[, , i] <- .horizon(dtmf, (i - 
                                                      1) * 15)
    } else {
      hor <- .is(hor)
    } 
    nhor <- dim(hor)[3]
    i <- round(azi/(360/nhor)) + 1
    i[i == (nhor + 1)] <- 1
    hora <- hor[, , i]
    shadowmask <- hora * 0 + 1
    alt <- (90 - ze) * pi/180
    shadowmask[hora > tan(alt)] <- 0
    shadowmask[(90 - ze) < 0] <- 0
    svf <- .rta(.skyview(dtmf), dim(shadowmask)[3])
    drf <- dp * svf * swradfh
    swf <- (1 - dp) * shadowmask * swradfh + dp * svf * swradfh
    drf <- .rast(drf, dtmf)
    swf <- .rast(swf, dtmf)
    if (!hourly) {
      swf <- hourtodayR(.is(swf), "mean")
      drf <- hourtodayR(.is(drf), "mean")
      swf <- .rast(swf, dtmf)
      drf <- .rast(drf, dtmf)
    }
    terra::time(swf) <- tme
    terra::time(drf) <- tme
    out <- list(swf = swf, drf = drf)
  } else {
    terra::time(swradf) <- tme
    out <- list(swf = swradf, drf = NA)
  }
  return(out)
}

.skyview<-function(dtm,steps=36) {
  r<-dtm
  dtm[is.na(dtm)]<-0
  ha <- array(0, dim(dtm)[1:2])
  for (s in 1:steps) { # uses horizon angle in calc but places equal importance on each sector of sky
    ha<-ha+atan(.horizon(dtm,s*360/steps))
  }
  ha<-ha/steps
  ha<-tan(ha)
  svf<-0.5*cos(2*ha)+0.5
  svf<-.rast(svf,dtm)
  svf<-mask(svf,r)
  return(svf)
}

.solazi <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  hh<-0.5*pi-acos(sh)
  sazi<-cos(d)*sin(tt)/cos(hh)
  cazi<-(sin(lat*pi/180)*cos(d)*cos(tt)-cos(lat*pi/180)*sin(d))/
    sqrt((cos(d)*sin(tt))^2+(sin(lat*pi/180)*cos(d)*cos(tt)-cos(lat*pi/180)*sin(d))^2)
  sqt <- 1-sazi^2
  sqt[sqt<0]<-0
  solz<-180+(180*atan(sazi/sqrt(sqt)))/pi
  solz[cazi<0 & sazi<0]<-180-solz[cazi<0 & sazi<0]
  solz[cazi<0 & sazi>=0]<-540-solz[cazi<0 & sazi>=0]
  solz
}

.horizon <- function(dtm, azimuth) {
  reso<-res(dtm)[1]
  dtm<-.is(dtm)
  dtm[is.na(dtm)]<-0
  dtm<-dtm/reso
  azi<-azimuth*pi/180
  horizon<-array(0,dim(dtm))
  dtm3<-array(0,dim(dtm)+200)
  x<-dim(dtm)[1]
  y<-dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)]<-dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2),na.rm=T)
  }
  horizon
}

.simpatch<-function(varf,af,n=dim(varf)[3],varn="precip") {
  if (varn == "precip") {
    varf<-log(varf+1)
    sdr<- -0.14934+0.1777*log(af)
    if (sdr> 1) sdr<-1
  }
  if (varn == "csfrac") {
    vf<-.is(varf)
    s0<-which(vf==0)
    s1<-which(vf==1)
    vf<-suppressWarnings(log(vf/(1-vf)))
    vf[is.infinite(vf)]<-NA
    vf[s0]<-min(vf,na.rm=TRUE)
    vf[s1]<-max(vf,na.rm=TRUE)
    varf<-.rast(vf,varf)
    sdr<- 0.012950+0.175718*log(af)
    if (sdr> 1) sdr<-1
  }
  # Generate n simulations
  xy <- expand.grid(1:dim(varf)[2], 1:dim(varf)[1])
  names(xy) <- c('x','y')
  g1<-gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=0,
            model=gstat::vgm(psill=1,
                      range=af/2,
                      nugget=0,
                      model='Sph'), nmax = 40)
  ni<-dim(varf)[3]
  hps<-round(ni/(n-1),0)
  if (n != ni) {
    days<-ni/hps
  } else days <- ni+1
  yy1 <- predict(g1,newdata=xy,nsim=days-1,debug.level=0)
  anom<-.rast(yy1,varf[[1]])
  anom<-resample(aggregate(anom,2),anom)
  anom<-mask(anom,varf[[1]])
  ni<-dim(varf)[3]
  # if n not equal to ni generate ni anomoly layers by averaging
  # allows for temporal autocorrelation in simulations
  if (n != ni) {
    ad<-.is(anom)
    wgts<-c(0:(hps-1))/hps
    am<-array(NA,dim=dim(varf))
    ij<-1
    for (i in 1:days) {
      for (j in 1:hps) {
        am[,,ij]<-(1-wgts[j])*ad[,,i]+wgts[j]*ad[,,ifelse(i == days, i, i+1)]
        ij<-ij+1
      }
    }
    anom<-.rast(am,anom[[1]])
  }
  # adjust to ensure anomaly centered on zero and with correct sd
  me<-apply(.is(anom),3,mean,na.rm=T)
  sds<-apply(.is(anom),3,sd,na.rm=T)
  sda<-apply(.is(varf),3,sd,na.rm=T)
  mu <- (sda*sdr)/sds
  a<-.is(anom)
  a<-(a-.vta(me,varf[[1]]))*.vta(mu,varf[[1]])
  anom<-.rast(a,anom[[1]])
  varnf<-varf+anom
  if (varn == "precip") {
    varnf<-exp(varnf)-1
    varnf[varnf<0]<-0
  }
  if (varn == "csfrac") {
    varnf<-1/(1+exp(-varnf))
  }
  return(varnf)
}

difpropmR <- function(swrad, jd, lt, lat, lon) {
  # This function calculates the diffuse proportion of global radiation for a matrix of locations and times.
  # It is a direct translation of the provided C++ Rcpp code.
  # It assumes a function 'difpropvR' exists.
  
  # Check for correct dimensions
  if (nrow(swrad) != length(lat) || ncol(swrad) != length(jd)) {
    stop("Dimensions of 'swrad' matrix do not match the lengths of 'lat' and 'jd'.")
  }
  
  # Get matrix dimensions
  rows <- length(lat)
  cols <- length(jd)
  
  # Initialize a list to store the results
  d_list <- vector("list", rows)
  
  # Loop through each location (row in the swrad matrix)
  for (i in 1:rows) {
    # Check if the second value for the current location's radiation vector is not NA
    # The C++ code checks swradc[i][1] (second element, 0-indexed).
    if (!is.na(swrad[i, 2])) {
      # Pass the entire row of radiation data to the vector function
      d_list[[i]] <- difpropvR(swrad[i, ], jd, lt, lat[i], lon[i])
    } else {
      # If the check fails, fill the corresponding row with NAs
      d_list[[i]] <- rep(NA_real_, cols)
    }
  }
  
  # Convert the list of vectors into a matrix and return it
  dr <- do.call(rbind, d_list)
  
  return(dr)
}

difpropvR <- function(swrad, jd, lt, lat, lon) {
  # This function calculates the diffuse proportion for each value in a vector of solar radiation measurements.
  # It is a direct translation of the provided C++ code.
  # It assumes a function 'difpropR' exists.
  
  # Check if the input vectors have the same length
  if (length(swrad) != length(jd) || length(swrad) != length(lt)) {
    stop("Input vectors 'swrad', 'jd', and 'lt' must have the same length.")
  }
  
  # Use mapply for efficient vectorization. It applies the 'difpropR'
  # function to corresponding elements of the input vectors.
  d <- mapply(difpropR, swrad, jd, lt,
              MoreArgs = list(lat = lat, lon = lon))
  
  return(d)
}

difpropR <- function(swrad, jd, lt, lat, lon) {
  # This function calculates the diffuse proportion of global solar radiation.
  # It is a direct translation of the provided C++ code.
  # It assumes 'soltimeR' and 'solzenR' functions exist.
  
  pi <- base::pi
  d <- 1.0 # Initialize diffuse proportion to 1 (totally diffuse)
  
  # Only perform calculations if there is solar radiation
  if (swrad > 0) {
    # Get solar time and solar zenith angle
    st <- soltimeR(jd, lt, lon)
    z <- solzenR(jd, st, lat, lon)
    
    # Only calculate if the sun is above the horizon
    if (z < pi / 2) {
      # Convert zenith angle to degrees
      zd <- z * 180 / pi
      
      # Calculate empirical parameters
      k1 <- 0.83 - 0.56 * exp(-0.06 * (90 - zd))
      
      # Calculate the cosine of the zenith angle for the si variable
      si <- 0.0
      if (z <= pi / 2) si <- cos(z)
      
      # Calculate the clearness index (k)
      k <- 0.0
      if (si > 0) k <- swrad / (1352 * si)
      
      # Clamp k to a reasonable range
      if (k > k1) k <- k1
      if (k < 0) k <- 0
      
      rho <- k / k1
      sigma3 <- 0
      
      # Calculate sigma3 based on rho
      if (rho > 1.04) {
        sigma3 <- 0.12 + 0.65 * (rho - 1.04)
      } else {
        sigma3 <- 0.021 + 0.397 * rho - 0.231 * rho^2 - 0.13 * exp(-1 * ((rho - 0.931) / 0.134)^2 * 0.834)
      }
      
      k2 <- 0.95 * k1
      d1 <- 1.0
      
      # Calculate d1 based on zenith angle
      if (zd < 88.6) d1 <- 0.07 + 0.046 * zd / (93 - zd)
      
      K <- 0.5 * (1 + sin(pi * (k - 0.22) / (k1 - 0.22) - pi / 2))
      d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 * K^2))
      d3 <- (d2 * k2) * (1 - k) / (k * (1 - k2))
      
      alpha <- (1 / cos(z))^0.6
      kbmax <- 0.81^alpha
      kmax <- (kbmax + d2 * k2 / (1 - k2)) / (1 + d2 * k2 / (1 - k2))
      dmax <- (d2 * k2) * (1 - kmax) / (kmax * (1 - k2))
      
      # Determine the diffuse proportion based on k and other parameters
      d <- 1 - kmax * (1 - dmax) / k
      
      if (k <= kmax) d <- d3
      if (k <= k2) d <- d2
      if (k <= 0.22) d <- 1
      
      kX <- 0.56 - 0.32 * exp(-0.06 * (90 - zd))
      kL <- (k - 0.14) / (kX - 0.14)
      kR <- (k - kX) / 0.71
      
      delta <- 0
      if (k >= 0.14 && k < kX) {
        delta <- -3 * kL^2 * (1 - kL) * sigma3^1.3
      } else if (k >= kX && k < (kX + 0.71)) {
        delta <- 3 * kR * (1 - kR)^2 * sigma3^0.6
      }
      
      if (sigma3 > 0.01) d <- d + delta
    }
  }
  return(d)
}

hourtodayR <- function(a, fun = "mean") {
  # This function aggregates hourly data in a 3D array to daily data.
  # It is a direct translation of the provided C++ Rcpp code.
  
  # Check if the input is a 3D array
  if (length(dim(a)) != 3) {
    stop("Input 'a' must be a 3D array.")
  }
  
  dims <- dim(a)
  dim1 <- dims[1]
  dim2 <- dims[2]
  dim3 <- dims[3]
  
  # Ensure the third dimension is a multiple of 24
  if (dim3 %% 24 != 0) {
    stop("The third dimension of the array must be a multiple of 24.")
  }
  
  ndays <- dim3 / 24
  
  # Reshape the 3D array into a 4D array for easier day-wise operations
  # The new dimensions will be (dim1, dim2, 24, ndays)
  a_reshaped <- array(a, dim = c(dim1, dim2, 24, ndays))
  
  # Initialize the output array with the correct dimensions
  daily <- array(NA_real_, dim = c(dim1, dim2, ndays))
  
  # Apply the specified function to the data
  if (fun == "mean") {
    daily <- apply(a_reshaped, c(1, 2, 4), mean, na.rm = TRUE)
  } else if (fun == "sum") {
    daily <- apply(a_reshaped, c(1, 2, 4), sum, na.rm = TRUE)
  } else if (fun == "max") {
    daily <- apply(a_reshaped, c(1, 2, 4), max, na.rm = TRUE)
  } else if (fun == "min") {
    daily <- apply(a_reshaped, c(1, 2, 4), min, na.rm = TRUE)
  } else {
    stop("Invalid function specified. Choose 'mean', 'sum', 'max', or 'min'.")
  }
  
  return(daily)
}

.rta <- function(r,n) {
  m<-.is(r)
  a<-array(rep(m,n),dim=c(dim(r)[1:2],n))
  a
}

juldayR <- function(year, month, day) {
  # This function calculates the Julian Day number from a given date.
  # The formula is based on the algorithm found in "Numerical Recipes".
  
  # Adjust month and year for the formula
  madj <- month + (month < 3) * 12
  yadj <- year + (month < 3) * -1
  
  # Calculate the Julian Day
  j <- floor(365.25 * (yadj + 4716)) + floor(30.6001 * (madj + 1)) + day - 1524.5
  
  # Adjust for the Gregorian calendar
  b <- 2 - floor(yadj / 100) + floor(floor(yadj / 100) / 4)
  
  # Apply the Gregorian adjustment if the date is after the adoption date
  jd <- j + (j > 2299160) * b
  
  # Return the integer part of the Julian Day
  return(as.integer(jd))
}

clearskyradhR <- function(jd, lt, lat, lon) {
  # This function calculates hourly clear-sky radiation for a vector of Julian days and local times.
  
  # Corrected call to clearskyradR with the right argument names (tc, rh, pk).
  Ic <- mapply(clearskyradR, jd, lt, MoreArgs = list(lat = lat, lon = lon,
                                                     tc = 15.0, rh = 80.0, pk = 101.3))
  
  return(Ic)
}

clearskyraddR <- function(jd, lat) {
  # This function calculates the average daily clear-sky radiation for a vector of Julian days.
  
  minute_times <- seq(0, 1439) / 60
  
  Ic <- sapply(jd, function(current_jd) {
    # Corrected call to clearskyradR with the right argument names (tc, rh, pk).
    minute_rads <- mapply(clearskyradR, current_jd, minute_times, MoreArgs = list(lat = lat, lon = 0.0,
                                                                                  tc = 15.0, rh = 80.0, pk = 101.3))
    
    return(sum(minute_rads) / 1440)
  })
  
  return(Ic)
}

clearskyradR <- function(jd, lt, lat, lon, tc = 15.0, rh = 80.0, pk = 101.3) {
  # This function calculates the clear-sky radiation based on the C++ algorithm.
  # It assumes that soltimeR and solzenR functions exist.
  
  # The value of M_PI in C++ is equivalent to pi in R.
  pi <- base::pi
  
  # Calculate solar time and solar zenith angle.
  # Assuming soltimeR and solzenR are the R equivalents of the C++ functions.
  st <- soltimeR(jd, lt, lon)
  z <- solzenR(jd, st, lat, lon)
  
  # Check if the solar zenith angle is less than or equal to 90 degrees (pi/2 radians).
  # If it is, proceed with the radiation calculation. Otherwise, radiation is 0.
  if (z <= pi / 2) {
    # Calculate the air mass ratio (m)
    m <- 35 * cos(z) * (1224 * cos(z)^2 + 1)^(-0.5)
    
    # Calculate the Transmittance of Rayleigh scattering and permanent gases (TrTpg)
    TrTpg <- 1.021 - 0.084 * sqrt(m * 0.00949 * pk + 0.051)
    
    # Calculate the dew point temperature (Td)
    xx <- log(rh / 100) + ((17.27 * tc) / (237.3 + tc))
    Td <- (237.3 * xx) / (17.27 - xx)
    
    # Calculate the equivalent atmospheric water vapor content (u)
    u <- exp(0.1133 - log(3.78) + 0.0393 * Td)
    
    # Calculate the Transmittance of water vapour (Tw)
    Tw <- 1 - 0.077 * (u * m)^0.3
    
    # Calculate the Transmittance of aerosols (Ta)
    Ta <- 0.935^m
    
    # Calculate the total atmospheric optical depth (od)
    od <- TrTpg * Tw * Ta
    
    # Calculate the clear-sky radiation (Ic)
    Ic <- 1352.778 * cos(z) * od
  } else {
    # When the sun is below the horizon, clear-sky radiation is zero.
    Ic <- 0.0
  }
  
  return(Ic)
}

clearskyradmR <- function(jd, lt, lat, lon, hourly = TRUE) {
  # This function calculates a matrix of clear-sky radiation for multiple locations and dates.
  # It is a direct translation of the C++ function and uses the 'clearskyradhR' and 'clearskyraddR' functions.
  
  if (length(lat) != length(lon)) {
    stop("Vectors 'lat' and 'lon' must have the same length.")
  }
  
  # Pre-allocate a list to store the results for each location.
  Ic_list <- vector("list", length(lat))
  
  # Iterate over each location (latitude and longitude pair).
  for (i in seq_along(lat)) {
    # Check for NaN values, similar to the C++ code's isnan check.
    if (!is.na(lat[i])) {
      if (hourly) {
        Ic_list[[i]] <- clearskyradhR(jd, lt, lat[i], lon[i])
      } else {
        Ic_list[[i]] <- clearskyraddR(jd, lat[i])
      }
    } else {
      # If latitude is NA, fill the row with NAs.
      if (hourly) {
        Ic_list[[i]] <- rep(NA_real_, length(jd))
      } else {
        Ic_list[[i]] <- rep(NA_real_, length(jd))
      }
    }
  }
  
  # Convert the list of vectors into a matrix, with each vector as a column.
  # Use 'do.call(rbind, ...)' to combine the vectors into a matrix.
  # The C++ code's `NumericMatrix ICm` and `convertoRmatrix` likely transpose the matrix,
  # so `t(do.call(rbind, Ic_list))` is the most accurate translation.
  # The number of rows is lat.size() and cols is jd.size(), so rbind is the right approach.
  result_matrix <- do.call(rbind, Ic_list)
  
  return(result_matrix)
}

clearskyradmR <- function(jd, lt, lat, lon, hourly = TRUE) {
  # This function calculates a matrix of clear-sky radiation values.
  # It is a direct translation of the provided C++ Rcpp code.
  
  # Check if the latitude and longitude vectors have the same length
  if (length(lat) != length(lon)) {
    stop("Vectors 'lat' and 'lon' must have the same length.")
  }
  
  # Determine the dimensions of the final matrix
  rows <- length(lat)
  cols <- length(jd)
  
  # Initialize a list to store the results for each location
  Ic_list <- vector("list", rows)
  
  # Loop through each location (latitude/longitude pair)
  for (i in 1:rows) {
    # Check if the latitude is not a missing value (NaN in C++, NA in R)
    if (!is.na(lat[i])) {
      if (hourly) {
        # Call the hourly calculation function for this location
        Ic_list[[i]] <- clearskyradhR(jd, lt, lat[i], lon[i])
      } else {
        # Call the daily calculation function for this location
        Ic_list[[i]] <- clearskyraddR(jd, lat[i])
      }
    } else {
      # If latitude is NA, fill the corresponding row with NAs
      Ic_list[[i]] <- rep(NA_real_, cols)
    }
  }
  
  # Convert the list of vectors into a matrix
  # The C++ code's `NumericMatrix` is typically row-major,
  # meaning each row corresponds to a location.
  result_matrix <- do.call(rbind, Ic_list)
  
  return(result_matrix)
}

soltimeR <- function(jd, lt, lon) {
  # Calculate solar time from Julian day, local time, and longitude.
  # The formula approximates the Equation of Time (eot).
  
  m <- 6.24004077 + 0.01720197 * (jd - 2451545)
  eot <- -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932)
  st <- lt + (4 * lon + eot) / 60
  
  return(st)
}

daylengthR <- function(jd, lat) {
  # Calculate the length of the day in hours.
  
  pi <- base::pi
  
  # Calculate solar declination in radians.
  declin <- (pi * 23.5 / 180) * cos(2 * pi * ((jd - 159.5) / 365.25))
  
  # Convert latitude to radians.
  latr <- lat * pi / 180
  
  # Calculate hour angle term (hc).
  hc <- -0.01453808 / (cos(latr) * cos(declin)) - tan(latr) * tan(declin)
  
  # Handle special cases for polar regions.
  if (hc < -1) {
    dl <- 24  # 24 hours of daylight (polar day).
  } else if (hc > 1) {
    dl <- 0   # 0 hours of daylight (polar night).
  } else {
    # Calculate hour angle (ha) and use the equation of time for sunrise/sunset times.
    ha <- (acos(hc)) * 180 / pi
    m <- 6.24004077 + 0.01720197 * (jd - 2451545)
    eot <- -7.659 * sin(m) + 9.863 * sin(2 * m + 3.5932)
    sr <- (720 - 4 * ha - eot) / 60  # Sunrise time in hours
    ss <- (720 + 4 * ha - eot) / 60  # Sunset time in hours
    dl <- ss - sr
  }
  
  return(dl)
}

solzenR <- function(jd, st, lat, lon) {
  # Calculate the solar zenith angle in radians.
  
  pi <- base::pi
  
  # Convert latitude to radians.
  latr <- lat * pi / 180
  
  # Calculate hour angle term (tt).
  tt <- 0.261799 * (st - 12)
  
  # Calculate solar declination in radians.
  dec <- (pi * 23.5 / 180) * cos(2 * pi * ((jd - 159.5) / 365.25))
  
  # Calculate the cosine of the solar zenith angle.
  coh <- sin(dec) * sin(latr) + cos(dec) * cos(latr) * cos(tt)
  
  # Calculate the solar zenith angle (z) by taking the arc cosine.
  # Ensure the argument to acos is within [-1, 1] to avoid NaNs.
  coh <- pmax(-1, pmin(1, coh))
  
  z <- acos(coh)
  
  return(z)
}

solzenmR <- function(jd, lt, lat, lon) {
  # This function calculates a matrix of solar zenith angles in degrees for multiple locations and times.
  # It is a direct translation of the C++ Rcpp code.
  
  # Check for matching dimensions
  if (length(lat) != length(lon)) {
    stop("Vectors 'lat' and 'lon' must have the same length.")
  }
  
  # Pre-allocate a list to store the results
  z_list <- vector("list", length(lat))
  
  # Loop through each location
  for (i in seq_along(lat)) {
    if (!is.na(lat[i])) {
      z_list[[i]] <- solzenvR(jd, lt, lat[i], lon[i])
    } else {
      z_list[[i]] <- rep(NA_real_, length(jd))
    }
  }
  
  # Combine the list of vectors into a matrix
  zd <- do.call(rbind, z_list)
  
  return(zd)
}

solzenvR <- function(jd, lt, lat, lon) {
  # This function calculates the solar zenith angle in degrees for each time point in a vector.
  # It assumes 'soltimeR' and 'solzenR' functions exist.
  
  # Use mapply for vectorized calculation
  # It iterates through each element of jd and lt
  st <- mapply(soltimeR, jd, lt, MoreArgs = list(lon = lon))
  z_rad <- mapply(solzenR, jd, st, MoreArgs = list(lat = lat, lon = lon))
  
  # Convert radians to degrees
  z_deg <- z_rad * 180 / base::pi
  
  return(z_deg)
}

tempintdayR <- function(tmn, tmnn, tmx, dl, stt, lat, lon, srte = 0.09) {
  # Calculates 24 hourly temperatures for a single day based on min/max temperatures.
  
  pi <- base::pi
  
  # Calculate predicted night fraction
  ngtp <- 0.04187957 * ((tmx - tmn) * (1 - dl / 24)) + 0.4372056
  ngtp <- pmax(0.01, pmin(0.99, ngtp)) # Clamp the value
  
  # Calculate sunrise time
  sr <- 12 - 0.5 * dl
  
  # Create a vector for 24 hours of a day
  thour <- numeric(24)
  
  for (lt in 0:23) {
    st <- lt + stt - sr # Solar time after sunrise
    
    if (st < 0) st <- st + 24
    if (st > 24) st <- st - 24
    
    # Handle polar day (24 hours of light)
    if (dl == 24) {
      thour[lt + 1] <- (tmx - tmn) * sin((pi * st) / 28) + tmn
      gr <- (tmnn - tmn) * (st / 24)
      thour[lt + 1] <- thour[lt + 1] + gr
    } else if (dl == 0) {
      # Handle polar night (0 hours of light)
      st <- lt + stt
      if (st < 0) st <- st + 24
      if (st > 24) st <- st - 24
      thour[lt + 1] <- (tmx - tmn) * sin((pi * st) / 24) + tmn
      gr <- (tmnn - tmn) * (st / 24)
      thour[lt + 1] <- thour[lt + 1] + gr
    } else {
      # Normal day (day and night cycle)
      k <- -(24 - dl) / log(srte / ngtp)
      ph <- -0.5 * dl * ((pi / (asin(ngtp) - pi)) + 1)
      rho <- dl + 2 * ph
      
      if (st > dl) {
        thour[lt + 1] <- ngtp * exp(-(st - dl) / k)
      } else {
        thour[lt + 1] <- sin((pi * st) / rho)
      }
      
      # Adjust by tmax and tmin
      thour[lt + 1] <- (tmx - tmn) * thour[lt + 1] + tmn
      
      # Apply gradient
      if (lt + stt > dl) {
        gr <- (tmnn - tmn) * (st - dl) / (24 - dl)
        thour[lt + 1] <- thour[lt + 1] + gr
      }
    }
  }
  
  # Adjust to ensure tmax and tmin match the calculated values
  ptmx <- max(thour, na.rm = TRUE)
  ptmn <- min(thour, na.rm = TRUE)
  
  b <- (ptmx - ptmn) / (tmx - tmn)
  a <- b * tmn - ptmn
  
  # Final adjustment
  for (lt in 0:23) {
    thour[lt + 1] <- (thour[lt + 1] + a) / b
  }
  
  return(thour)
}

.satvap <- function(tc) {
  e0<-(tc<0)*610.78/1000+(tc>=0)*611.2/1000
  L <- (tc<0)*2.834*10^6+(tc>=0)*((2.501*10^6)-(2340*tc))
  T0<-(tc<0)*273.15+(tc>=0)*273.15
  estl<-e0*exp((L/461.5)*(1/T0-1/(tc+273.15)))
  estl
}

.jday <- function(tme) {
  yr<-tme$year+1900
  mth<-tme$mon+1
  dd<-tme$mday+(tme$hour+(tme$min+tme$sec/60)/60)/24
  madj<-mth+(mth<3)*12
  yadj<-yr+(mth<3)*-1
  jd<-trunc(365.25*(yadj+4716))+trunc(30.6001*(madj+1))+dd-1524.5
  b<-(2-trunc(yadj/100)+trunc(trunc(yadj/100)/4))
  jd<-jd+(jd>2299160)*b
  jd
}

.soltime <- function(localtime, long, jd, merid = 0, dst = 0) {
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-localtime+(4*(long-merid)+eot)/60-dst
  st
}

.solalt <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  sa<-(180*atan(sh/sqrt(1-sh^2)))/pi
  sa
}

.solazi <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-159.5)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  hh<-0.5*pi-acos(sh)
  sazi<-cos(d)*sin(tt)/cos(hh)
  cazi<-(sin(lat*pi/180)*cos(d)*cos(tt)-cos(lat*pi/180)*sin(d))/
    sqrt((cos(d)*sin(tt))^2+(sin(lat*pi/180)*cos(d)*cos(tt)-cos(lat*pi/180)*sin(d))^2)
  sqt <- 1-sazi^2
  sqt[sqt<0]<-0
  solz<-180+(180*atan(sazi/sqrt(sqt)))/pi
  solz[cazi<0 & sazi<0]<-180-solz[cazi<0 & sazi<0]
  solz[cazi<0 & sazi>=0]<-540-solz[cazi<0 & sazi>=0]
  solz
}

.cfcelev<-function(csfc,dtmf,dtmc) {
  if (crs(dtmc) != crs(dtmf)) dtmc<-project(dtmc,crs(dtmf))
  if (crs(csfc) != crs(dtmf)) csfc<-project(csfc,crs(dtmf))
  # Calculate mean
  ca<-.is(csfc)
  cmean<-apply(ca,c(1,2),mean,na.rm=TRUE)
  # logit transform
  s0<-which(cmean==0)
  s1<-which(cmean==1)
  lc<-log(cmean/(1-cmean))
  lc[is.infinite(lc)]<-NA
  lc[s0]<-min(lc,na.rm=TRUE)
  lc[s1]<-max(lc,na.rm=TRUE)
  # resample
  lc<-.resample(.rast(lc,dtmc),dtmf)
  # Apply elevation adjustment to mean
  dc<-.is(.resample(dtmc,dtmf))
  ddif<-.is(dtmf)-dc
  ldif<- -9.938e-03-6.550e-04*ddif+4.966e-05*dc+1.358e-06*dc*ddif
  lcn<-.rast(.is(lc)+ldif,dtmf)
  lcn<-mask(lcn,dtmf)
  # Calculate difference in total between elevation adjusted and resampled / resampled
  mu<-.is(lcn-lc)
  mu<-.rta(mu,dim(csfc)[3])
  # Calculated individual logit transoformed radiation
  cf<-.resample(csfc,dtmf)
  cf<-.is(cf)
  s0<-which(cf==0)
  s1<-which(cf==1)
  cf<-log((cf)/(1-cf))
  cf[is.infinite(cf)]<-NA
  cf[s0]<-min(cf,na.rm=TRUE)
  cf[s1]<-max(cf,na.rm=TRUE)
  cf<-cf+mu
  # Back transform
  cf<-1/(1+exp(-cf))
  cf<-.rast(cf,dtmf)
  cf<-mask(cf,dtmf)
  return(cf)
}

.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}

.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}

.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}

.windhgt<-function (wspeed, zi, zo) {
  if (zo < 0.2 & zo > (5.42/67.8)){
    warning("Wind-height profile function performs poorly below 20 cm so output height converted to 20 cm")
    zo <- 0.2
  }
  return(wspeed * log(67.8 * zo - 5.42)/log(67.8 * zi - 5.42))
}

.is <- function(r) {
  if (class(r)[1] == "PackedSpatRaster") r<-terra::rast(r)
  if (class(r)[1] == "SpatRaster") {
    if (dim(r)[3] == 1) {
      m<-as.matrix(r,wide=TRUE)
    } else m<-as.array(r)
  } else {
    m<-r
  }
  return(m)
}

.tempelev <- function(tc, lrf, lrc, dtmf, dtmc = NA) {
  if (class(dtmc)[1] == "logical")  dtmc<-.resample(dtmf,tc)
  if (class(tc)[1] == "array") tc<-.rast(tc,dtmc)
  
  # Convert NA to elevation of 0 in dtmc
  dtmc<-ifel(is.na(dtmc),0,dtmc)
  
  # Lapse rate multipliers x elev
  #lrc<-as.array(lrc)*.rta(dtmc,n)
  #lrf<-as.array(lrf)*.rta(dtmf,n)
  lrc<-lrc*dtmc
  lrf<-lrf*dtmf
  
  # Sea-level temperature from coarse resolution lapse rate
  #stc<-.rast(.is(tc)+lrc,dtmc)
  stc<-tc+lrc
  if (crs(dtmc) != crs(dtmf)) stc<-project(stc,crs(dtmf))
  stc<-.resample(stc,dtmf)
  
  # Actual temperature from resampled sea temps and fine resolution lapse rate
  #tcf<-suppressWarnings(stc-.rast(lrf,dtmf))
  tcf<-suppressWarnings(stc-lrf)
  return(tcf)
}

.cad_multiplier<-function(dtmf, basins = NA, refhgt = 2){
  # Calculate elevation difference between basin height point and pixel
  if (class(basins) == "logical"){ basins<-basindelin(dtmf,refhgt)}
  b<-.is(basins)
  d<-.is(dtmf)
  u<-unique(as.vector(b))
  u<-u[is.na(u)==FALSE]
  bmx<-b*0
  for (i in 1:length(u)) {
    s<-which(b==u[i])
    mx<-max(d[s],na.rm=TRUE)
    bmx[s]<-mx
  }
  edif<-bmx-.is(dtmf)
  edif<-.rast(edif,dtmf)
  
  # Calculate lapse-rate multiplication factor (basin elev difference)
  mu<-edif*.cadpotential(dtmf,basins,refhgt)
  return(mu)
}
# dtm <- dtmf
basindelin <- function (dtm, boundary = 0) 
{
  dm <- dim(dtm)
  bsn <- .basindelin(dtm, boundary)
  return(bsn)
}
.edge<-function(v) {
  o<-0
  if (is.na(v[1]) == FALSE) {
    if (max(v,na.rm=TRUE) > v[1]) o<-1
  }
  o
}
.edgec<-function(v) {
  o<-v*0
  if (is.na(v[1]) == FALSE) {
    s<-which(v>v[1])
    o[s]<-1
  }
  o
}
.asign3<-function(bm2,bea,rw,cl) {
  b3<-bm2[rw:(rw+2),cl:(cl+2)]
  v<-bea[rw,cl,]
  if (is.na(v[2])==FALSE & v[2] > 0) b3[2,1]<-b3[2,2]
  if (is.na(v[3])==FALSE & v[3] > 0)  b3[2,3]<-b3[2,2]
  if (is.na(v[4])==FALSE & v[4] > 0)  b3[1,2]<-b3[2,2]
  if (is.na(v[5])==FALSE & v[5] > 0)  b3[1,1]<-b3[2,2]
  if (is.na(v[6])==FALSE & v[6] > 0)  b3[1,3]<-b3[2,2]
  if (is.na(v[7])==FALSE & v[7] > 0) b3[3,2]<-b3[2,2]
  if (is.na(v[8])==FALSE & v[8] > 0)  b3[3,1]<-b3[2,2]
  if (is.na(v[9])==FALSE & v[9] > 0)  b3[3,3]<-b3[2,2]
  b3
}
.basindelin_big<-function(dtm, boundary = 0, tilesize = 100, plotprogress = FALSE) {
  # chop into tiles
  e<-ext(dtm)
  reso<-res(dtm)
  xmxs<-as.numeric(ceiling((e$xmax-e$xmin)/reso[1]/tilesize))-1
  bma<-.docolumn(dtm,tilesize,boundary,0)
  for (x in 1:xmxs) {
    ta<-suppressWarnings(max(as.vector(bma),na.rm=T))
    if (is.infinite(ta)) ta<-0
    bo<-.docolumn(dtm,tilesize,boundary,x)+ta
    ed<-Sys.time()
    bma<-.basinmosaic(bma,bo)
    if (plotprogress) plot(bma,main=x)
  }
  m<-.is(bma)
  # renumber basins
  u<-unique(as.vector(m))
  u<-u[is.na(u) == FALSE]
  u<-u[order(u)]
  m<-array(renumberbasin(m,u),dim=dim(m))
  bout<-.rast(m,dtm)
  return(bout)
}

.basindelin<-function(dtm, boundary = 0) {
  # Delineate basins
  dm<-dim(dtm)
  me<-mean(as.vector(dtm),na.rm=TRUE)
  if (is.na(me) == FALSE) {
    bsn<-.basindelinR(dtm)
    # Merge basins if boundary > 0
    if (boundary > 0) {
      mx<-max(as.vector(bsn),na.rm=T)
      tst<-1
      while (tst == 1) {
        u<-unique(as.vector(bsn))
        u<-u[is.na(u) == F]
        if (length(u) > 1) {
          bsn<-.basinmerge(dtm,bsn,boundary)
          u<-unique(as.vector(bsn))
          u<-u[is.na(u) == F]
          if (length(u) > 1) {
            bsn<-.basinmerge(dtm,bsn,boundary)
            u<-unique(as.vector(bsn))
            u<-u[is.na(u) == F]
          }
          u<-unique(as.vector(bsn))
          u<-u[is.na(u) == F]
        }
        if (length(u) == 1) tst<-0
        mx2<-max(as.vector(bsn),na.rm=T)
        if (mx2 ==  mx) {
          tst<-0
        } else mx<-mx2
      } # end while
    } # end if boundary
  } else bsn<-dtm # end if boundary
  return(bsn)
}

.basindelinR<-function(dtm) {
  dm<-.is(dtm)
  dm[is.na(dm)]<-9999
  dm2<-array(9999,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # (2) create blank basin file
  bsn<-dm2*NA
  dun<-array(0,dim=dim(bsn))
  bsn<-basinR(dm2, bsn, dun)
  dd<-dim(bsn)
  bsn<-bsn[2:(dd[1]-1),2:(dd[2]-1)]
  if(class(bsn)[1]!='matrix') bsn<-matrix(bsn,ncol=ncol(dtm),nrow=nrow(dtm))
  r<-.rast(bsn,dtm)
  return(r)
}

.basinmosaic<-function(b1,b2) {
  e1<-ext(b1)
  e2<-ext(b2)
  reso<-res(b1)
  # *********** Do this if the tiles are vertically adjoined  *************** #
  if (abs(e1$ymax-e2$ymax) > reso[1]) {
    if (e2$ymax > e1$ymax) {  # b2 above b1
      m1<-.is(b1)
      m2<-.is(b2)
    } else {  # b1 above b2
      m1<-.is(b2)
      m2<-.is(b1)
    }
    for (itr in 1:3) {
      # merge based on top row of b1
      v1<-m1[1,] # top row of b1
      n<-dim(m2)[1] # mumber of rows
      v2<-m2[n,] # bottom row of b2
      # Create unique pairs matrix
      mup<-as.matrix(cbind(v1,v2))
      mup<-unique(mup)
      s<-which(is.na(mup[,1])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      s<-which(is.na(mup[,2])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      # Create vector of unique v1s
      u1<-unique(v1)
      u1<-u1[is.na(u1)==FALSE]
      u1<-u1[order(u1)]
      ras2<-list() # list of basins in m2 that should be re-asigned for each basin in u1
      ras1<-list() # list of basins in m1 that should be re-asigned for each basin in u1
      if (length(u1) > 0) {
        for (i in 1:length(u1)) {
          s<-which(mup[,1]==u1[i])
          u2<-mup[s,2] # list of basins in m2 that need reassigned
          u2<-u2[order(u2)]
          ras2[[i]]<-u2
          # list of basins in m1 that need reassinged
          s<-which(mup[,2]==u2[1])
          u1n<-mup[s,1] # list of basins in m2 that need reassigned
          if (length(u2) > 1) {
            for (j in 2:length(u2)) {
              s<-which(mup[,2]==u2[j])
              u1n<-c(u1n,mup[s,1]) # list of basins in m2 that need reassigned
            }
          }
          u1n<-unique(u1n)
          u1n<-u1n[u1n>u1[i]]
          ras1[[i]]<-u1n[order(u1n)]
          u2<-ras2[[i]]
          # Reassign basins in m2
          if (length(u2) > 0) for (j in 1:length(u2)) m2[m2==u2[j]]<-u1[i]
          u1n<-ras1[[i]]
          # Reassign basins in m1
          if (length(u1n) > 0) for (j in 1:length(u1n)) m1[m1==u1n[j]]<-u1[i]
        } # end for u1
      } # end if u1
    } # end iter
    # Convert back to SpatRasts
    if (e2$ymax > e1$ymax) {  # b2 above b1
      b1n<-.rast(m1,b1)
      b2n<-.rast(m2,b2)
    } else {  # b1 above b2
      b1n<-.rast(m2,b1)
      b2n<-.rast(m1,b2)
    }
  } else {# end do this if the tiles are vertically adjoined
    # *********** Do this if the tiles are horizontally adjoined  ************** #
    if (e2$xmax > e1$xmax) {  # b2 right of b1
      m1<-.is(b1)
      m2<-.is(b2)
    } else {  # b2 left of b1
      m1<-.is(b2)
      m2<-.is(b1)
    }
    for (itr in 1:3) {
      # merge based on right hand column of b1
      n<-dim(m1)[2]
      v1<-m1[,n] # right-hand column of b1
      v2<-m2[,1] # left-hand column of b2
      # Create unique pairs matrix
      mup<-as.matrix(cbind(v1,v2))
      mup<-unique(mup)
      s<-which(is.na(mup[,1])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      s<-which(is.na(mup[,2])==FALSE)
      mup<-mup[s,]
      if (class(mup)[1] != "matrix") mup<-t(as.matrix(mup))
      # Create vector of unique v1s
      u1<-unique(v1)
      u1<-u1[is.na(u1)==FALSE]
      u1<-u1[order(u1)]
      ras2<-list() # list of basins in m2 that should be re-asigned for each basin in u1
      ras1<-list() # list of basins in m1 that should be re-asigned for each basin in u1
      if (length(u1) > 0) {
        for (i in 1:length(u1)) {
          s<-which(mup[,1]==u1[i])
          u2<-mup[s,2] # list of basins in m2 that need reassigned
          u2<-u2[order(u2)]
          ras2[[i]]<-u2
          # list of basins in m1 that need reassinged
          s<-which(mup[,2]==u2[1])
          u1n<-mup[s,1] # list of basins in m2 that need reassigned
          if (length(u2) > 1) {
            for (j in 2:length(u2)) {
              s<-which(mup[,2]==u2[j])
              u1n<-c(u1n,mup[s,1]) # list of basins in m2 that need reassigned
            }
          }
          u1n<-unique(u1n)
          u1n<-u1n[u1n>u1[i]]
          ras1[[i]]<-u1n[order(u1n)]
          u2<-ras2[[i]]
          # Reassign basins in m2
          if (length(u2) > 0) for (j in 1:length(u2)) m2[m2==u2[j]]<-u1[i]
          u1n<-ras1[[i]]
          # Reassign basins in m1
          if (length(u1n) > 0) for (j in 1:length(u1n)) m1[m1==u1n[j]]<-u1[i]
        } # end for u1
      } # end if u1
    } # end iter
    # Convert back to SpatRasts
    if (e2$xmax > e1$xmax) {  # b2 above b1
      b1n<-.rast(m1,b1)
      b2n<-.rast(m2,b2)
    } else {  # b1 above b2
      b1n<-.rast(m2,b1)
      b2n<-.rast(m1,b2)
    }
  }
  # ********************************** Mosaic ******************************* #
  bout<-mosaic(b1n,b2n)
  return(bout)
}

.docolumn<-function(dtm,tilesize,boundary,x) {
  e<-ext(dtm)
  reso<-res(dtm)
  ymxs<-as.numeric(ceiling((e$ymax-e$ymin)/reso[2]/tilesize))-1
  xmn<-as.numeric(e$xmin)+reso[1]*tilesize*x
  xmx<-xmn+reso[1]*tilesize
  ymn<-as.numeric(e$ymin)+reso[2]*tilesize*0
  ymx<-ymn+reso[2]*tilesize
  if (xmx > e$xmax) xmx<-e$xmax
  if (ymx > e$ymax) ymx<-e$ymax
  ec<-ext(xmn,xmx,ymn,ymx)
  dc<-crop(dtm,ec)
  bma<-basindelin(dc,boundary)
  # delineate basins for columns
  for (y in 1:ymxs) {
    xmn<-as.numeric(e$xmin)+reso[1]*tilesize*x
    xmx<-xmn+reso[1]*tilesize
    ymn<-as.numeric(e$ymin)+reso[2]*tilesize*y
    ymx<-ymn+reso[2]*tilesize
    if (xmx > e$xmax) xmx<-e$xmax
    if (ymx > e$ymax) ymx<-e$ymax
    ec<-ext(xmn,xmx,ymn,ymx)
    dc<-crop(dtm,ec)
    ta<-suppressWarnings(max(as.vector(bma),na.rm=T))
    if (is.infinite(ta)) ta<-0
    bo<-basindelin(dc,boundary)+ta
    bma<-.basinmosaic(bma,bo)
  } # end y
  return(bma)
}

.cadpotential <- function(dtm, basins = NA, refhgt = 2) {
  if (class(basins) == "logical") basins<-basindelin(dtm,refhgt)
  fa<-flowacc(dtm, basins) - 1
  # Calculate basin size
  fre<-freq(basins)
  b<-.is(basins)
  bsize<-b*0
  for (i in 1:length(fre$value))  {
    s<-which(b==fre$value[i])
    bsize[s]<-fre$count[i]
  }
  cadfr<-.is(fa)/bsize
  cadfr<-.rast(cadfr,dtm)
  cadfr[cadfr>1]<-1
  return(cadfr)
}
.cad_conditions<-function(tc,u2,swrad,lwrad,dtmc,dtmf,refhgt = 2){
  # determine whether cold-air drainage conditions exist
  d<-0.65*0.12
  zm<-0.1*0.12
  uf<-(0.4*u2)/log((refhgt-d)/zm)
  H<-(swrad+lwrad-(5.67*10^-8*0.97*(tc+273.15)^4))*0.5
  st<- -(0.4*9.81*(refhgt-d)*H)/(1241*(tc+273.15)*uf^3)
  if(!inherits(st,"SpatRaster")) st<-.rast(st,dtmc)
  if (crs(st) != crs(dtmf)) st<-project(st,crs(dtmf))
  st<-.resample(st,dtmf)
  st<-ifel(st>1,1,st)
  st<-ifel(st<1,0,st)
  return(st)
}
.apply_cad<-function(lrf,mu,st){
  cad<-lrf*-mu
  ce<-cad*st
  return(ce)
}
.basinmerge<-function(dtm,bsn,boundary=0.25) {
  # Put buffer around basin and dtn
  bm<-.is(bsn)
  bm2<-array(NA,dim=c(dim(bm)[1]+2,dim(bm)[2]+2))
  bm2[2:(dim(bm)[1]+1),2:(dim(bm)[2]+1)]<-bm
  dm<-.is(dtm)
  dm2<-array(NA,dim=c(dim(dm)[1]+2,dim(dm)[2]+2))
  dm2[2:(dim(dm)[1]+1),2:(dim(dm)[2]+1)]<-dm
  # Create 3D array of  basin numbers  with adjoining cells
  bma<-array(NA,dim=c(dim(bm),9))
  bma[,,1]<-bm # rw, cl
  bma[,,2]<-bm2[2:(dim(bm)[1]+1),1:dim(bm)[2]] # rw, cl-1
  bma[,,3]<-bm2[2:(dim(bm)[1]+1),3:(dim(bm)[2]+2)] # rw, cl+1
  bma[,,4]<-bm2[1:dim(bm)[1],2:(dim(bm)[2]+1)] # rw-1, cl
  bma[,,5]<-bm2[1:dim(bm)[1],1:dim(bm)[2]] # rw-1, cl-1
  bma[,,6]<-bm2[1:dim(bm)[1],3:(dim(bm)[2]+2)] # rw-1, cl+1
  bma[,,7]<-bm2[3:(dim(bm)[1]+2),2:(dim(bm)[2]+1)] # rw+1, cl
  bma[,,8]<-bm2[3:(dim(bm)[1]+2),1:dim(bm)[2]] # rw+1, cl-1
  bma[,,9]<-bm2[3:(dim(bm)[1]+2),3:(dim(bm)[2]+2)] # rw+1, cl+1
  # Create 3D array of elevation differences with adjoining cells
  dma<-array(NA,dim=c(dim(dm),9))
  dma[,,1]<-dm # rw, cl
  dma[,,2]<-dm2[2:(dim(dm)[1]+1),1:dim(dm)[2]]-dm # rw, cl-1
  dma[,,3]<-dm2[2:(dim(dm)[1]+1),3:(dim(dm)[2]+2)]-dm  # rw, cl+1
  dma[,,4]<-dm2[1:dim(dm)[1],2:(dim(dm)[2]+1)]-dm  # rw-1, cl
  dma[,,5]<-dm2[1:dim(dm)[1],1:dim(dm)[2]]-dm  # rw-1, cl-1
  dma[,,6]<-dm2[1:dim(dm)[1],3:(dim(dm)[2]+2)]-dm  # rw-1, cl+1
  dma[,,7]<-dm2[3:(dim(dm)[1]+2),2:(dim(dm)[2]+1)]-dm  # rw+1, cl
  dma[,,8]<-dm2[3:(dim(dm)[1]+2),1:dim(dm)[2]]-dm  # rw+1, cl-1
  dma[,,9]<-dm2[3:(dim(dm)[1]+2),3:(dim(dm)[2]+2)]-dm  # rw+1, cl+1
  dma2<-dma*0
  dma2[abs(dma)<boundary]<-1
  bma<-bma*dma2
  bma[,,1]<-bm
  # identify edge and basin merge cells
  be<-apply(bma,c(1,2),.edge)
  bea<-aperm(apply(bma,c(1,2),.edgec),c(2,3,1))
  s<-which(be>0,arr.ind=TRUE)
  for (i in 1:dim(s)[1]) {
    rw<-as.numeric(s[i,1])
    cl<-as.numeric(s[i,2])
    b3<-.asign3(bm2,bea,rw,cl)
    bm2[rw:(rw+2),cl:(cl+2)]<-b3
  }
  # reassign basin number
  u<-unique(as.vector(bm2))
  u<-u[is.na(u)==FALSE]
  u<-u[order(u)]
  bm3<-bm2
  for (i in 1:length(u)) {
    s<-which(bm2==u[i])
    bm3[s]<-i
  }
  dd<-dim(bm3)
  bsn<-bm3[2:(dd[1]-1),2:(dd[2]-1)]
  # case where single column/row and bsn is a vector
  if(class(bsn)[1]!='matrix') bsn<-matrix(bsn,ncol=ncol(dtm),nrow=nrow(dtm))
  r<-.rast(bsn,dtm)
}

.rast <- function(m,tem) {
  # Throw warning if dim of m do not match dim of tem
  if(any(dim(m)[1:2]!=dim(tem)[1:2])){
    warning("In .rast dimensions of matrix/array do not match dimensions of rast template!!!")
    #print(dim(m))
    #print(dim(tem))
  }
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}

whichminR <- function(dm2, dun) {
  # Identifies the 1-based row and column index of the lowest value in dm2
  # where the corresponding value in dun (done) is 0 (or < 1).
  
  # Create a mask: TRUE where dun < 1
  mask <- dun < 1
  
  # Apply the mask to dm2, setting 'done' values to Inf so they're ignored by min()
  masked_dm2 <- dm2
  masked_dm2[!mask] <- Inf
  
  # Find the minimum value in the masked matrix
  d <- min(masked_dm2)
  
  # If the minimum is still Inf, it means all pixels are done
  if (is.infinite(d)) {
    return(c(-1, -1)) # Return the C++ sentinel value
  }
  
  # Find the 1-based index (row, column) of the first occurrence of the minimum value
  # which(matrix, arr.ind = TRUE) returns a matrix of indices
  idx <- which(masked_dm2 == d, arr.ind = TRUE)[1, ]
  
  # The C++ function returns 0-based indices, but standard R returns 1-based.
  # We convert to 0-based to match the C++ output structure.
  rw <- idx[1] - 1
  cl <- idx[2] - 1
  
  return(c(rw, cl))
}

whichmin2R <- function(dm2, b, dun, bn) {
  # Identifies the 1-based row and column index of the lowest value in dm2
  # where the pixel belongs to basin 'bn' and is 'undone' (dun == 0).
  
  # Create the combined mask: (basin is bn) AND (not done)
  mask <- (b == bn) & (dun == 0)
  
  # Apply the mask: set irrelevant values to Inf
  masked_dm2 <- dm2
  masked_dm2[!mask] <- Inf
  
  # Find the minimum value
  d <- min(masked_dm2)
  
  # Check if all pixels for this basin are done
  if (is.infinite(d)) {
    return(c(-1, -1))
  }
  
  # Find the 1-based index (row, column) of the first occurrence of the minimum value
  idx <- which(masked_dm2 == d, arr.ind = TRUE)[1, ]
  
  # Convert to 0-based indices to match C++ output structure
  rw <- idx[1] - 1
  cl <- idx[2] - 1
  
  return(c(rw, cl))
}

sel3R <- function(dm, s) {
  # Selects a 3x3 matrix/array surrounding the focal cell given by 0-based index s.
  # R matrices are 1-based, so the C++ indices (s[0], s[1]) correspond to R indices (s[0]+1, s[1]+1).
  
  # R indices for the focal cell
  r_row <- s[1] + 1
  r_col <- s[2] + 1
  
  # R index ranges for the 3x3 neighborhood
  row_range <- (r_row - 1):(r_row + 1)
  col_range <- (r_col - 1):(r_col + 1)
  
  # Extract the 3x3 neighborhood
  m3 <- dm[row_range, col_range]
  
  return(m3)
}
# Since the logic for NumericMatrix and IntegerMatrix is the same, one function suffices.
# The type of 'dm' (numeric or integer matrix) will determine the output type.

assignhigherR <- function(m3, b3, bn) {
  # Assigns grid cells in the 3x3 neighborhood (b3) to basin 'bn' if their height
  # (in m3) is >= the focal cell's height (m3[2,2]) and they aren't already assigned.
  # Assumes m3 and b3 are 3x3 matrices.
  
  focal_height <- m3[2, 2] # R uses 1-based indexing: (2, 2) is the center
  
  # Create a mask: (height >= focal height) AND (not already assigned, b3 != 0) AND (not NoData, m3 != 9999)
  # The original C++ has b3 != 0, implying 0 is the unassigned basin value.
  mask <- (m3 >= focal_height) & (b3 != 0) & (m3 != 9999)
  
  # Apply the assignment
  b3[mask] <- bn
  
  return(b3)
}

slotinR <- function(b, b3, s) {
  # Slots the 3x3 neighborhood matrix (b3) back into the full basin matrix (b)
  # based on the 0-based index s.
  
  # R indices for the top-left corner of the 3x3 neighborhood
  r_start_row <- s[1]
  r_start_col <- s[2]
  
  # R index ranges for the 3x3 neighborhood in the full matrix (1-based)
  row_range <- r_start_row:(r_start_row + 2)
  col_range <- r_start_col:(r_start_col + 2)
  
  # Slot the 3x3 matrix back into the full matrix
  b[row_range, col_range] <- b3
  
  return(b)
}

# Helper function to add 1-cell padding of NAs/0s
pad_matrix <- function(M, padding_value) {
  rows <- nrow(M)
  cols <- ncol(M)
  
  # Create a matrix with 2 extra rows and columns for the padding
  M_padded <- matrix(padding_value, nrow = rows + 2, ncol = cols + 2)
  
  # Slot the original matrix into the center
  M_padded[2:(rows + 1), 2:(cols + 1)] <- M
  
  return(M_padded)
}

basinR <- function(dm2, bsn, dun) {
  # NOTE: This fixed function REQUIRES the input matrices to be padded first
  # and the resulting basin matrix to be unpadded at the end.
  
  # --- PADDING ADJUSTMENT ---
  # Pad the input matrices to handle boundary conditions safely.
  # Use a large number (9999) for DM2 padding to exclude it from 'whichmin'.
  # Use 0 for basin (unassigned) and 1 for 'dun' (already done/ignored).
  
  # dm2: Padded with a large number (or Inf) so it's never the minimum.
  dm2_p <- pad_matrix(dm2, padding_value = 9999) 
  
  # bsn: Padded with 0 (unassigned basin ID).
  bsn_p <- pad_matrix(bsn, padding_value = 0)
  
  # dun: Padded with 1 (done/ignore) so the border is never considered.
  dun_p <- pad_matrix(dun, padding_value = 1) 
  
  # Store original dimensions for un-padding later
  orig_rows <- nrow(dm2)
  orig_cols <- ncol(dm2)
  
  # The rest of the logic follows the C++ code, but we adjust the loops to
  # iterate over the *original* matrix's content by starting the index s at 1,1
  # and ending at N,N, which will correspond to the padded matrix's 2 to N+1 range.
  
  bn <- 1
  tsta <- 1
  while (tsta == 1) {
    # Find the lowest undone pixel in the PADDED matrix.
    # Note: whichminR must be modified to accept the PADDED dimensions!
    s <- whichminR(dm2_p, dun_p) # returns 0-based index [rw, cl]
    
    if (s[2] > -1) {
      s_r <- s + 1 # Convert to 1-based index for R matrix access
      
      # Assign basin and mark as done in PADDED matrices
      bsn_p[s_r[1], s_r[2]] <- bn
      dun_p[s_r[1], s_r[2]] <- 1
      
      # Select 3x3 matrix around dm2_p and bsn_p
      m3 <- sel3R(dm2_p, s)
      b3 <- sel3R(bsn_p, s)
      
      # Assign higher/equal cells
      b3 <- assignhigherR(m3, b3, bn)
      
      # Slot the updated b3 back into the full basin matrix
      bsn_p <- slotinR(bsn_p, b3, s)
    } else {
      tsta <- 0
    }
    
    # Subsequent iteration
    tst <- 1
    while (tst == 1) {
      s <- whichmin2R(dm2_p, bsn_p, dun_p, bn) # lowest undone pixel in basin
      
      if (s[2] > -1) {
        s_r <- s + 1
        
        # Assign basin and mark as done in PADDED matrices
        bsn_p[s_r[1], s_r[2]] <- bn
        dun_p[s_r[1], s_r[2]] <- 1
        
        # Select 3x3 matrix
        m3 <- sel3R(dm2_p, s)
        b3 <- sel3R(bsn_p, s)
        
        # Assign higher/equal cells
        b3 <- assignhigherR(m3, b3, bn)
        
        # Slot the updated b3 back into the full basin matrix
        bsn_p <- slotinR(bsn_p, b3, s)
      } else {
        tst <- 0
      }
    }
    bn <- bn + 1
  }
  
  # --- UN-PADDING ---
  # Return only the original content of the basin matrix
  bsn_final <- bsn_p[2:(orig_rows + 1), 2:(orig_cols + 1)]
  
  return(bsn_final)
}

renumberbasinR <- function(m) {
  # Renumbers the basin IDs in the vector m sequentially.
  # The C++ code takes a vector 'u' of unique IDs, but in R we can derive this.
  
  # Find unique, non-zero basin IDs (assuming 0 is 'No Data' or 'Unassigned')
  unique_ids <- sort(unique(m[m != 0]))
  
  # Create a named vector (a lookup table) for renumbering
  # Old ID -> New ID (1, 2, 3, ...)
  renumber_map <- seq_along(unique_ids)
  names(renumber_map) <- unique_ids
  
  # Create the new vector
  m_new <- m
  
  # Apply the renumbering using the lookup table
  # Loop through the unique IDs and replace them
  for (i in seq_along(unique_ids)) {
    old_id <- unique_ids[i]
    new_id <- i
    m_new[m_new == old_id] <- new_id
  }
  
  return(m_new)
}

extract_timedim <- function(nc) {
  # Extract time dimension
  # Specifically: pull out the first dimension that has 'time' in its name
  return(nc$dim[grepl("time", names(nc$dim))][[1]])
}

.applynotna<-function(a,fun) {
  if(inherits(a,c('numeric','integer'))|length(dim(a))==1) a<-array(a,dim=c(1,1,length(a)))
  m<-matrix(a,ncol=dim(a)[1]*dim(a)[2],byrow=T)
  sel<-which(is.na(m[1,])==F)
  r<-apply(m[,sel,drop=FALSE],2,fun)
  n<-dim(r)[1]
  ao<-array(NA,dim=c(dim(a)[1:2],n))
  sel<-which(is.na(a[,,1:n])==F)
  ao[sel]<-aperm(r,c(2,1))
  ao
}

.satvap <- function(tc) {
  e0<-(tc<0)*610.78/1000+(tc>=0)*611.2/1000
  L <- (tc<0)*2.834*10^6+(tc>=0)*((2.501*10^6)-(2340*tc))
  T0<-(tc<0)*273.15+(tc>=0)*273.15
  estl<-e0*exp((L/461.5)*(1/T0-1/(tc+273.15)))
  estl
}

.hourtoday<-function(a,fun=mean) {
  # Convert input a to 3D array
  if(inherits(a,"SpatRaster")){
    tem<-a[[1]]
    a<-.is(a)
    toSpatRaster<-TRUE
  } else toSpatRaster<-FALSE
  if(inherits(a,c('integer','numeric')) | length(dim(a))<3){
    a<-array(a,dim=c(1,1,length(a)))
    toVector<-TRUE
  } else toVector<-FALSE
  
  .htd<-function(x) {
    y<-matrix(x,ncol=24,byrow=T)
    apply(y,1,fun,na.rm=T)
  }
  d<-.applynotna(a,fun=.htd)
  if (toVector) d<-as.vector(d)
  if (toSpatRaster) d<-.rast(d,tem)
  return(d)
}

.ehr<-function(a) {
  n<-dim(a)[1]*dim(a)[2]
  o1<-rep(c(1:n),24*dim(a)[3])
  o2<-rep(c(1:dim(a)[3]),each=24*n)-1
  o2<-o2*max(o1,na.rm=T)
  o<-o1+o2
  ah<-rep(a,24)
  ah<-ah[o]
  ah<-array(ah,dim=c(dim(a)[1:2],dim(a)[3]*24))
  ah
}

humfromdew <- function(tdew, tc, p) {
  pk <- p / 1000
  ea <- 0.6108 * exp(17.27 * tdew / (tdew + 237.3))
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  s <- 0.622 * e0 / pk
  hr <- (ea / e0) * 100
  hs <- (hr / 100) * s
  hs
}
