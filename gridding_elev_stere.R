# gridding functions

library('Rgrib2')
library('sp')
library('fastgrid')

gdomain1 <- function (gribhandle) {
  LonEps <- 10^(-5)
  gridtype <- Ginfo(gribhandle, StrPar = "gridType")$gridType
  earthshape <- Ginfo(gribhandle, IntPar = c("editionNumber", 
                                             "shapeOfTheEarth", "scaledValueOfRadiusOfSphericalEarth", 
                                             "scaleFactorOfRadiusOfSphericalEarth", "scaleFactorOfEarthMajorAxis", 
                                             "scaledValueOfEarthMajorAxis", "scaleFactorOfEarthMinorAxis", 
                                             "scaledValueOfEarthMinorAxis"))
  earthproj <- switch(as.character(earthshape$shapeOfTheEarth), 
                      `0` = list(R = 6367470), `1` = list(R = 10^-earthshape$scaleFactorOfRadiusOfSphericalEarth * 
                                                            earthshape$scaledValueOfRadiusOfSphericalEarth), 
                      `2` = list(a = 6378160, b = 6356775), `3` = list(a = 10^(3 - 
                                                                                 earthshape$scaleFactorOfEarthMajorAxis) * earthshape$scaledValueOfEarthMajorAxis, 
                                                                       b = 10^(3 - earthshape$scaleFactorOfEarthMinorAxis) * 
                                                                         earthshape$scaledValueOfEarthMinorAxis), `4` = list(ellps = "GRS80"), 
                      `5` = list(ellps = "WGS84"), `6` = list(R = 6371229), 
                      `7` = list(a = 10^(-earthshape$scaleFactorOfEarthMajorAxis) * 
                                   earthshape$scaledValueOfEarthMajorAxis, b = 10^(-earthshape$scaleFactorOfEarthMinorAxis) * 
                                   earthshape$scaledValueOfEarthMinorAxis), list(R = 6371200))
  if (earthshape$shapeOfTheEarth >= 8) 
    warning(paste("This earth shape is not yet fully implemented. Defaulting to sphere with radius", 
                  6371200))
  if (gridtype == "regular_ll") {
    info <- gridtype
    projection <- list(proj = "latlong")
    ggg <- Ginfo(gribhandle, IntPar = c("Nx", "Ny", "iScansNegatively", 
                                        "jScansPositively", "jPointsAreConsecutive", "alternativeRowScanning"), 
                 DblPar = c("latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees", 
                            "latitudeOfLastGridPointInDegrees", "longitudeOfLastGridPointInDegrees", 
                            "iDirectionIncrementInDegrees", "jDirectionIncrementInDegrees"))
    if (ggg$iScansNegatively) {
      Lon2 <- ggg$longitudeOfFirstGridPointInDegrees
      Lon1 <- ggg$longitudeOfLastGridPointInDegrees
    }
    else {
      Lon1 <- ggg$longitudeOfFirstGridPointInDegrees
      Lon2 <- ggg$longitudeOfLastGridPointInDegrees
    }
    Lon1 <- Lon1%%360
    if (Lon1 > 180) 
      Lon1 <- Lon1 - 360
    Lon2 <- Lon2%%360
    if (Lon2 > 180) 
      Lon2 <- Lon2 - 360
    if (Lon1 > Lon2) {
      if (Lon2 > 0) 
        Lon1 <- Lon1 - 360
      else Lon2 <- Lon2 + 360
      if (Lon1 >= 0) 
        lon0 <- 180
      else if (Lon2 <= 0) 
        lon0 <- -180
      else lon0 <- (Lon1 + Lon2)/2
    }
    else {
      lon0 <- 0
    }
    projection$lon0 <- lon0
    if (abs((Lon2 - Lon1)/(ggg$Nx - 1) - ggg$iDirectionIncrementInDegrees) > 
        LonEps) {
      warning(paste("Longitudes may be inconsistent: Lon1=", 
                    Lon1, "Lon2=", Lon2, "Nx=", ggg$Nx, "Dx=", ggg$iDirectionIncrementInDegrees))
      delx <- (Lon2 - Lon1)/(ggg$Nx - 1)
    }
    else {
      delx <- ggg$iDirectionIncrementInDegrees
    }
    if (ggg$jScansPositively) {
      Lat1 <- ggg$latitudeOfFirstGridPointInDegrees
      Lat2 <- ggg$latitudeOfLastGridPointInDegrees
    }
    else {
      Lat2 <- ggg$latitudeOfFirstGridPointInDegrees
      Lat1 <- ggg$latitudeOfLastGridPointInDegrees
    }
    if (Lat1 > Lat2) 
      warning(paste("Inconsistent Lat1=", Lat1, "Lat2=", 
                    Lat2))
    if (abs((Lat2 - Lat1)/(ggg$Ny - 1) - ggg$jDirectionIncrementInDegrees) > 
        LonEps) {
      warning(paste("Latitudes may be inconsistent: Lat1=", 
                    Lat1, "Lat2=", Lat2, "Ny=", ggg$Ny, "Dy=", ggg$jDirectionIncrementInDegrees))
      dely <- (Lat2 - Lat1)/(ggg$Ny - 1)
    }
    else {
      dely <- ggg$jDirectionIncrementInDegrees
    }
    SW <- c(Lon1, Lat1)
    NE <- c(Lon2, Lat2)
    result <- list(projection = projection, nx = ggg$Nx, 
                   ny = ggg$Ny, SW = SW, NE = NE, dx = delx, dy = dely)
    class(result) <- "geodomain"
    result
  }
  else if (gridtype == "lambert") {
    info <- gridtype
    ggg <- Ginfo(gribhandle, IntPar = c("Nx", "Ny", "iScansNegatively", 
                                        "jScansPositively", "jPointsAreConsecutive", "alternativeRowScanning"), 
                 DblPar = c("Latin1InDegrees", "Latin2InDegrees", 
                            "LoVInDegrees", "DxInMetres", "DyInMetres", "latitudeOfFirstGridPointInDegrees", 
                            "longitudeOfFirstGridPointInDegrees"))
    La1 <- ggg$latitudeOfFirstGridPointInDegrees
    Lo1 <- ggg$longitudeOfFirstGridPointInDegrees
    delx <- ggg$DxInMetres
    dely <- ggg$DyInMetres
    nx <- ggg$Nx
    ny <- ggg$Ny
    if (Lo1 > 180) 
      Lo1 <- Lo1 - 360
    rlat1 <- ggg$Latin1InDegrees
    rlat2 <- ggg$Latin2InDegrees
    rlon <- ggg$LoVInDegrees
    projection <- c(list(proj = "lcc", lon_0 = rlon, lat_1 = rlat1, 
                         lat_2 = rlat2), earthproj)
    SW <- c(Lo1, La1)
    if (!requireNamespace("meteogrid", quietly = TRUE)) {
      warning("meteogrid is not available! Some grid details can not be computed.")
      NE <- NULL
    }
    else {
      xy <- meteogrid::project(c(Lo1, La1), proj = projection, 
                               inv = FALSE)
      x0 <- xy$x[1]
      y0 <- xy$y[1]
      x1 <- x0 + (nx - 1) * delx
      y1 <- y0 + (ny - 1) * dely
      xy <- meteogrid::project(c(x1, y1), proj = projection, 
                               inv = TRUE)
      NE <- c(xy$x, xy$y)
    }
    result <- list(projection = projection, nx = nx, ny = ny, 
                   SW = SW, NE = NE, dx = delx, dy = dely)
    class(result) <- "geodomain"
    result
  }
  else if (gridtype == "polar_stereographic") {
    info <- gridtype
    ggg <- Ginfo(gribhandle, IntPar = c("Nx", "Ny", "iScansNegatively", 
                                        "jScansPositively", "jPointsAreConsecutive", "alternativeRowScanning"), 
                 #DblPar = c("LoVInDegrees", "latitudeOfFirstGridPointInDegrees", 
                 #            "longitudeOfFirstGridPointInDegrees", "DxInDegrees", 
                 #            "DyInDegrees"))
                 DblPar = c("latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees", 
                            "DxInMetres", "DyInMetres", "LaDInDegrees"))
    La1 <- ggg$latitudeOfFirstGridPointInDegrees
    Lo1 <- ggg$longitudeOfFirstGridPointInDegrees
    Lo1 <- Lo1%%360
    if (Lo1 > 180) 
      Lo1 <- Lo1 - 360
    nx <- ggg$Nx
    ny <- ggg$Ny
    rlon <- ggg$LoVInDegrees
    delx <- ggg$DxInMetres
    dely <- ggg$DyInMetres
    projection <- c(list(proj = "stere", lon_0 = 0, lat_0 = 90), 
                    earthproj)
    #projection <- c(list(proj = "merc", lat_ts = rlat), earthproj)
    SW <- c(Lo1, La1)
    if (!requireNamespace("meteogrid", quietly = TRUE)) {
      warning("meteogrid is not available! Some grid details can not be computed.")
      NE <- NULL
    }
    else {
      xy <- meteogrid::project(x = Lo1, y = La1, proj = projection, 
                               inv = FALSE)
      x0 <- xy$x[1]
      y0 <- xy$y[1]
      x1 <- x0 + (nx - 1) * delx
      y1 <- y0 + (ny - 1) * dely
      xy <- meteogrid::project(list(x = x1, y = y1), proj = projection, 
                               inv = TRUE)
      NE <- c(xy$x, xy$y)
    }
    result <- list(projection = projection, nx = nx, ny = ny, 
                   SW = SW, NE = NE, dx = delx, dy = dely)
    class(result) <- "geodomain"
    result
  }
  else if (gridtype == "mercator") {
    info <- gridtype
    ggg <- Ginfo(gribhandle, IntPar = c("Nx", "Ny", "iScansNegatively", 
                                        "jScansPositively", "jPointsAreConsecutive", "alternativeRowScanning"), 
                 DblPar = c("latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees", 
                            "DxInMetres", "DyInMetres", "LaDInDegrees"))
    La1 <- ggg$latitudeOfFirstGridPointInDegrees
    Lo1 <- ggg$longitudeOfFirstGridPointInDegrees
    if (Lo1 > 180) 
      Lo1 <- Lo1 - 360
    nx <- ggg$Nx
    ny <- ggg$Ny
    delx <- ggg$DxInMetres
    dely <- ggg$DyInMetres
    rlat <- ggg$LaDInDegrees
    projection <- c(list(proj = "merc", lat_ts = rlat), earthproj)
    SW <- c(Lo1, La1)
    if (!requireNamespace("meteogrid", quietly = TRUE)) {
      warning("meteogrid is not available! Some grid details can not be computed.")
      NE <- NULL
    }
    else {
      xy <- meteogrid::project(list(x = Lo1, y = La1), 
                               proj = projection, inv = FALSE)
      x0 <- xy$x[1]
      y0 <- xy$y[1]
      x1 <- x0 + (nx - 1) * delx
      y1 <- y0 + (ny - 1) * dely
      xy <- meteogrid::project(list(x = x1, y = y1), proj = projection, 
                               inv = TRUE)
      NE <- c(xy$x, xy$y)
    }
    result <- list(projection = projection, nx = nx, ny = ny, 
                   SW = SW, NE = NE, dx = delx, dy = dely)
    class(result) <- "geodomain"
    result
  }
  else if (gridtype == "rotated_ll") {
    info <- gridtype
    ggg <- Ginfo(gribhandle, IntPar = c("Nx", "Ny", "iScansNegatively", 
                                        "jScansPositively", "jPointsAreConsecutive", "alternativeRowScanning"), 
                 DblPar = c("latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees", 
                            "latitudeOfLastGridPointInDegrees", "longitudeOfLastGridPointInDegrees", 
                            "iDirectionIncrementInDegrees", "jDirectionIncrementInDegrees", 
                            "angleOfRotationInDegrees", "latitudeOfSouthernPoleInDegrees", 
                            "longitudeOfSouthernPoleInDegrees"))
    nx <- ggg$Nx
    ny <- ggg$Ny
    if (ggg$iScansNegatively) {
      Lon2 <- ggg$longitudeOfFirstGridPointInDegrees
      Lon1 <- ggg$longitudeOfLastGridPointInDegrees
    }
    else {
      Lon1 <- ggg$longitudeOfFirstGridPointInDegrees
      Lon2 <- ggg$longitudeOfLastGridPointInDegrees
    }
    if (ggg$jScansPositively) {
      Lat1 <- ggg$latitudeOfFirstGridPointInDegrees
      Lat2 <- ggg$latitudeOfLastGridPointInDegrees
    }
    else {
      Lat2 <- ggg$latitudeOfFirstGridPointInDegrees
      Lat1 <- ggg$latitudeOfLastGridPointInDegrees
    }
    if (Lat1 > Lat2) 
      warning(paste("Inconsistent Lat1=", Lat1, "Lat2=", 
                    Lat2))
    SW <- c(Lon1, Lat1)
    NE <- c(Lon2, Lat2)
    info <- "Rotated LatLon grid"
    SPlat <- ggg$latitudeOfSouthernPoleInDegrees
    SPlon <- ggg$longitudeOfSouthernPoleInDegrees
    SPangle <- ggg$angleOfRotationInDegrees
    if (SPangle != 0) 
      warning("Rotated LatLon with SPangle not supported yet.")
    projection <- list(proj = "ob_tran", o_proj = "latlong", 
                       o_lat_p = -SPlat, o_lon_p = 0, lon_0 = SPlon)
    if (!requireNamespace("meteogrid", quietly = TRUE)) {
      warning("meteogrid is not available! Some grid details can not be computed.")
      SW <- NULL
      NE <- NULL
    }
    else {
      RR <- meteogrid::project(list(x = c(Lon1, Lon2) * 
                                      pi/180, y = c(Lat1, Lat2) * pi/180), proj = projection, 
                               inv = TRUE)
      SW <- c(RR$x[1], RR$y[1])
      NE <- c(RR$x[2], RR$y[2])
    }
    result <- list(projection = projection, nx = nx, ny = ny, 
                   SW = SW, NE = NE, dx = ggg$iDirectionIncrementInDegrees, 
                   dy = ggg$jDirectionIncrementInDegrees)
    class(result) <- "geodomain"
    result
  }
  else if (gridtype == "reduced_gg") {
    info <- "Reduced gaussian grid (experimental!)"
    N <- Ginfo(gribhandle, IntPar = c("N"))$N
    Nggg <- paste0("N", N)
    RGtable <- eval(parse(text = Nggg))
    nlon <- RGtable$reduced
    result <- list(name = Nggg, nlon = nlon, latlist = RGtable$latitude)
  }
  else {
    info <- "Unimplemented grid"
    projection <- list(proj = "unknown")
    warning(paste("Unknown grid:", gridtype))
    ggg <- Ginfo(gribhandle, IntPar = c("Nx", "Ny", "iScansNegatively", 
                                        "jScansPositively", "jPointsAreConsecutive", "alternativeRowScanning"))
    nx <- ggg$Nx
    ny <- ggg$Ny
    x0 <- 0
    y0 <- 0
    delx <- 1/(nx - 1)
    dely <- 1/(ny - 1)
    x1 <- 1
    y1 <- 1
    asp <- 1
    return(list(info = info, projection = projection, grid = gridtype, 
                nx = nx, ny = ny, x0 = 0, y0 = 0, x1 = 1, y1 = 1, 
                delx = delx, dely = dely))
  }
}



# GRIB Lambert projection
crs.lambert <- CRS('+ellps=WGS84 +proj=lcc +lon_0=15.0 +lat_0=0.0 +x_0=0.0 +y_0=0 +lat_1=63.3 +lat_2=63.3 +no_defs')
crs.lonlat <- CRS("+init=epsg:4326")

crs.polster <- CRS('+proj=stere +lat_0=90 +lat_ts=60 +lon_0=20 +k=1 +x_0=0 +y_0=0 +a=6371220 +b=6371220 +units=m +no_defs')


# Transform one x y point from coordinate system p1 to p2
transpoint<- function(x,y,p1,p2) {
  xy <- SpatialPoints(matrix(c(x,y),nrow=1),proj4string = p1)
  coordinates(spTransform(xy,p2))
}

# Read valid times from a grib file and return as POSIXct
getfcdates <- function(f1,tz='UTC') {
  g <- Gopen(f1)
  # take 3 number time values into account
  a <- format(strptime(sprintf('%04d', g$validityTime), format='%H%M'), '%H%M') 
  validtime <- as.POSIXct(paste(g$validityDate,a),tz=tz,format="%Y%m%d %H%M")
  return(validtime)
}

# return the grib parameter name
getfcparam <- function(f1) {
  g <- Gopen(f1)
  param <- g$parameterNumber[1]
  return(param)
}

# read one field from MEPS grib file in Lambert projection
# returns output as SpatialGridDataFrame
readgrib <- function(f1,msg=1,gcrs=crs.polster,ocrs=crs.lonlat) {
  
  g <- Gopen(f1)
  gh <- Rgrib2::Ghandle(g,msg)
  gdat  <- Rgrib2::Gdec(gh)
  gdom <- gdomain1(gh)
  
  IntPar <- c("Nx", "Ny", "iScansNegatively", "jScansPositively",
              "jPointsAreConsecutive",
              "alternativeRowScanning","missingValue", "numberOfMissing")
  StrPar <- c("units")
  
  ginf <- Rgrib2::Ginfo(gh, IntPar = IntPar, StrPar = StrPar)
  
  dx = gdom$dx
  dy = gdom$dy
  nx = gdom$nx
  ny = gdom$ny
  lon0 <- gdom$SW[1]
  lat0 <- gdom$SW[2]
  if (lon0 > 180) lon0 <-lon0-360
  
  xy0 <- transpoint(lon0,lat0,ocrs,gcrs)
  x0<-xy0[1]
  y0<-xy0[2]
  
  #x1 <- x0+(nx-1)*dx
  #y1 <- y0+(ny-1)*dy
  #ll2 = transpoint(x1,y1,p1,p0)
  #x <- seq(x0,x1,len=nx)
  #y <- seq(y0,y1,len=ny)
  
  if (ginf$jScansPositively==1){ # check this!!!
    gdat <- gdat[,dim(gdat)[2]:1]
  }
  
  gt <- sp::GridTopology(cellcentre.offset = c(x0,y0),
                         cellsize = c(dx,dy),
                         cells.dim = c(nx,ny))
  
  out<-sp::SpatialGridDataFrame(gt, data.frame(VAR1=as.vector(gdat)),proj4string = gcrs)
  sp::coordnames(out) <- c('x','y')
  sp::gridded(out)<-TRUE
  sp::fullgrid(out) <- TRUE
  
  GhandleFree(gh)
  
  return(out)
}

# copy msg from gin to gout and replace data with newdata
#Gmod(gh,data=newdata,IntPar=list(centre=89,productDefinitionTemplateNumber=0))
savegrib <- function(gin,gout,msg=1,newdata=NULL,append=FALSE) {
  g <- Gopen(gin) # open the original grib file
  gh <- Ghandle(g,msg)
  if (!is.null(newdata)) {
    dim(newdata) <- c(g$Nx[msg],g$Ny[msg]) # reshape
    newdata <- newdata[,dim(newdata)[2]:1] # swap column order
    Gmod(gh,data=newdata,IntPar=list(centre=86,productDefinitionTemplateNumber=70,generatingProcessIdentifier=207)) # replace grib data with new
    # Modify grib meta-data generatingProcessIndentifier=207 is "MNWC Himan" producer
  }
  Gwrite(gh,filename = gout, append = append) # save to file
  GhandleFree(gh)
  invisible(NULL)
}

# gridding point observations to background grid using 'fastgrid' package
gridobs <- function(obs,grid, variable="VAR1", sigma = 5.0, nugget=0.2, clen=1000*10, lsm=NULL, melev=NULL) {
  
  covpars <- c(sigma,clen,nugget)
  trend_model <- as.formula(paste(variable,'~ -1'))
  coordnames(grid) <- c('longitude','latitude')
  var.pred <- fastkriege(trend_model = trend_model, obs, grid, cov.pars = covpars, 
                         lsm=lsm,lsmy=obs$lsm, alt=melev, alty=obs$elevation,
                         bg=grid, variable=variable, LapseRate = 0.0,
                         method='nearest')
  
  return(var.pred)
}

# quality check for corrected forecast 
qcheck <- function(fcvalue,param,grib_paramnb) {
  if(param==grib_paramnb$parnumb[2]) { # rh
    fcvalue[fcvalue>1]<-1
    fcvalue[fcvalue<0]<-0
  }
  if(param==grib_paramnb$parnumb[3]) { # ws
    fcvalue[fcvalue<0]<-0
  }
  return(fcvalue)
}

# quality check for error correction, limit the max difference allowed to BC
diff_qcheck <- function(diffvalue,param,grib_paramnb) {
  if(param==grib_paramnb$parnumb[2]) { # rh
    diffvalue[diffvalue>0.2]<-0.2
    diffvalue[diffvalue<(-0.2)]<-(-0.2)
  }
  # if(param==grib_paramnb$parnumb[3]) { # ws
  #  fcvalue[fcvalue<0]<-0
  # }
  return(diffvalue)
}

# obs processing for NetAtmo data
# obsNetAx <- obs_prepareNetA(obsNetA,t1,raster::extent(out),LSM,grid=var.pred,melev=MELEV)
# obsN <- obsNetA
# fc_time<-t1
# e <- raster::extent(out)
# lsm <- LSM
# grid=var.pred
# melev=MELEV
obs_prepareNetA <- function(obsN, fc_time, e=NULL, lsm=NULL, grid=NULL, melev=NULL) {
  # convert to MEPS coordinate system
  obsxPP <- spTransform(obsN,crs.lambert)
  obsxPP$VAR1 <- obsN$observation # name has to match that of the FC field
  #obsx <- obsx[obsx$time==fc_time,] # select ftime only not used for NetAtmo
  obsxPP <- obsxPP[complete.cases(obsxPP$VAR1),] # remove NaN values
  obsxPP <- obsxPP[complete.cases(obsxPP$elevation),]
  # Thinning NetAtmo observations by taking the mean of possible obs values in a raster grid
  R <- raster(grid) 
  r <- rasterize(obsxPP, R, c('VAR1','elevation'), fun=mean) 
  # Interpolate aggregated NetAtmo data back from raster to points
  NetAtmoPoints <- rasterToPoints(r,fun=NULL,spatial=TRUE)
  coordnames(NetAtmoPoints) <- c('longitude','latitude')
  #names(NetAtmoPoints) <- 'VAR1'
  # Interpolate points from SYNOP+model field and calculate the difference to NetAtmo data for QC
  MNWCpoints <- grid2points(var.pred,NetAtmoPoints,method='bilinear',variable='VAR1')[,1]
  # Calculate the difference between NetAtmo aggregated point values and model point values 
  # remove the values that differ more than |5| degrees  
  NetAtmoPoints$MNWCdiff <- NetAtmoPoints$VAR1 - MNWCpoints
  NetAtmoPoints <- NetAtmoPoints[(abs(NetAtmoPoints$MNWCdiff)<=5),]
  if (!is.null(e))
    NetAtmoPoints <- raster::crop(NetAtmoPoints,e)
  # add lsm to observation locations
  if (!is.null(lsm))
    NetAtmoPoints$lsm <- grid2points(lsm,NetAtmoPoints,method='nearest',variable='lsm')[,1]
  if (!is.null(melev))
    NetAtmoPoints$melev <- grid2points(melev,NetAtmoPoints,method='bilinear',variable='melev')[,1]
  return(NetAtmoPoints)
}

# obs processing
obs_prepare <- function(obs, fc_time, e=NULL, lsm=NULL, melev=NULL) {
  # convert to MEPS coordinate system
  obsx <- spTransform(obs,crs.polster)
  # obsx <- obsx[obsx$time==fc_time,] # select ftime only
  obsx$VAR1 <- obsx$observation # name has to match that of the FC field
  obsx <- obsx[complete.cases(obsx$VAR1),] # remove NaN values
  # crop needs the raster package
  # e <- raster::extent(out)
  if (!is.null(e))
    obsx <- raster::crop(obsx,e)
  # add lsm to observation locations
  if (!is.null(lsm))
    obsx$lsm <- grid2points(lsm,obsx,method='nearest',variable='lsm')
  if (!is.null(melev))
    obsx$melev <- grid2points(melev,obsx,method='bilinear',variable='melev')
  return(obsx)
}
