# Error correction of 1-5h forecasts for Metcoop nowcast (meps control) model for parameters t2m, rh2m and ws10m
# model data: 2t in Kelvins, Rh 0-1
# obs: 2t in Celsius, Rh 0-100
# return input fields, if no obs available 
#
# library(devtools)
# install_github("harphub/Rgrib2")
library('Rgrib2')
library('httr')
library('sp')
library('fastgrid')
library('rgeos')
#library('raster')
#library('outliers')
source('gridding_elev_stere.R') # subroutines
source('smartmet_obs_elev.R') # obs fetch from smartmet server
#source('file-access.R') # wrapper for reading/writing files from/to s3
library('data.table')

#args = commandArgs(trailingOnly=TRUE)
#args[1] = input-file
#args[2] = output-file
#args[3] = Netatmo TRUE/FALSE

#if (length(args) == 0) {
#  quit(status=1)
#}

#g_out <- strip_protocol(args[2])

#use_NetA <- args[3] # TRUE if NetAtmo data is used, FALSE if not

#g_in = '/data/hietal/pot.grib2'
#g_in <- '/data/hietal/snwc_dev.grib2'
#g_out = 'joku.grib2'

#the latest origtime used to create POT analysis field (just the metadata is copied)
#otime_old <- as.character(as.numeric(format(as.POSIXlt((Sys.time()-60*60), tz = "UTC"), "%Y%m%d%H"))) # time in utc -1h
otime <- as.character(as.numeric(format(as.POSIXlt(Sys.time(), tz = "UTC"), "%Y%m%d%H"))) # current uts in hours
g_out <- paste('/data/hietal/POT/analysis/',otime,'_0h.grib2',sep='')
#Ukkosdata 
txt <- paste("https://routines-data.lake.fmi.fi/smartmet-nwc/development/2.5/",otime,"/fcst_tstm.grib2",sep="")
destf <- paste('/data/hietal/POT/analysis/',otime,'.grib2',sep='')
download.file(url=txt,destfile=destf,quiet=TRUE)#, method="",extra = '--no-proxy --no-check-certificate',quiet=TRUE)

# t1<- as.character(as.numeric(format(as.POSIXlt((Sys.time()-0*60), tz = "UTC"), "%Y%m%d%H")))
# t2<- as.character(as.numeric(format(as.POSIXlt((Sys.time()+60*60), tz = "UTC"), "%Y%m%d%H")))
# t3<- as.character(as.numeric(format(as.POSIXlt((Sys.time()+120*60), tz = "UTC"), "%Y%m%d%H")))
# t4<- as.character(as.numeric(format(as.POSIXlt((Sys.time()+180*60), tz = "UTC"), "%Y%m%d%H")))
# #Current local hour 
# HOD <- strftime(as.POSIXct(Sys.time()), format="%H")
# tt <- c(paste("/data/hietal/POT/",HOD,"/",t1,"00.grib2",sep=""),paste("/data/hietal/POT/",HOD,"/",t2,"00.grib2",sep=""),paste("/data/hietal/POT/",HOD,"/",t3,"00.grib2",sep=""),paste("/data/hietal/POT/",HOD,"/",t4,"00.grib2",sep=""))
# first save the individual tpratet msg-files for motion vector using previous SNWC analysis time 
g_in <- destf

# data.frame to map grib input-file parameterNumber to corresponding smartmet server OBS name and error weights 
grib_paramnb <- data.frame('parname'=c('2t','rh','ws','tstm'),'parnumb'=c(0,192,1,2),'obs_parm'=c('TA_PT1M_AVG','RH_PT1M_AVG','WS_PT10M_AVG','flash'),'w1'=c(0.9,0.9,0.8,NaN),'w2'=c(0.9,0.8,0.7,NaN),'w3'=c(0.8,0.7,0.6,NaN),'w4'=c(0.7,0.7,0.6,NaN),stringsAsFactors=FALSE)

# get the parameter based on grib-input parameterNumber
param <- getfcparam(g_in) # t2m=0,rh=192,ws=1, tstm=2
# select the row with extra info needed for correct parameter
par.row <- which(grib_paramnb$parnumb==param)

# Use this if time is set based on sys.time 
# # define lates full hour in UTC (the time error is calculated from)
# #local.time <- Sys.time()
# local.time <- as.POSIXct('2020-07-28 13:00:00 UTC') # use this for test data! 
# attr(local.time, "tzone") <- "UTC"
# current.hour <- round(as.POSIXct(local.time, format="%H:%M:%S", tz="UTC"), units="hours")
# rm(local.time)
# fc_hours <- current.hour + c(0,1,2,3,4)*60*60

# for test data
# g_in <- 'smnwc-T-K.grib2' # temperature test dataset  
# g_in <- 'smnwc-RH-0TO1.grib2' # relative humidity test dataset
# g_in <- 'smnwc-FF-MS.grib2' # wind speed test dataset

# correlation length km
#clen <- 30 # limited for SYNOP from 50 --> 30km 
#clenx <- 10 # smaller correlation length used for NetAtmo data
#############
clenf <- 20 # correlation length for flash strikes
# land sea mask used
#LSM <- readRDS('MEPS_lsm.Rds') # only sea (+Vänern and Wettern)
#coordnames(LSM) <- c('longitude','latitude') # needed for gridding (should be fixed!)
# Model elevation
#MELEV <- readRDS('MEPS_alt.Rds')
#coordnames(MELEV) <- c('longitude','latitude')

#DEM <- raster("/data/statcal/projects/meps_calibration/data/lcc/HARP_dem_new_meps_100m.tif")
#DEM1 <- raster("/data/statcal/projects/meps_calibration/data/wgs84/HARP_dem_new_100m.tif")
# Forecast leadtimes 15min/1h (depends on the data) steps in Metcoop nowcast
mtimes <- getfcdates(g_in) # in grib files times <10 with 3 numbers, this function works also when 00:00 --> 0!

atime <- mtimes[1] # analysis times
fc_hours <- atime + c(1,2,3,4,5)*60*60 # 1h,2h,3h,4h,5h

# find indexes (=message numbers) in grib forecast times
msgs <- sapply(fc_hours, function(x) which(mtimes==x)) 

# start with 1h fc and do correction against observations
m <- msgs[1] 
t1 <- fc_hours[1] # 1h forecast of Metcoop NWC to used in error correction
out <- readgrib(g_in,msg=m)
coordnames(out) <- c('longitude','latitude') # needed for gridding

# define corresponding obs parameter name for smartmet server
obs_param <- grib_paramnb$obs_parm[par.row]

if (obs_param == 'flash') {
  t1 <- atime  
  startt <- t1-60*40 # flash data retrieved from 60min interval 40min to and 20min past the hour
  endt <- t1+60*20
  obs <- tryCatch(fread(paste("http://smartmet.fmi.fi/timeseries?producer=flash&tz=gmt&starttime=",startt,"&endtime=",endt,"&param=peak_current,flash_id,latitude,longitude,utctime&bbox=5,52,33,71&format=ascii",sep="")), error=function(e) NA)
  print(obs)
  print(startt)
  print(endt)
  print(length(obs))
  #obs <- fread(paste("http://smartmet.fmi.fi/timeseries?producer=flash&starttime=2021-04-22T01:00:00&endtime=2021-04-22T12:00:00&param=peak_current,flash_id,latitude,longitude,utctime&bbox=10,60,20,53&format=ascii",sep=""))
  if (nrow(obs)>1) {
    names(obs)<-c('name','fmisid','latitude','longitude','time')
    obs$observation <- 1
    obs <- obs_spatial(obs)
    obsx <- obs_prepare(obs,raster::extent(out))
    var.pred <- gridobs(obsx,out,clen=1000*clenf)# lsm=LSM$lsm)
    # diff = modified - original --> error correction 
    fcerr <- var.pred$diff  # error correction
    #fc.mod <- qcheck(var.pred$VAR1,param,grib_paramnb)
    # plot corrected & diff field 
    #var.pred$VAR1 <- var.pred$VAR1
    MOSplotting::MOS_plot_field(var.pred,layer = "VAR1", shapetrans = TRUE,cmin=0, cmax = 1,
                              main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_",m,".png",sep=""))
    
    #   MOSplotting::MOS_plot_field(var.pred,layer = "diff", shapetrans = TRUE,cmin=(0), cmax = 100,
    #                               stations = obsx,main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_",m,".png",sep=""))
    #
    #   MOSplotting::MOS_plot_field(out,layer = "VAR1", shapetrans = TRUE,cmin=0, cmax = 100,
    #                               stations = obsx,main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_",m,".png",sep=""))
    #   
    
    #savegrib(g_in,g_out,msg = m, newdata = fc.mod, append=FALSE)
    savegrib(g_in,g_out,msg = m, newdata = var.pred$VAR1, append=FALSE)
    } else if (file.exists(destf)) { # if obs == NULL remove downloaded grib-files  
        file.remove(destf) # delete the original grib file if no thunder
  #  savegrib(g_in,g_out,msg = m, newdata = out$VAR1, append=FALSE)
  }
  
  
}



# # obs from smartmet server
# # keyword for smartmet server obs retrieval: snwc includes the stations in ErrCorrectStations.csv
# obs <- readobs_all(t1,t1,obs_param,spatial = TRUE)
# 
# # if wind speed, then use potential wind speed opervations when available (for Finland)
# if(param==1)  {
#   # retrieve WSP obs data
#   obsPT <- readobs_all(t1,t1,"WSP_PT10M_AVG",spatial = TRUE)
#   # define the indexes for which there's WSP available 
#   ind_fmisid <- which(obs$fmisid %in% obsPT$fmisid)
#   # extra check in case there are WSPT values but no WS values (this has happened for station 101061) 
#   obsPT <- obsPT[which(obsPT$fmisid %in% obs$fmisid),]
#   # replace original obs values with WSP when available
#   obs$observation <- replace(obs$observation,ind_fmisid,obsPT$observation)
# }
# 
# if(!is.null(obs)) { # do error correction if there's obs available, if not return input fields 
#   # prepare obs
#   obsx <- obs_prepare(obs,t1,raster::extent(out), LSM, MELEV)
#   # obs for 2t in C --> K in model data
#   # obs for rh 0-100 --> 0-1 in model data
#   if (param==grib_paramnb$parnumb[1]){ # temperature from K --> C
#     obsx$VAR1 <- obsx$VAR1 + 273.15
#   }
#   if (param==grib_paramnb$parnumb[2]){ # RH from 0-1 --> 0-100
#     obsx$VAR1 <- obsx$VAR1/100
#   }
# 
#   # gridding
#   var.pred <- gridobs(obsx,out,clen=1000*clen, lsm=LSM$lsm, melev = MELEV$melev)
#   
#   # diff = modified - original --> error correction 
#   fcerr <- var.pred$diff  # error correction
# 
#   # quality control:
#   # RH<0-->0 & RH>1-->1
#   # ws<0-->0
#   fc.mod <- qcheck(var.pred$VAR1,param,grib_paramnb)
#   # limit the modifications BC can do
#   # RH diff field max +-20% 
#   var.pred$diff <- diff_qcheck(var.pred$diff,param,grib_paramnb)
# 
#   # plot corrected & diff field 
#   var.pred$VAR1 <- fc.mod
#   MOSplotting::MOS_plot_field(var.pred,layer = "VAR1", shapetrans = TRUE,cmin=min(out$VAR1), cmax = max(out$VAR1),
#                                 stations = obsx,main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_var",m,".png",sep=""))
#   
#   MOSplotting::MOS_plot_field(var.pred,layer = "diff", shapetrans = TRUE,cmin=(-5), cmax = 5,
#                               stations = obsx,main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("SYN_",m,".png",sep=""))
#   
#   # if temperature then add NetAtmo QC corrected obs to bias correction
#   # gridding point observations to background grid using 'fastgrid' package
#   gridobs1 <- function(obs,grid, variable="VAR1", sigma = 5.0, nugget=0.2, clen=1000*10, lsm=NULL) {
#     
#     covpars <- c(sigma,clen,nugget)
#     trend_model <- as.formula(paste(variable,'~ -1'))
#     coordnames(grid) <- c('longitude','latitude')
#     var.pred <- fastkriege(trend_model = trend_model, obs, grid, cov.pars = covpars, 
#                            lsm=lsm,lsmy=obs$lsm, alt=NULL, alty=NULL,
#                            bg=grid, variable=variable, LapseRate = 0.0,
#                            method='nearest')
#     
#     return(var.pred)
#   }
#   
#   if (use_NetA == TRUE & param==grib_paramnb$parnumb[1]){ 
#     var.pred$VAR1 <- fc.mod 
#     alku <- t1-(10*60) # NetAtmo data is taken from 10 min time interval xx:50-t1
#     obsNetA <- readobs_all(alku,t1,parname = 'NetAtmo',spatial = TRUE) # fetch NetAtm QC corrected obs
#     if(!is.null(obsNetA)) { # do error correction if there's NetAtmo obs available, if not return SYNOP corrected fields  
#       names(obsNetA)<-c('name','fmisid','time','elevation','observation')
#       obsNetA$observation <- obsNetA$observation + 273.15
#       # prepare obs
#       obsNetAx <- obs_prepareNetA(obsNetA,t1,raster::extent(out),LSM,grid=var.pred,melev=NULL)
#       # interpolate NetAtmo elevation from DEM 100m
#       obsNetAx$elevation <- extract(DEM,obsNetAx)
#       #tes <- extract(DEM1,obsNetAx)
#       # coordnames(obsNetAx) <- c('longitude','latitude')
#       var.predX <- gridobs(obsNetAx,var.pred,clen=1000*clenx,lsm=LSM$lsm, melev = MELEV$melev)
#       a <- Sys.time()
#       var.predX1 <- gridobs1(obsNetAx,var.pred,clen=1000*clenx,lsm=LSM$lsm)
#       b <- Sys.time()-a
#       # diff = modified - original --> error correction 
#       # error correction "fcerr" is now the SYNOP+NETATMO corrected field minus original model data   
#       fcerr <- var.predX$VAR1 - out$VAR1 # var.predX$diff  # error correction
#       fc.mod <- qcheck(var.predX$VAR1,param,grib_paramnb)
#       
#   MOSplotting::MOS_plot_field(var.predX,layer = "diff", shapetrans = TRUE,cmin=(-5), cmax = 5,
#                           main=paste(fc_hours[1], 'clen =',clenx,'km'),zoom=c(-1000,300000,8800000,9200000),pngfile=paste("mod_NA",m,".png",sep=""))
#   MOSplotting::MOS_plot_field(var.predX1,layer = "diff", shapetrans = TRUE,cmin=(-5), cmax = 5,
#                               main=paste(fc_hours[1], 'clen =',clenx,'km'),zoom=c(-1000,300000,8800000,9200000),pngfile=paste("mod_NA",m,".png",sep=""))
#   # MOSplotting::MOS_plot_field(var.predX,layer = "VAR1", shapetrans = TRUE,cmin=(250), cmax = 290,
#   #                            stations = obsNetAx,main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_VAR_NA",m,".png",sep=""))
#   # MOSplotting::MOS_plot_field(var.pred,layer = "VAR1", shapetrans = TRUE,cmin=(250), cmax = 290,
#   #                            main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_",m,".png",sep=""))
#   # MOSplotting::MOS_plot_field(out,layer = "VAR1", shapetrans = TRUE,cmin=(250), cmax = 290,
#   #                                main=paste(fc_hours[1], 'clen =',clen,'km'),pngfile=paste("mod_",m,".png",sep=""))
#       
#     }
#   }
#   savegrib(g_in,g_out,msg = m, newdata = fc.mod, append=FALSE)
# 
#   # now add the 1h correction to 2h,3h,4h and 5h
#   wcorr <- as.vector(t(grib_paramnb[par.row,c('w1','w2','w3','w4')])) # choose the relevant error weights 
#   for (m in 1:(length(msgs)-1)) {
#     mm <- msgs[2:length(msgs)][m] # index of grib message
#     out2 <- readgrib(g_in,msg=mm) # raw model data 
#     tmp_fc.mod2 <- out2$VAR1 + wcorr[m] * fcerr # error correction using varying weights for different forecast times
#     fc.mod2 <- qcheck(tmp_fc.mod2,param,grib_paramnb)
#     savegrib(g_in,g_out,msg = mm, newdata = fc.mod2, append=TRUE)
#   }
# } else { # if obs == NULL save just the input fields 
#   for (m in 1:(length(msgs))) {
#     mm <- msgs[m] # index of grib message
#     out2 <- readgrib(g_in,msg=mm) # raw model data 
#     fc.mod2 <- out2$VAR1
#     #fc.mod2 <- qcheck(tmp_fc.mod2,param,grib_paramnb)
#     savegrib(g_in,g_out,msg = mm, newdata = fc.mod2, append=TRUE)
#   }
# }
# 
# write_s3(args[2])


# extra plotting
# out2$mod <- fc.mod2 # corrected forecast
# out2$diff <- out2$mod - out2$VAR1 # difference to original
# 
# MOSplotting::MOS_plot_field(out2,layer = "mod", shapetrans = TRUE,cmin=min(out2$VAR1), cmax = max(out2$VAR1),
#                            stations = obsx,main=paste(fc_hours[2], 'clen =',clen,'km'),pngfile="testi.png")
# 
# MOSplotting::MOS_plot_field(out2,layer = "diff", shapetrans = TRUE,cmin=min(out2$diff), cmax = max(out2$diff),
#                            stations = obsx,main=paste(fc_hours[2], 'clen =',clen,'km'))
