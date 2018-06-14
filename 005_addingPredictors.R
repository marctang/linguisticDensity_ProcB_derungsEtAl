##################################################
## Name: 005_addingPredictors  
## Script purpose: Associationg grid points with environmental information, retrieved for different spatial resolutions.
## Date: 2018
## Author: Curdin Derungs
##################################################

library(rgdal)
library(raster)
library(maptools)
rm(list=ls())

##defining all resolutions
resolutions<-c(1171,295,3267)

##path setting for raster data
rasterOptions(tmpdir="input/environmentalVar/")

##read bio geotiff (www.worldclim.org)
#data from 6000bc
per_cv_6bc <- raster(readGDAL("input/environmentalVar/bcmidbi15.tif"),layer="band1")
per_wetqu_6bc<-raster(readGDAL("input/environmentalVar/bcmidbi16.tif"),layer="band1")
temp_mean_6bc<-raster(readGDAL("input/environmentalVar/bcmidbi1.tif"),layer="band1")
temp_warmqu_6bc<-raster(readGDAL("input/environmentalVar/bcmidbi8.tif"),layer="band1")

#define spatial projection
proj4string(per_cv_6bc)<-CRS("+proj=longlat +datum=WGS84")
proj4string(per_wetqu_6bc)<-CRS("+proj=longlat +datum=WGS84")
proj4string(temp_mean_6bc)<-CRS("+proj=longlat +datum=WGS84")
proj4string(temp_warmqu_6bc)<-CRS("+proj=longlat +datum=WGS84")

#data from 2000ad
per_cv_2ad <- raster(readGDAL("input/environmentalVar/bio_2000_15.tif"),layer="band1")
per_wetqu_2ad<-raster(readGDAL("input/environmentalVar/bio_2000_16.tif"),layer="band1")
temp_mean_2ad<-raster(readGDAL("input/environmentalVar/bio_2000_1.tif"),layer="band1")
temp_warmqu_2ad<-raster(readGDAL("input/environmentalVar/bio_2000_8.tif"),layer="band1")

#define spatial projection
proj4string(per_cv_2ad)<-CRS("+proj=longlat +datum=WGS84")
proj4string(per_wetqu_2ad)<-CRS("+proj=longlat +datum=WGS84")
proj4string(temp_mean_2ad)<-CRS("+proj=longlat +datum=WGS84")
proj4string(temp_warmqu_2ad)<-CRS("+proj=longlat +datum=WGS84")

##read landcover data (source: www.landcover.org)
#read data from 6000bc and 2000ad
crop_2ad <- raster(readAsciiGrid("input/environmentalVar/crop2000AD.asc"))
crop_6bc <- raster(readAsciiGrid("input/environmentalVar/crop6000BC.asc"))
gras_2ad <- raster(readAsciiGrid("input/environmentalVar/gras2000AD.asc"))
gras_6bc <- raster(readAsciiGrid("input/environmentalVar/pasture6000BC.asc"))

#define spatial projection
proj4string(crop_2ad)<-CRS("+proj=longlat +datum=WGS84")
proj4string(crop_6bc)<-CRS("+proj=longlat +datum=WGS84")
proj4string(gras_2ad)<-CRS("+proj=longlat +datum=WGS84")
proj4string(gras_6bc)<-CRS("+proj=longlat +datum=WGS84")

##read population data (source: http://themasites.pbl.nl/tridion/en/themasites/hyde/index.html)
#for 2000ad
pop_2ad <- raster(readAsciiGrid("input/environmentalVar/popd_2000AD.asc"))
proj4string(pop_2ad)<-CRS("+proj=longlat +datum=WGS84")

#for 6000bc
pop_6bc <- raster(readAsciiGrid("input/environmentalVar/popd_6000BC.asc"))
proj4string(pop_6bc)<-CRS("+proj=longlat +datum=WGS84")

##read a series of climatic data (source: http://www.worldclim.org/)
load("input/environmentalVar/warmMonthsRaster_2000.Rdata")
n_warm_month_2ad<-temps
proj4string(n_warm_month_2ad)<-CRS("+proj=longlat +datum=WGS84")

load("input/environmentalVar/warmMonthsRaster_6000bc.Rdata")
n_warm_month_6bc<-temps
proj4string(n_warm_month_6bc)<-CRS("+proj=longlat +datum=WGS84")
rm(list="temps")

##read distances to oceans and rivers
#both data sets are pre-processed with GIS-Software
load("input/environmentalVar/distOcean.Rdata")
load("input/environmentalVar/distRiver.Rdata")


##read digital elevation model (source: https://www.ngdc.noaa.gov/mgg/global/global.html)
dhm <- raster(read.asciigrid("input/environmentalVar/etopo_lowres.asc"))
proj4string(dhm)<-CRS("+proj=longlat +datum=WGS84")

#below sea level = 0
dhm@data@values[dhm@data@values<0]<-0 #no below sealevel
dhm@data@values[is.na(dhm@data@values)]<-0 #no below sealevel

#iterating through resolutions
for(i in 1:length(resolutions)){
  
  resolution<-resolutions[i]
  print(resolution)
  
  ##load the output from the previous script (004_random_NS.R)
  load(paste("output/004_random_NS/randPtsLangPoiss_",resolution,".Rdata",sep=""))
  
  ##intersect each grid point with all environmental variables
  #meso scale
  reg.spdf$per_cv_2ad<-extract(per_cv_2ad,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_mean_2ad<-extract(temp_mean_2ad,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$per_wetqu_2ad<-extract(per_wetqu_2ad,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_warmqu_2ad<-extract(temp_warmqu_2ad,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$n_warm_month_2ad<-extract(n_warm_month_2ad,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$gras_2ad<-extract(gras_2ad,reg.spdf,buffer=10000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$crop_2ad<-extract(crop_2ad,reg.spdf,buffer=10000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$pop_2ad<-extract(pop_2ad,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$per_cv_6bc<-extract(per_cv_6bc,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_mean_6bc<-extract(temp_mean_6bc,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$per_wetqu_6bc<-extract(per_wetqu_6bc,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_warmqu_6bc<-extract(temp_warmqu_6bc,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$n_warm_month_6bc<-extract(n_warm_month_6bc,reg.spdf,buffer=10000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$gras_6bc<-extract(gras_6bc,reg.spdf,buffer=10000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$crop_6bc<-extract(crop_6bc,reg.spdf,buffer=10000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$pop_6bc<-extract(pop_6bc,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$dist_river<-extract(distRiver,reg.spdf,buffer=10000,fun=function(x){min(x,na.rm=T)})
  reg.spdf$dist_ocean<-extract(distOcean,reg.spdf,buffer=10000,fun=function(x){min(x,na.rm=T)})
  reg.spdf$elev<-extract(dhm,reg.spdf,buffer=10000,fun=function(x){max(x,na.rm=T)})
  reg.spdf$elevSD<-extract(dhm,reg.spdf,buffer=100000,fun=function(x){sd(x,na.rm=T)})
  
  #micro scale
  reg.spdf$per_cv_2ad_0<-extract(per_cv_2ad,reg.spdf)
  reg.spdf$temp_mean_2ad_0<-extract(temp_mean_2ad,reg.spdf)
  reg.spdf$per_wetqu_2ad_0<-extract(per_wetqu_2ad,reg.spdf)
  reg.spdf$temp_warmqu_2ad_0<-extract(temp_warmqu_2ad,reg.spdf)
  reg.spdf$n_warm_month_2ad_0<-extract(n_warm_month_2ad,reg.spdf)
  reg.spdf$gras_2ad_0<-extract(gras_2ad,reg.spdf)
  reg.spdf$crop_2ad_0<-extract(crop_2ad,reg.spdf)
  reg.spdf$pop_2ad_0<-extract(pop_2ad,reg.spdf)
  reg.spdf$per_cv_6bc_0<-extract(per_cv_6bc,reg.spdf)
  reg.spdf$temp_mean_6bc_0<-extract(temp_mean_6bc,reg.spdf)
  reg.spdf$per_wetqu_6bc_0<-extract(per_wetqu_6bc,reg.spdf)
  reg.spdf$temp_warmqu_6bc_0<-extract(temp_warmqu_6bc,reg.spdf)
  reg.spdf$n_warm_month_6bc_0<-extract(n_warm_month_6bc,reg.spdf)
  reg.spdf$gras_6bc_0<-extract(gras_6bc,reg.spdf)
  reg.spdf$crop_6bc_0<-extract(crop_6bc,reg.spdf)
  reg.spdf$pop_6bc_0<-extract(pop_6bc,reg.spdf)
  reg.spdf$dist_river_0<-extract(distRiver,reg.spdf)
  reg.spdf$dist_ocean_0<-extract(distOcean,reg.spdf)
  reg.spdf$elev_0<-extract(dhm,reg.spdf)
  reg.spdf$elevSD_0<-extract(dhm,reg.spdf)
  
  #macro scale
  reg.spdf$per_cv_2ad_02<-extract(per_cv_2ad,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_mean_2ad_02<-extract(temp_mean_2ad,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$per_wetqu_2ad_02<-extract(per_wetqu_2ad,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_warmqu_2ad_02<-extract(temp_warmqu_2ad,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$n_warm_month_2ad_02<-extract(n_warm_month_2ad,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$gras_2ad_02<-extract(gras_2ad,reg.spdf,buffer=100000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$crop_2ad_02<-extract(crop_2ad,reg.spdf,buffer=100000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$pop_2ad_02<-extract(pop_2ad,reg.spdf,buffer=1000000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$per_cv_6bc_02<-extract(per_cv_6bc,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_mean_6bc_02<-extract(temp_mean_6bc,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$per_wetqu_6bc_02<-extract(per_wetqu_6bc,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$temp_warmqu_6bc_02<-extract(temp_warmqu_6bc,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$n_warm_month_6bc_02<-extract(n_warm_month_6bc,reg.spdf,buffer=100000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$gras_6bc_02<-extract(gras_6bc,reg.spdf,buffer=100000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$crop_6bc_02<-extract(crop_6bc,reg.spdf,buffer=100000,fun=function(x){modal(x,na.rm=T)})
  reg.spdf$pop_6bc_02<-extract(pop_6bc,reg.spdf,buffer=1000000,fun=function(x){median(x,na.rm=T)})
  reg.spdf$dist_river_02<-extract(distRiver,reg.spdf,buffer=100000,fun=function(x){min(x,na.rm=T)})
  reg.spdf$dist_ocean_02<-extract(distOcean,reg.spdf,buffer=100000,fun=function(x){min(x,na.rm=T)})
  reg.spdf$elev_02<-extract(dhm,reg.spdf,buffer=50000,fun=function(x){max(x,na.rm=T)})
  reg.spdf$elevSD_02<-extract(dhm,reg.spdf,buffer=500000,fun=function(x){sd(x,na.rm=T)})

  #save the grid points with the environmental information for all resolutions
  save(reg.spdf,file=paste("output/005_addingPredictors/randPtsLangPoissEnv_",resolution,".Rdata",sep=""))
}

