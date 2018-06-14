##################################################
## Name: 001_regularGrid  
## Script purpose: Creating an even distribution of grid points for different spatial resolutions
## Date: 2018
## Author: Curdin Derungs
##################################################

library(geosphere)
library(rgdal)
library(plyr)
rm(list=ls())

##three spatial resolutions are defined as numbers of regions used in the following function
#the script has the be run for each resolution seperately
res<-30 #appr. 1000 grid points
# res<-15 #300
# res<-50 #3000

##rading continental shapes
#even grid points are only distributed over continental land masses
cont <- readOGR("input", "continents_simple")
proj4string(cont)<-CRS("+proj=longlat +datum=WGS84")


##creating regular points using the regularCoordinates() function from the geosphere package
#the function solves the Thompson Problem
reg.coords<-regularCoordinates(res)

##converting coordinates to spatial points
reg.spdf<-SpatialPointsDataFrame(SpatialPoints(reg.coords),data.frame(id=1:nrow(reg.coords)))

#defining the spatial projection
proj4string(reg.spdf)<-CRS("+proj=longlat +datum=WGS84")

##intersecting evenly distributed points with the continental polygons
ov<-over(reg.spdf,cont)
#filtering for points on continents
reg.spdf<-reg.spdf[!is.na(ov$OBJECTID),]

#save the grids
save(reg.spdf,file=paste("output/001_regularGrid/randPts_",nrow(reg.spdf),".Rdata",sep=""))
