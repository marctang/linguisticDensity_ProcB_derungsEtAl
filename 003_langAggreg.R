##################################################
## Name: 003_langAggreg  
## Script purpose: Count languages and langauge families for FP and HG at the grid point locations
## Date: 2018
## Author: Curdin Derungs
##################################################

library(rgdal)
library(geosphere)
rm(list=ls())

##defining the three spatial resolutions for the grid points, created with 001_regularGrid.R
#the script has to be run for each resolution seperately (uncomment resp. resolution)
resolution<-1171
# resolution<-295
# resolution<-3267

##loading grid points
load(paste("output/001_regularGrid/randPts_",resolution,".Rdata",sep=""))

##loading langauage points
load("output/002_langPrep/allLang.Rdata")


##intersecting grid points with languages
#computing distances between all grid points
reg.dist<-spDists(reg.spdf,longlat = T)
#diagonal is NA
diag(reg.dist)<-NA
#retrieving distance between each grid point and its nearest neighbour
reg.dist.min<-apply(reg.dist,2,function(x){min(x,na.rm=T)})

#delete large distance matrix
rm(list="reg.dist")

#compute distance between grid points and all FP and HG languages
fp.reg.dist<-spDists(reg.spdf@coords,subset(allLang,hunter_gatherer=="not hg")@coords,longlat = T)
hg.reg.dist<-spDists(reg.spdf@coords,subset(allLang,hunter_gatherer=="hg")@coords,longlat = T)

#compute minimum distance for each language to its neares neighbour grid point
fp.min<-apply(fp.reg.dist,2,min)
hg.min<-apply(hg.reg.dist,2,min)

#iterate through all grid points and count nearest languages
ind<-rep(NA,nrow(allLang))
fp<-NULL
fpFam<-NULL
hg<-NULL
hgFam<-NULL

#iteration over grid points
for(i in 1:nrow(reg.spdf)){
  
  #get index of nearest languages given the distance
  ind[fp.reg.dist[i,]==fp.min]<-i
  
  #count languages that are neares to grid point
  fp<-c(fp,as.numeric(table(fp.reg.dist[i,]==fp.min)["TRUE"]))
  hg<-c(hg,as.numeric(table(hg.reg.dist[i,]==hg.min)["TRUE"]))
  
  #count families that are nearest to grid point
  fams<-length(levels(factor(subset(allLang,hunter_gatherer=="not hg")$glottolog.stock[fp.reg.dist[i,]==fp.min])))
  fpFam<-c(fpFam,fams)
  fams<-length(levels(factor(subset(allLang,hunter_gatherer=="hg")$glottolog.stock[hg.reg.dist[i,]==hg.min])))
  hgFam<-c(hgFam,fams)
}

#combine all information
ns<-data.frame(fp=fp, fpFam=fpFam, hg=hg, hgFam=hgFam)
ns$fp[is.na(ns$fp)]<-0
ns$fpFam[is.na(ns$fpFam)]<-0
ns$hg[is.na(ns$hg)]<-0
ns$hgFam[is.na(ns$hgFam)]<-0

#add language counts to grid points
reg.spdf$fp<-ns$fp
reg.spdf$fpFam<-ns$fpFam
reg.spdf$hg<-ns$hg
reg.spdf$hgFam<-ns$hgFam

#save grid points with language counts
save(reg.spdf,file=paste("output/003_langAggreg/randPtsLang_",resolution,".Rdata",sep=""))
