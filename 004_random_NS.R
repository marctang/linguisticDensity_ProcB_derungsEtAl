##################################################
## Name: 004_random_NS  
## Script purpose: Simulating random (i.e. expected) language counts for each grid point
## Date: 2018
## Author: Curdin Derungs
##################################################

library(rgdal)
library(geosphere)
rm(list=ls())

#automatic iteration over all grid resolutions
resolutions<-c(1171,295,3267)

for(i in 1:length(resolutions)){
  
  resolution<-resolutions[i]
  print(resolution)
  
  ##loading grid points
  load(paste("output/003_langAggreg/randPtsLang_",resolution,".Rdata",sep=""))
  
  ##loading language points
  load("output/002_langPrep/allLang.Rdata")
  
  ##reading a large dataset with global random points 
  #points were created with randomCoordinates() from geosphere, again incorporating the Thompson Problem
  #ArcMap was used to preprocess the data and only keeping continental points from some 1M random global points
  #This step could be automized in R but requires long computation (no spatial index)
  randCoords.cont<-readOGR(dsn = "input",layer =  "randPointsCont")
  randCoords.cont.df<-cbind(randCoords.cont@data,randCoords.cont@coords)
  
  #distance between grid points points (similar to previous script)
  reg.dist<-spDists(reg.spdf,longlat = T)
  diag(reg.dist)<-NA
  regMin<-apply(reg.dist,2,function(x){min(x,na.rm=T)})
  rm(list="reg.dist")
  
  #computing numbers or FP and HG and languages and language families, required for the previous computations
  #numbers of language LD
  nFP<-nrow(allLang[allLang$hunter_gatherer=="not hg",])
  nHG<-nrow(allLang[allLang$hunter_gatherer=="hg",])
  
  #numbers of family LD
  fpFam<-as.character(allLang[allLang$hunter_gatherer=="not hg",]$glottolog.stock)
  fpFamTab<-as.data.frame(table(fpFam))
  nFPFam<-length(unique(fpFam))
  hgFam<-as.character(allLang[allLang$hunter_gatherer=="hg",]$glottolog.stock)
  hgFamTab<-as.data.frame(table(hgFam))
  nHGFam<-length(unique(hgFam))
  
  #first iteration for FP
  fp.all<-numeric(nrow(reg.spdf))
  fpFam.all<-numeric(nrow(reg.spdf))
  
  ##creating 500 random distributions
  for(j in 1:500){
    print(j)
    
    ##taking a random sample from the random points of the size of FP languages
    randCoords.sample<-randCoords.cont.df[sample(1:nrow(randCoords.cont.df),size=nFP),]
    
    ##the random distribution of families is linear over latitude and longitude
    ##this means that families are to be distributed over "geographically ordered" grid points
    #four types of geographic orders are distinguished:
    #sums of lat and long 
    coorRank<-randCoords.sample$coords.x1+randCoords.sample$coords.x2
    #difference between lat and long 
    coorRank<-randCoords.sample$coords.x1-randCoords.sample$coords.x2
    #only lat
    coorRank<-cbind(coorRank,randCoords.sample$coords.x2)
    #only long
    coorRank<-cbind(coorRank,randCoords.sample$coords.x1)
    
    #one of the four options is selected
    coorRank<-coorRank[,sample(1:3,size=1)]
    
    #now it is randomly decided if the geographic order is used in increasing or decreasing order
    up<-sample(0:1,size=1)
    
    #the given geographic order is used to sort all random points
    if(up==1){
      randCoords.sample<-randCoords.sample[order(coorRank,decreasing = F),]
    }else{
      randCoords.sample<-randCoords.sample[order(coorRank,decreasing = T),]
    }
    
    #compute distance between grid points and random points
    lang.rand.dist<-spDists(reg.spdf@coords,as.matrix(randCoords.sample[,4:5]),longlat = T)
    
    #families are randomized and distributed over the geographically ordered random points
    fpFamTab<-fpFamTab[sample(1:nrow(fpFamTab),size=nrow(fpFamTab),replace=F),]
    famRand<-rep(fpFamTab$fpFam,fpFamTab$Freq)
    
    randCoords.sample$randFPfam<-famRand
    
    #minimum distances are used to define nearest neighbours
    minDist<-apply(lang.rand.dist,2,min)
    minDist[minDist>min(regMin)]<-NA
    
    #find closest points an count
    fp.temp<-NULL
    fp.fam.temp<-NULL
    
    for(i in 1:nrow(reg.spdf)){
      fp.temp<-c(fp.temp,sum(lang.rand.dist[i,]==minDist,na.rm=T))
      
      fams<-length(levels(factor(randCoords.sample$randFPfam[lang.rand.dist[i,]==minDist])))
      fp.fam.temp<-c(fp.fam.temp,fams)
    }
    
    fp.all<-rbind(fp.all,fp.temp)
    fpFam.all<-rbind(fpFam.all,fp.fam.temp)
    
  }
  
  #aggregage all information
  fp.all<-fp.all[-1,]
  fp.all[is.na(fp.all)]<-0
  
  fpFam.all<-fpFam.all[-1,]
  fpFam.all[is.na(fpFam.all)]<-0
  
  ##count the number of times the observed language count is higher than the expected language count for each grid point
  fp.smaller<-numeric(nrow(reg.spdf))
  for(i in 1:nrow(reg.spdf)){
    fp.smaller[i]<-sum(reg.spdf$fp[i]>fp.all[,i])/nrow(fp.all)
  }
  fp.smaller[is.na(fp.smaller)]<-0
  
  fpFam.smaller<-numeric(nrow(reg.spdf))
  for(i in 1:nrow(reg.spdf)){
    fpFam.smaller[i]<-sum(reg.spdf$fpFam[i]>fpFam.all[,i])/nrow(fpFam.all)
  }
  fpFam.smaller[is.na(fpFam.smaller)]<-0
  
  
  ##second iteration for HG
  hg.all<-numeric(nrow(reg.spdf))
  hgFam.all<-numeric(nrow(reg.spdf))
  
  for(j in 1:500){
    print(j)
    
    randCoords.sample<-randCoords.cont.df[sample(1:nrow(randCoords.cont.df),size=nHG),]
    
    coorRank<-randCoords.sample$coords.x1+randCoords.sample$coords.x2
    coorRank<-cbind(coorRank,randCoords.sample$coords.x2-randCoords.sample$coords.x1)
    coorRank<-cbind(coorRank,randCoords.sample$coords.x2)
    coorRank<-cbind(coorRank,randCoords.sample$coords.x1)
    
    coorRank<-coorRank[,sample(1:4,size=1)]
    
    up<-sample(0:1,size=1)
    
    if(up==1){
      randCoords.sample<-randCoords.sample[order(coorRank,decreasing = F),]
    }else{
      randCoords.sample<-randCoords.sample[order(coorRank,decreasing = T),]
    }
    
    lang.rand.dist<-spDists(reg.spdf@coords,as.matrix(randCoords.sample[,4:5]),longlat = T)
    
    #randomly distribute families
    hgFamTab<-hgFamTab[sample(1:nrow(hgFamTab),size=nrow(hgFamTab),replace=F),]
    famRand<-rep(hgFamTab$hgFam,hgFamTab$Freq)
    
    randCoords.sample$randHGfam<-famRand
    
    minDist<-apply(lang.rand.dist,2,min)
    
    #find closest points an count
    hg.temp<-NULL
    hg.fam.temp<-NULL
    
    for(i in 1:nrow(reg.spdf)){
      hg.temp<-c(hg.temp,sum(lang.rand.dist[i,]==minDist,na.rm=T))
      
      fams<-length(levels(factor(randCoords.sample$randHGfam[lang.rand.dist[i,]==minDist])))
      hg.fam.temp<-c(hg.fam.temp,fams)
    }
    
    hg.all<-rbind(hg.all,hg.temp)
    hgFam.all<-rbind(hgFam.all,hg.fam.temp)
    
  }
  
  hg.all<-hg.all[-1,]
  hg.all[is.na(hg.all)]<-0
  
  hgFam.all<-hgFam.all[-1,]
  hgFam.all[is.na(hgFam.all)]<-0
  
  ##count the number of times the observed language count is higher than the expected language count for each grid point
  hg.smaller<-numeric(nrow(reg.spdf))
  for(i in 1:nrow(reg.spdf)){
    hg.smaller[i]<-sum(reg.spdf$hg[i]>hg.all[,i])/nrow(hg.all)
  }
  hg.smaller[is.na(hg.smaller)]<-0
  
  hgFam.smaller<-numeric(nrow(reg.spdf))
  for(i in 1:nrow(reg.spdf)){
    hgFam.smaller[i]<-sum(reg.spdf$hgFam[i]>hgFam.all[,i])/nrow(hgFam.all)
  }
  hgFam.smaller[is.na(hgFam.smaller)]<-0
  
  
  #add the information to the grid points
  reg.spdf$hgFamRelSmall<-hgFam.smaller
  reg.spdf$hgRelSmall<-hg.smaller
  reg.spdf$fpFamRelSmall<-fpFam.smaller
  reg.spdf$fpRelSmall<-fp.smaller

  save(reg.spdf,file=paste("output/004_random_NS/randPtsLangPoiss_",resolution,".Rdata",sep=""))
}
