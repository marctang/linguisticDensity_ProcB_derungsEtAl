##################################################
## Name: 002_langPrep  
## Script purpose: Preparing, cleaning and merging the input language data as provided by glottolog and Tom Güldemann
## Date: 2018
## Author: Curdin Derungs
##################################################

library(sp)
rm(list=ls())

##reading the glottolog language points (souce: http://glottolog.org/glottolog/language) 
allLang<-read.csv("input/glottologLangFam.csv")

##reading a list of hunter-gatherer languages as provided by Tom Güldemann
hg<-read.csv("input/hg2017list.csv")
hg.id<-data.frame(iso=hg$ISO,hunter_gatherer=rep("hg",nrow(hg)))

##filtering for languages that have geographic coordinates
allLang<-allLang[!is.na(allLang$longitude),]
allLang<-allLang[!is.na(allLang$latitude),]

##only attested languages are used in the analysis
allLang<-allLang[allLang$glottolog.stock!="Unattested [unat1236]",]

##merging the glottolog and hunter-gather lists using the iso nr.
allLang<-merge(allLang,hg.id,by.y="iso",by.x="ISO.g",all.x=T)

##preparing variables
allLang$hunter_gatherer<-as.character(allLang$hunter_gatherer)
allLang$hunter_gatherer[is.na(allLang$hunter_gatherer)]<-"not hg"
allLang$hunter_gatherer<-factor(allLang$hunter_gatherer)

##add binfords data that is missing in our hunter-gatherer list
#A preprocessing step was required to filter the Binford list as provided by the R Package <Binford>)
binf<-read.csv("input/binford.csv")
allLang<-rbind(allLang,binf)

##splitting the data into "hg" and "not hg"
fp<-allLang[allLang$hunter_gatherer=="not hg",]
hg<-allLang[allLang$hunter_gatherer!="not hg",]

##convert languages to spatial point data
coordinates(allLang)<-~longitude+latitude
proj4string(allLang)<-CRS("+proj=longlat +datum=WGS84")

#save spatial language data
save(allLang,file="output/allLang.Rdata")
