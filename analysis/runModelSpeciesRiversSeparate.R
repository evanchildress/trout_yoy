source('~/trout_yoy/createModelFunctionSpeciesRiversSeparate.R')
source("~/trout_yoy/analysis/runModelFunctionSpeciesRiversSeparate.R")
for(species in c('bkt','bnt')){
  for(r in 1:4){
    for(sr in c(T,F)){
      for(e in c('none',"mean","extreme","pca")){
        randomYear<-F
        if(!sr & e=='none'){randomYear<-T}
        if(species == 'bnt' & r %in% 2:3){
          cat(species,r,sr,e)
          next}
        cat(species,r,ifelse(sr,"sr","no sr"),e,"\n","\n")
        runModel(stockRecruit=sr,env=e,randomYear,
                 species=species,river=r,highResTemp=T)
        
      }
    }
  }
}

source('~/trout_yoy/createModelFunctionSpeciesRiversSeparate.R')
source("~/trout_yoy/analysis/runModelFunctionSpeciesRiversSeparate.R")
for(species in c('bkt')){
  for(r in c(1,4)){
    for(sr in c(T)){
      for(e in c("extreme")){
        randomYear<-F
        if(!sr & e=='none'){randomYear<-T}
        if(species == 'bnt' & r %in% 2:3){
          cat(species,r,sr,e)
          next}
        cat(species,r,ifelse(sr,"sr","no sr"),e,"\n","\n")
        runModel(stockRecruit=sr,env=e,randomYear,
                 species=species,river=r,highResTemp=T)
        
      }
    }
  }
}