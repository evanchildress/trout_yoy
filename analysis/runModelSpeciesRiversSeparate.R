source('~/trout_yoy/createModelFunctionSpeciesRiversSeparate.R')
source("~/trout_yoy/analysis/runModelFunctionSpeciesRiversSeparate.R")
for(species in c('bkt','bnt')){
  for(r in 1:4){
    for(sr in c(T,F)){
      for(e in c('none',"mean","extreme","pca")){
        randomYear<-F
        if(!sr & e=='none'){randomYear<-T}
        runModel(stockRecruit=sr,env=e,randomYear,species=species,river=r)
      }
    }
  }
}
