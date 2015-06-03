source('~/trout_yoy/createModelFunctionSpeciesSeparate.R')
source("~/trout_yoy/analysis/runModelFunctionSpeciesSeparate.R")
for(species in c('bkt','bnt')){
  for(sr in c(T,F)){
    for(e in c('none',"mean","extreme","pca")){
      randomYear<-F
      if(!sr & e=='none'){randomYear<-T}
      runModel(stockRecruit=sr,env=e,randomYear,species=species)
    }
  }
}
