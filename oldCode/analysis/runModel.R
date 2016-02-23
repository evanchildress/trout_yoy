source('~/trout_yoy/createModelFunction.R')
source("~/trout_yoy/analysis/runModelFunction.R")

for(sr in c(T,F)){
  for(e in c('none',"mean","extreme","pca")){
    randomYear<-F
    if(!sr & e=='none'){randomYear<-T}
    runModel(stockRecruit=sr,env=e,randomYear)
  }
}

