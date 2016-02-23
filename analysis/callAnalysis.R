# rmseOut<-list()
# for(bla in c("flowExtension","dailyDischarge")){
qType="dailyDischarge"
whichSpecies="bnt"

source("analysis/compareRmse.R")
source("analysis/loadDataAndCjsResults.R")
source("analysis/covariateInputs.R")
source("analysis/createCovariates.R")
# if(!exists("jagsData")){
#   jagsData<-readRDS("cjsInputs/jagsData.rds")
# }
source("analysis/createModel.R")
source("analysis/runModelFunction.R")

for(species in c(whichSpecies)){
  for(sr in c(T,F)){
      for(e in c("extreme","mean","none")){
        randomYear<-F
        if(!sr & e=='none'){randomYear<-T}
        if(species == whichSpecies){
          cat(species,sr,e)
          next}
        cat(species,ifelse(sr,"sr","no sr"),e,"\n","\n")
        runModel(jagsData,stockRecruit=sr,env=e,randomYear,
                 species=species)
    }
  }
}

for(e in 1:length(toEvaluate)){
  rmseSummary[e,]<-apply(rmse(toEvaluate[e]),1,mean)
}

