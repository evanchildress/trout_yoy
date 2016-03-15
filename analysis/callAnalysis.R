require(tidyr)
require(getWBData)
require(data.table)
require(dplyr)

# rmseOut<-list()
# for(bla in c("flowExtension","dailyDischarge")){
qType="dailyDischarge"
for(whichSpecies in c("bkt","bnt")){
source("analysis/loadDataAndCjsResults.R")
source("analysis/covariateInputs.R")
source("analysis/createCovariates.R")
}

for(whichSpecies in c("bkt","bnt")){
jagsData<-readRDS(paste0("cjsInputs/jagsData",
                            toupper(substr(whichSpecies,1,1)),
                            substr(whichSpecies,2,nchar(whichSpecies)),
                            ".rds"))
source("analysis/createModel.R")
source("analysis/runModelFunction.R")


  for(sr in c(T,F)){
      for(e in c("extreme","mean","none")){
        randomYear<-F
        if(!sr & e=='none'){randomYear<-T}
        cat(whichSpecies,ifelse(sr,"sr","no sr"),e,"\n","\n")
        runModel(jagsData,stockRecruit=sr,env=e,randomYear,
                 species=whichSpecies)
    }
  }
}
source("analysis/compareRmse.R")
source("results/figures/manuscript/envCovariates.r")
