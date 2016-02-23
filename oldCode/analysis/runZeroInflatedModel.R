source("~/trout_yoy/zeroInflatedModelMitchellBnt.R")  
     
  load('~/trout_yoy/abundanceData.RData')
  
  sumTempMean<-data.table(covariates[,,9])
  sumTempFill<-lm(V4~V3+V2+V1,data=sumTempMean)
  covariates[6:9,4,9]<-predict(sumTempFill,data.frame(sumTempMean[6:9,list(V3,V2,V1)]))
  rm(sumTempMean)
  rm(sumTempFill)
  
  y<-y[3:15,2,2,1]
  
  covariates<-covariates[1:13,2,]
  detection<-detection[1:13,2,2]
  A<-A[1:14,2,2]
  
  N<-y #use naive estimate as initial value
  
  win.data<-list(yDATA=y,
                 nsamples=length(y),
                 covariates=covariates,
                 p=detection,
                 adultDATA=A)
  
  inits<-function(){
        list(N=N
        )
  }
    

  
  

  params<-c("N","c","beta","zbeta","fit","fitNew","yExp",'pCheck')
  
  ni=15000
  nt=5
  nb=10000
  nc=3
  
  resultsDir<-"~/trout_yoy/results"
  out<-jags(win.data,inits,params,"~/trout_yoy/abundanceModel.txt",n.chains=nc,n.iter=ni,
            n.thin=nt,n.burnin=nb)
  saveRDS(out,file.path(resultsDir,paste0("modelOutput/SpeciesRiverSeparate/","extremeBntMitchell",".rds")))
  
  sims<-out$BUGSoutput$sims.list
