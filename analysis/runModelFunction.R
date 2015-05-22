runModel<-function(stockRecruit=F,env='none',randomYear=F){
  resultsDir<-"~/trout_yoy/results/"
  
  outName<-ifelse(randomYear,'randomYear',
                  ifelse(env=='none','stockRecruit',
                         ifelse(stockRecruit,paste0(env,'StockRecruit'),env)))

  if(!exists("createModel")){source('~/trout_yoy/createModelFunction.R')}
  createModel(stockRecruit=stockRecruit,env=env,randomYear=randomYear)
  
     
  load('~/trout_yoy/abundanceData.RData')
  
  sumTempMean<-data.table(covariates[,,9])
  sumTempFill<-lm(V4~V3+V2+V1,data=sumTempMean)
  covariates[6:9,4,9]<-predict(sumTempFill,data.frame(sumTempMean[6:9,list(V3,V2,V1)]))
  rm(sumTempMean)
  rm(sumTempFill)
  
  y<-y[3:15,,,1]
  
  covariates<-covariates[1:13,,]
  detection<-detection[1:13,,]
  A<-A[1:14,,]
  
  N<-y #use naive estimate as initial value
  
  win.data<-list(yDATA=y,
                 nsamples=dim(y)[1],
                 covariates=covariates,
                 p=detection,
                 adultDATA=A)
  
  inits<-function(){
  
    if(stockRecruit & !env=='none'){
      return(
        list(c=array(c(runif(8,0,5),rep(0,16)),dim=c(4,2,3)),
                           beta=array(rnorm(64,0,0.01),dim=c(4,8,2)),
                           N=N
        )
      )
    }
    
    if(stockRecruit & env=='none'){
      return(
        list(c=array(c(runif(8,0,5),rep(0,16)),dim=c(4,2,3)),
                          N=N
        )
      )
    }
    
    if(!stockRecruit & env!='none'){
      return(
        list(beta=array(rnorm(64,0,0.01),dim=c(4,8,2)),
              N=N
        )
      )
    }
    
    if(randomYear){return(list(N=N))}
  }
  
  
  if(stockRecruit & env!='none'){
    params<-c("N","c","beta","fit","fitNew","yExp")
    } else {
      if(stockRecruit & env=='none'){
        params<-c("N","c","fit","fitNew","yExp")
      } else {
        if(!stockRecruit & env!='none'){
          params<-c("N","beta","fit","fitNew","yExp")
        } else {
          if(randomYear){
            params<-c("N","eps","fit","fitNew","yExp")
          }
        }
      }
    }
  
  ni=15000
  nt=5
  nb=10000
  nc=3
  
  out<-jags(win.data,inits,params,"abundanceModel.txt",n.chains=nc,n.iter=ni,
            n.thin=nt,n.burnin=nb,working.directory=getwd())
  saveRDS(out,file.path(resultsDir,paste0("modelOutput/",outName,".rds")))
  
  sims<-out$BUGSoutput$sims.list
  
  
  tiff.par(file.path(resultsDir,paste0("figures/bayesP_",outName,".tif")),
           mfrow=c(1,2),width=8,height=4,mar=c(2.5,2.5,1,0))
  plot(sims$fitNew[,,1]~sims$fit[,,1],
       xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
       bty='l',main='bkt')
  for(r in 1:4){
    points(sims$fitNew[,r,1]~sims$fit[,r,1],col=palette()[r])
  }
  abline(0,1,lwd=2)
  text(max(sims$fit[,,1])*0.9,max(sims$fitNew[,,1])*0.9,
       round(mean(sims$fitNew[,,1]>sims$fit[,,1]),3))
  
  plot(sims$fitNew[,,2]~sims$fit[,,2],
       xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
       bty='l',main='bnt',pch=NA)
  for(r in 1:4){
    points(sims$fitNew[,r,2]~sims$fit[,r,2],col=palette()[r])
  }
  abline(0,1,lwd=2)
  text(max(sims$fit[,,2])*0.9,max(sims$fitNew[,,2])*0.9,
              round(mean(sims$fitNew[,,2]>sims$fit[,,2]),3))
  dev.off()
}