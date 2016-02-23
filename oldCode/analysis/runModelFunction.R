runModel<-function(stockRecruit=F,env='none',randomYear=F){
  resultsDir<-"~/trout_yoy/results/"
  
  outName<-ifelse(randomYear,'randomYear',
                  ifelse(env=='none','stockRecruit',
                         ifelse(stockRecruit,paste0(env,'StockRecruit'),env)))

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
  
  covPCA<-array(NA,dim=dim(covariates))
  if(env=='pca'){
    for(r in 1:4){
    pca<-prcomp(covariates[,r,],center=T,scale=T)
    scores<-pca$x
#     loadings<-pca$rotation
#     loadings[,1]<-loadings[,1]
    covPCA[,r,]<-scores
    }
    
  covariates<-covPCA
#     tiff.par(file="~/trout_yoy/results/figures/envPca.tif",
#          height=4,width=8,mai=c(.8,.8,0,0),mar=c(3,3,0,0),
#          mfrow=c(1,2))
#     
#       plot(scores[,1], scores[,2], xlab="", ylab="",
#       type="n", xlim=c(-8,8),
#       ylim=c(-8,5),bty='l')
#       arrows(0,0,loadings[,1]*10,loadings[,2]*10, length=0.1,
#       angle=20, col="red")
#       title(xlab="PC1",ylab="PC2",line=2)
#       text(loadings[,1]*10*1.2,loadings[,2]*10*1.2,
#       rownames(loadings), col="red", cex=0.7)
#       points(scores[,1],scores[,2], col="blue",cex=0.7,type="p")
#     
#       plot(scores[,3], scores[,4], xlab="", ylab="",
#       type="n", xlim=c(-8,10),
#       ylim=c(-6,8),bty='l')
#       arrows(0,0,loadings[,3]*10,loadings[,4]*10, length=0.1,
#       angle=20, col="red")
#       title(xlab="PC3",ylab="PC4",line=2)
#       text(loadings[,3]*10*1.2,loadings[,4]*10*1.2,
#       rownames(loadings), col="red", cex=0.7)
#       points(scores[,3],scores[,4], col="blue",cex=0.7,type="p")
#     
#     dev.off()
    
  }
  
  N<-y #use naive estimate as initial value
  
  win.data<-list(yDATA=y,
                 nsamples=dim(y)[1],
                 covariates=covariates,
                 p=detection,
                 adultDATA=A)
  
  inits<-function(){
  
    if(stockRecruit & env %in% c('mean','extreme')){
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
    
    if(stockRecruit & env=='pca'){
      return(
        list(c=array(c(runif(8,0,5),rep(0,16)),dim=c(4,2,3)),
             beta=array(rnorm(72,0,0.01),dim=c(4,9,2)),
             N=N
        )
      )
    }
    
    if(!stockRecruit & env=='pca'){
      return(
        list(beta=array(rnorm(72,0,0.01),dim=c(4,9,2)),
             N=N
        )
      )
    }
    
    if(!stockRecruit & env %in% c('mean','extreme')){
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
          params<-c("N","c","beta","fit","fitNew","yExp")
        } else {
          if(randomYear){
            params<-c("N","c","eps","fit","fitNew","yExp")
          }
        }
      }
    }
  
  ni=15000
  nt=5
  nb=10000
  nc=3
  
  out<-jags(win.data,inits,params,"~/trout_yoy/abundanceModel.txt",n.chains=nc,n.iter=ni,
            n.thin=nt,n.burnin=nb,working.directory=getwd())
  saveRDS(out,file.path(resultsDir,paste0("modelOutput/",outName,".rds")))
  
  sims<-out$BUGSoutput$sims.list
  
################### make posterior predictive check figure ###################
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
#####################################################################  
}