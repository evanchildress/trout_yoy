setwd('~/trout_yoy/results')

rmse<-function(dataName){
  data<-readRDS(paste0('modelOutput',dataName,'.rds'))$BUGSoutput$sims.list
  
  nYears<-13
  nRivers<-4
  yExpBkt<-data$yExp[,,,1]
  yExpBnt<-data$yExp[,,,2]
  n<-dim(yExpBkt)[1]
  
  byRow<-function(x,sp=1){
    result<-NULL
    for(r in 1:nRivers){
      result[r]<-sqrt(sum((x[,r]-melt(y[,r,sp])$value)^2)/(nYears))
    }
    return(result)
  }
  rmseBkt<-apply(yExpBkt,1,byRow)
  rmseBnt<-apply(yExpBnt,1,byRow,sp=2)
  
  #rmseBkt<-sqrt(sum((yExpBkt-melt(y[,,1])$value)^2)/(nYears*nRivers))
  #rmseBnt<-sqrt(sum((yExpBnt-melt(y[,,2])$value)^2)/(nYears*nRivers))
  
  rmse<-list('bkt'=rmseBkt,'bnt'=rmseBnt)
  return(rmse)
}

toEvaluate<-c('Random','StockRecruit','Extreme','MeanEnv')
rmseOut<-list()
for(d in toEvaluate){
  rmseOut[[d]]<-rmse(d)
}
#rmseOut$model<-toEvaluate
rivers<-c('jimmy','mitchell','obear','west brook')
limits<-c(20,20,50,80)

tiff.par("~/trout_yoy/results/figures/rmse.tif",
         mfrow=c(2,2),mar=c(2.5,2.5,1,0))
for(r in 1:4){
  plot(NA,xlim=c(0,limits[r]),ylim=c(0,0.5),
       main=rivers[r],xlab='rmse (# yoy)',
       ylab='density')
  for(d in toEvaluate){
    points(density(rmseOut[[d]]$bkt[r,]),type='l',
         col=palette()[which(d==toEvaluate)])
}
}
dev.off()
