setwd('~/trout_yoy/results/modelOutput')

rivers<-c('jimmy','mitchell','obear','west brook')

rmse<-function(fileName,fileDir){
  data<-readRDS(file.path(fileDir,fileName))$BUGSoutput$sims.list
  
  load('~/trout_yoy/abundanceData.RData')
  y<-y[3:15,,,1]
  
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

dic<-function(fileName){
  dic<-readRDS(fileName)$BUGSoutput$DIC
  return(dic)
}

pCheck<-function(fileName){
  data<-readRDS(fileName)$BUGSoutput$sims.list$pCheck
  return(mean(data))
}

toEvaluate<-list.files(getwd())
rmseOut<-list()
dicOut<-NULL
for(d in toEvaluate){
  rmseOut[[d]]<-rmse(d)
  dicOut[which(d==toEvaluate)]<-dic(d)
  pCheckOut[which(d==toEvalute)]<-pCheck(d)
}
#rmseOut$model<-toEvaluate

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

rmseSummary<-array(dim=c(length(toEvaluate),4,2))
summarize<-function(x){
  return(c(mean(x),quantile(x,probs=c(0.025,0.975))))
}
for(d in toEvaluate){
  rmseSummary[which(d==toEvaluate),,1]<-apply(rmseOut[[d]]$bkt,1,mean)
  rmseSummary[which(d==toEvaluate),,2]<-apply(rmseOut[[d]]$bnt,1,mean)
}
dimnames(rmseSummary)<-list(strsplit(toEvaluate,'.rds'),rivers,c('bkt','bnt'))
