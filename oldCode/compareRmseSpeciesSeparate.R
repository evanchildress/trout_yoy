fileDir<-'~/trout_yoy/results/modelOutput/speciesSeparate'
setwd(fileDir)

rivers<-c('jimmy','mitchell','obear','west brook')

rmse<-function(fileName){
  data<-readRDS(fileName)$BUGSoutput$sims.list
  
  load('~/trout_yoy/abundanceData.RData')
  y<-y[3:15,,ifelse(grepl('Bkt',fileName),1,2),1]
  
  nYears<-13
  nRivers<-4
  yExp<-data$yExp
  n<-dim(yExp)[1]
  
  byRow<-function(x){
    result<-NULL
    for(r in 1:nRivers){
      result[r]<-sqrt(sum((x[,r]-melt(y[,r])$value)^2)/(nYears))
    }
    return(result)
  }
  rmse<-apply(yExp,1,byRow)
  
  return(rmse)
}

dic<-function(fileName){
  dic<-readRDS(fileName)$BUGSoutput$DIC
  return(dic)
}

pCheck<-function(fileName){
  data<-readRDS(fileName)$BUGSoutput$sims.list$pCheck
  return(apply(data,2,mean))
}

toEvaluate<-list.files(fileDir)
rmseOut<-list()
dicOut<-NULL
pCheckOut<-array(NA,dim=c(length(toEvaluate),4))
for(d in toEvaluate){
  rmseOut[[d]]<-rmse(d)
  dicOut[which(d==toEvaluate)]<-dic(d)
  pCheckOut[which(d==toEvaluate),]<-pCheck(d)
}
#rmseOut$model<-toEvaluate

# limits<-c(20,20,50,80)

# tiff.par("~/trout_yoy/results/figures/rmse.tif",
#          mfrow=c(2,2),mar=c(2.5,2.5,1,0))
# for(r in 1:4){
#   plot(NA,xlim=c(0,limits[r]),ylim=c(0,0.5),
#        main=rivers[r],xlab='rmse (# yoy)',
#        ylab='density')
#   for(d in toEvaluate){
#     points(density(rmseOut[[d]]$bkt[r,]),type='l',
#          col=palette()[which(d==toEvaluate)])
# }
# }
# dev.off()

rmseSummary<-array(dim=c(length(toEvaluate),4))
summarize<-function(x){
  return(c(mean(x),quantile(x,probs=c(0.025,0.975))))
}
for(d in toEvaluate){
  if(!is.null(rmseOut[[d]])){
  rmseSummary[which(d==toEvaluate),]<-apply(rmseOut[[d]],1,mean)
}}
dimnames(rmseSummary)<-list(strsplit(toEvaluate,'.rds'),rivers)
dimnames(pCheckOut)<-list(strsplit(toEvaluate,'.rds'),rivers)



