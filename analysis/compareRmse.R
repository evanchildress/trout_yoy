nRivers<-c(4,3)[which(whichSpecies==c("bkt","bnt"))]

y<-readRDS(paste0("cjsInputs/jagsData",toupper(substr(whichSpecies,1,1)),substr(whichSpecies,2,nchar(whichSpecies)),".rds"))$y

rmse<-function(fileName,standardize=T){
  data<-readRDS(paste0("results/modelOutput/",fileName))$BUGSoutput$sims.list
  
  randomN<-readRDS(paste0("results/modelOutput/randomYear",whichSpecies,".rds"))$BUGSoutput$sims.list$N %>%
    apply(3,mean)
  
  nYears<-14
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
  if(standardize){rmse<-rmse/randomN}
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

# toEvaluate<-list.files(fileDir)
# rmseOut<-list()
# dicOut<-NULL
# pCheckOut<-array(NA,dim=c(length(toEvaluate),4))
# for(d in toEvaluate){
#   rmseOut[[d]]<-rmse(d)
#   dicOut[which(d==toEvaluate)]<-dic(d)
#   pCheckOut[which(d==toEvaluate),]<-pCheck(d)
# }
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
makeFileName<-function(models,species){
  species<-paste0(toupper(substr(species,1,1)),substr(species,2,nchar(species)))
  return(paste0(models,species,".rds"))
}

toEvaluate<-makeFileName(c("randomYear","mean","extreme","stockRecruit",
              "meanStockRecruit","extremeStockRecruit"),whichSpecies)

rmseSummary<-array(NA,dim=c(length(toEvaluate),nRivers+1))
rownames(rmseSummary)<-unlist(strsplit(toEvaluate,".rds"))
summarize<-function(x){
  return(c(mean(x),quantile(x,probs=c(0.025,0.975))))
}

dic<-array(NA,dim=c(length(toEvaluate),1))
rownames(dic)<-unlist(strsplit(toEvaluate,".rds"))

for(m in 1:length(toEvaluate)){
  overallRmse<- mean(readRDS(paste0("results/modelOutput/",toEvaluate[m]))$BUGSoutput$sims.list$overallRmse)/
    mean(apply(readRDS(paste0("results/modelOutput/",toEvaluate[m]))$BUGSoutput$sims.list$N,c(1,2),sum))
  rmseSummary[m,]<-c(apply(rmse(toEvaluate[m]),1,mean),overallRmse)
  dic[m]<-readRDS(paste0("results/modelOutput/",toEvaluate[m]))$BUGSoutput$DIC
}








