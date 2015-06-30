fileDir<-'~/trout_yoy/results/modelOutput/speciesRiverSeparate'
setwd(fileDir)

summarizeFit<-function(fileName){
  #load model output
  data<-readRDS(fileName)$BUGSoutput
  
  #identify the species, river, and model type its from
  sp<-ifelse(grepl('Bkt',fileName),1,2)
  species<-c("Bkt","Bnt")[sp]
  
  river<-strsplit(strsplit(fileName,species)[[1]][2],".rds")[[1]]
  r<-which(river==c("Jimmy","Mitchell","Obear","WestBrook"))
  
  modelType<-strsplit(fileName,species)[[1]][1]
  
  #load observed data and reduce to the appropriate river and species
  load('~/trout_yoy/abundanceData.RData')
  y<-y[3:15,r,sp,1]
  
  nYears<-13
  yExp<-data$sims.list$yExp
 
  #calculate mean rmse
  rmse<-mean(sqrt(apply(sweep(yExp,2,y)^2,1,sum)/nYears))
  
  #get DIC
  dic<-data$DIC
  
  #get pCheck
  pCheck<-mean(data$sims.list$pCheck)
  
  #collate and return results
  result<-data.frame(river,species,modelType,rmse,dic,pCheck)
  names(result)<-c("river","species","modelType","rmse","dic","pCheck")
  return(result)
}  

toEvaluate<-list.files(fileDir)

results<-data.table("river"=as.character(rep(NA,length(toEvaluate))),
                    "species"=as.character(rep(NA,length(toEvaluate))),
                    "modelType"=as.character(rep(NA,length(toEvaluate))),
                    "rmse"=as.numeric(rep(NA,length(toEvaluate))),
                    "dic"=as.numeric(rep(NA,length(toEvaluate))),
                    "pCheck"=as.numeric(rep(NA,length(toEvaluate))))

for(d in toEvaluate){
  results[which(d==toEvaluate),]<-summarizeFit(d)[1,]
}

setkey(results,species,river,dic)

forCoef<-c(grep('meanStockRecruit',toEvaluate),
  grep('extremeStockRecruit',toEvaluate))
