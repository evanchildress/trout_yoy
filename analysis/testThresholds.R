t<-readRDS("cjsInputs/tempData.rds")
q<-readRDS("cjsInputs/dischargeData.rds")
source("analysis/frequentistModelSelection.R")
start_year<-2002
end_year<-2015
nRivers<-4

rivers<-c("west brook","wb jimmy","wb mitchell","wb obear")

env.cov<-function(covariate,days,threshold=NA,high.low=NA,
                  freq.dur=NA,FUN=NULL){
  
  if(!is.null(FUN) & any(!is.na(threshold),!is.na(high.low),!is.na(freq.dur))){
    stop("cannot define both 'FUN' and threshold/high.low/freq.dur")
  }  
  
  if(!covariate %in% c('discharge','temperature')){
    stop("covariate must equal 'discharge' or 'temperature'")}
  
  if(covariate=='discharge'){data<-q}
  if(covariate=='temperature'){data<-t}
  
  data<-data[date >= as.Date(paste0(start_year-1,"-10-01")) & 
               date < as.Date(paste0(end_year,"-10-01"))]
  data[,year_of_effect:=year(date)]
  data[month(date)>=10,year_of_effect:=as.integer(year(date)+1)]
  data[,day:=as.numeric(date-as.Date(paste0(year_of_effect-1,"-09-30")))]
  
  result<-array(NA,dim=c(length(unique(year(data$date)))-1,4))
  
  if(is.null(FUN)){
    
    
    if(!high.low %in% c('high','low')){
      stop("high.low must equal 'high' or 'low' defining whether the event is a high or low extreme")}

    if(covariate == 'discharge' & (is.na(threshold) | threshold<0 | threshold > 1)){
      stop("for discharge covariates threshold must be a number between 0 and 1 representing the quantile definition of an extreme event")
    }
    
    num.sets<-function(x){
      if(length(x)==0){return(as.integer(0))} else
        xlag<-c(-1,x[1:(length(x)-1)])
      return(length(which(x-xlag>1)))
    }
    
    
    limit<-NULL
    for(r in 1:4){
      if(covariate=="discharge") {limit[r]<-quantile(data[river==rivers[r],get(covariate)],probs=threshold,na.rm=T)}
      if(covariate=="temperature") {limit[r]<-threshold}
    }
    
    for(r in 1:4){
      if(freq.dur=="frequency"){result[,r]<-data[river==rivers[r]& day %in% days,
                                                 ifelse(high.low=="high",
                                                        num.sets(which(get(covariate) >= limit[r])),
                                                        num.sets(which(get(covariate) <= limit[r]))),
                                                 by=year_of_effect]$V1
      }
      
      if(freq.dur=="duration"){result[,r]<-data[river==rivers[r]& day %in% days,
                                                ifelse(high.low=="high",
                                                       length(which(get(covariate) >= limit[r])),
                                                       length(which(get(covariate) <= limit[r]))),
                                                by=year_of_effect]$V1
      }
      if(freq.dur=="sumAbove"){result[,r]<-data[river==rivers[r]& day %in% days,
                                                ifelse(high.low=="high",
                                                       sum(get(covariate)[which(get(covariate) >= limit[r])]-limit[r]),
                                                       sum(limit[r]-get(covariate)[which(get(covariate) <= limit[r])])),
                                                by=year_of_effect]$V1
      }
    }} else{ #if FUN isn't null
      
      for(r in 1:4){
        result[,r]<-data[river == rivers[r] & day %in% days,
                         FUN(get(covariate)),
                         by = year_of_effect]$V1
      }
    }
  
  scale_simple<-function(x) {
    a<-scale(x)
    return(a[,1])}
  result<-apply(result,2,scale_simple)
  
  return(result)
}

for(dur in c(45)){
duration<-dur
startDays<-seq(1,366-duration,2)
thresholds<-seq(0.01,0.99,0.01)

rsq<-aic<-slope<-p<-array(NA,dim=c(length(startDays),length(thresholds),nRivers))
rsqMean<-aicMean<-slopeMean<-pValueMean<-array(NA,dim=c(length(startDays),nRivers))
pb<-txtProgressBar(0,nRivers*length(startDays),style=3)
for(riv in rivers){
  column<-which(riv==rivers)
  yoy<-data[river==riv,logPerCapita]
  adults<-jagsData$adultDATA[,column]
  
  for(start in startDays){
    days<-start:(start+duration)
    
    for(thresh in thresholds){
      e<-env.cov("discharge",threshold=thresh,high.low=ifelse(thresh>0.5,"high","low"),freq.dur="duration",days=days)
      if(all(is.na(e[,column]))) next
      a<-lm(yoy~e[,column]+adults)
      rsq[which(start==startDays),which(thresh==thresholds),column]<-summary(a)$r.squared
      aic[which(start==startDays),which(thresh==thresholds),column]<-AICc(a)
      slope[which(start==startDays),which(thresh==thresholds),column]<-coef(a)[2]
      p[which(start==startDays),which(thresh==thresholds),column]<-summary(a)$coefficients[11]
    }
    meanE<-env.cov("discharge",days=days,FUN=mean)
    meanModel<-lm(yoy~meanE[,column]+adults)
    rsqMean[which(start==startDays),column]<-summary(meanModel)$r.squared
    aicMean[which(start==startDays),column]<-AICc(meanModel)
    pValueMean[which(start==startDays),column]<-summary(meanModel)$coefficients[11]
   setTxtProgressBar(pb,which(start==startDays)+(column-1)*length(startDays))
  }
}

assign(paste0("pValue",dur),p)
assign(paste0("rsq",dur),rsq)
assign(paste0("aic",dur),aic)
assign(paste0("slope",dur),slope)
assign(paste0("rsqMean",dur),rsqMean)
assign(paste0("aicMean",dur),aicMean)
assign(paste0("slopeMean",dur),slopeMean)
}

rm(list=c("slope","rsq","aic","p"))
toSave<-c(ls()[grepl("rsq",ls())],ls()[grepl("aic",ls())],ls()[grepl("slope",ls())],
          ls()[grepl("pValue",ls())])
save(toSave,file="results/thresholdResults.rdata")



