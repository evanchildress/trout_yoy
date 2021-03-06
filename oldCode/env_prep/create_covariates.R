env.cov<-function(covariate,month,threshold=NA,high.low=NA,
                  freq.dur=NA,FUN=NULL){
  
  if(!is.null(FUN) & any(!is.na(threshold),!is.na(high.low),!is.na(freq.dur))){
    stop("cannot define both 'FUN' and threshold/high.low/freq.dur")
  }  

  if(!covariate %in% c('discharge','temperature')){
    stop("covariate must equal 'discharge' or 'temperature'")}
  
  if(covariate=='discharge'){data<-dailyDischarge}
  if(covariate=='temperature'){data<-ed}
  
  data<-data[date >= as.Date(paste0(start_year-1,"-10-01")) & 
         date < as.Date(paste0(end_year,"-10-01"))]
  data[,year_of_effect:=year(date)]
  data[month(date)>=10,year_of_effect:=year(date)+1]
  
  rivers<-unique(data$river)
  
  result<-array(NA,dim=c(length(unique(year(data$date)))-1,4))
  
  if(is.null(FUN)){


  if(!high.low %in% c('high','low')){
    stop("high.low must equal 'high' or 'low' defining whether the event is a high or low extreme")}
  
  if(!freq.dur %in% c('frequency','duration')){
    stop("freq.dur must equal 'frequency' or 'duration' defining whether the function returns the number of days when the threshold is crossed (duration) or the number of times(frequency)")}

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
    if(freq.dur=="frequency"){result[,r]<-data[river==rivers[r]&month(date) %in% month,
                                             ifelse(high.low=="high",
                                             num.sets(which(get(covariate) >= limit[r])),
                                             num.sets(which(get(covariate) <= limit[r]))),
                                             by=year_of_effect]$V1
    }
    
    if(freq.dur=="duration"){result[,r]<-data[river==rivers[r]&month(date) %in% month,
                                             ifelse(high.low=="high",
                                                    length(which(get(covariate) >= limit[r])),
                                                    length(which(get(covariate) <= limit[r]))),
                                             by=year_of_effect]$V1
    }
  }} else{ #if FUN isn't null
    
    for(r in 1:4){
      result[,r]<-data[river == rivers[r] & month(date) %in% month,
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

covInt<-function()

covariates<-array(dim=c(end_year-start_year+1,4,length(names(covariate_inputs))+2))
dimnames(covariates)<-list(start_year:end_year,
                           unique(ed$rivers),
                           c(names(covariate_inputs),
                             "summerTempHighRes",
                             "summerTempMeanHighRes"))

for(i in 1:length(names(covariate_inputs))){
  covariates[,,i]<-do.call(env.cov,covariate_inputs[[i]])
}

covariates[,,"summerTempHighRes"]<-apply(summerT,2,
                                   function(x){return(scale(x)[,1])})
covariates[,,"summerTempMeanHighRes"]<-apply(meanSummerT,2,
                                   function(x){return(scale(x)[,1])})


assign('covariates',covariates,env=shared_data)
saveRDS(covariates,"~/process-data/data_store/processed_data/covariates.rds")


