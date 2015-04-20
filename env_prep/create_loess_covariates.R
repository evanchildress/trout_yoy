covariate_inputs<-list(fall_discharge = list(covariate = "discharge",
                                             month = 10:12,
                                             FUN = median),
                       spring_floods =  list(covariate = "discharge",
                                             threshold = 0.99,
                                             high.low = "high",
                                             freq.dur = "duration",
                                             month = 2:3),
                       winter_floods =  list(covariate = "discharge",
                                             threshold = 0.99,
                                             high.low = "high",
                                             freq.dur = "duration",
                                             month = c(12,1)),
                       summer_low =     list(covariate = "discharge",
                                             threshold = 0.05,
                                             high.low = "low",
                                             freq.dur = "duration",
                                             month = 6:9),
                       summer_temp =    list(covariate = "temperature",
                                             threshold = 20,
                                             high.low = "high",
                                             freq.dur = "duration",
                                             month = 6:9),
                       spring_flow_mean =  list(covariate = "discharge",
                                             month = 2:3,
                                             FUN = mean),
                       winter_flow_mean =  list(covariate = "discharge",
                                             month = c(12,1),
                                             FUN = mean),
                       summer_flow_mean =  list(covariate = "discharge",
                                           month = c(6:8),
                                           FUN = mean),
                       summer_temp_mean =  list(covariate = "temperature",
                                           month = c(6:8),
                                           FUN = mean)
)

env.cov<-function(covariate, month, threshold=NA,
                  high.low=NA, freq.dur=NA, FUN=NULL){
#Creates covariates from environmental data
  
  if(!is.null(FUN) & any(!is.na(threshold), !is.na(high.low), !is.na(freq.dur))){
    stop("cannot define both 'FUN' and threshold/high.low/freq.dur")
  }  

  if(!covariate %in% c('discharge', 'temperature')){
    stop("covariate must equal 'discharge' or 'temperature'")}
  
  edLoess<-edLoess[date >= as.POSIXct(paste0(start_year-1,"-10-01")) & 
         date < as.POSIXct(paste0(end_year,"-10-01"))]
  edLoess[,year_of_effect:=year(date)]
  edLoess[month(date)>=10,year_of_effect:=year(date)+1]
  
  rivers<-unique(edLoess$river)
  
  result<-array(NA,dim=c(length(unique(year(edLoess$date)))-1,4))
  
  if(is.null(FUN)){


  if(!high.low %in% c('high','low')){
    stop("covariate must equal 'high' or 'low' defining whether the event is a high or low extreme")}
  
  if(!freq.dur %in% c('frequency','duration')){
    stop("covariate must equal 'frequency' or 'duration' defining whether the function returns the number of days when the threshold is crossed (duration) or the number of times(frequency)")}

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
    limit[r]<-quantile(edLoess[river==rivers[r],get(covariate)],probs=threshold,na.rm=T)
  }
  
  for(r in 1:4){
    if(freq.dur=="frequency"){result[,r]<-edLoess[river==rivers[r]&month(date) %in% month,
                                             ifelse(high.low=="high",
                                             num.sets(which(get(covariate) >= limit[r])),
                                             num.sets(which(get(covariate) <= limit[r]))),
                                             by=year_of_effect]$V1
    }
    
    if(freq.dur=="duration"){result[,r]<-edLoess[river==rivers[r]&month(date) %in% month,
                                             ifelse(high.low=="high",
                                                    length(which(get(covariate) >= limit[r])),
                                                    length(which(get(covariate) <= limit[r]))),
                                             by=year_of_effect]$V1
    }
  }} else{ #if FUN isn't nothing
    
    for(r in 1:4){
      result[,r]<-edLoess[river == rivers[r] & month(date) %in% month,
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


covariates<-array(dim=c(length(unique(year(edLoess$date)))-1,
                        4,
                        length(names(covariate_inputs))))

dimnames(covariates)<-list(unique(year(edLoess$date))[-1],
                           unique(edLoess$rivers),
                           names(covariate_inputs))

for(i in 1:length(names(covariate_inputs))){
  covariates[,,i]<-do.call(env.cov,covariate_inputs[[i]])
}

assign('covariatesLoess',covariates,env=shared_data)
saveRDS(covariates,"~/process-data/data_store/processed_data/covariates.rds")


