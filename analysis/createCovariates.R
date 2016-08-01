reconnect()
# t<-tbl(conDplyr,"data_hourly_temperature") %>% collect(n=Inf) %>% mutate(date=as.Date(datetime)) %>% data.table()
# saveRDS(t,"cjsInputs/tempData.rds")
t<-readRDS("cjsInputs/tempData.rds")
# if(qType=="flowExtension"){
#   q<-tbl(conDplyr,"data_flow_extension") %>% collect() %>% data.table() %>%
#     .[,date:=as.Date(date)]
# # just using the flow extension for all the rivers, but the function requires river specific discharge
#   q<-rbind(q,
#            q %>% mutate(river="wb mitchell"),
#            q %>% mutate(river="wb jimmy"),
#            q %>% mutate(river="wb obear")) %>%
#     setnames("qPredicted","discharge") %>%
#     mutate(date=as.Date(date))} else {
#       
#   q<-tbl(conDplyr,"data_daily_discharge") %>% collect(n=Inf) %>% data.table()
#     }

#saveRDS(q,"cjsInputs/dischargeData.rds")
q<-readRDS("cjsInputs/dischargeData.rds")


# q<-tbl(src_postgres("wb",user="postgres"),"data_predicted_discharge") %>% collect() %>% data.table() %>%
#   .[,river:=c("west brook","wb jimmy","wb mitchell","wb obear")[
#     match(river,c("wb","jimmy","mitchell","obear"))]]


start_year<-min(years)
end_year<-max(years)

rivers<-c("west brook","wb jimmy","wb mitchell","wb obear")

env.cov<-function(covariate,month,threshold=NA,high.low=NA,
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

covariates<-array(dim=c(end_year-start_year+1,4,length(names(covariate_inputs))))
dimnames(covariates)<-list(years,
                           rivers,
                           names(covariate_inputs))


for(i in c(1:length(names(covariate_inputs)))){
  covariates[,,i]<-do.call(env.cov,covariate_inputs[[i]])
}

#temporary while temps aren't available for 2015
# for(i in c(1:4,6:8,10:length(names(covariate_inputs)))){
#   covariates[,,i]<-do.call(env.cov,covariate_inputs[[i]])
# }
# for(i in c(5,9)){
#   end_year<-2014
#   covariates[1:13,,i]<-do.call(env.cov,covariate_inputs[[i]])
#   covariates[14,,i]<-colMeans(covariates[1:13,,i])
# }


jagsData$covariates<-covariates
saveRDS(jagsData,paste0("cjsInputs/jagsData",
                        toupper(substr(whichSpecies,1,1)),
                        substr(whichSpecies,2,nchar(whichSpecies)),
                        ".rds"))
