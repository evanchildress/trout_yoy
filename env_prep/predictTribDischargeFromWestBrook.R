tribFile<-"~/process-data/data_store/original_data/tribDischarge.xlsx"
measuredDischarge<-data.table(read_excel(tribFile))
measuredDischarge<-measuredDischarge[,list(river,date,discharge,stake)]
measuredDischarge[river=='mitchell',discharge:=discharge+0.001]
measuredDischarge[,date:=as.Date(date)]

k<-c(10,20,-3,0)
dailyDischarge[river=='mitchell'&log(depth+20)>3.8,depth:=NA]
dailyDischarge[,logDepth:=
                 log(depth+k[match(river,
                                  unique(dailyDischarge$river))])]

setkey(dailyDischarge,date)
wb<-dailyDischarge[river=='wb',list(date,discharge)]
setnames(wb,c('date','wbDischarge'))
setkey(wb,date)
dailyDischarge<-wb[dailyDischarge]
dailyDischarge[,logWbDischarge:=log(wbDischarge)]

fillDepth<-function(riverName){
  data<-dailyDischarge[!is.na(logDepth)&river==riverName]
  
  if(riverName=="mitchell"){
    fit1<-lm(logDepth~logWbDischarge,data=data[logWbDischarge< -4])
    fit2<-lm(logDepth~logWbDischarge,data=data[logWbDischarge> -4])
    plot(logDepth~logWbDischarge,data=data)
    abline(fit1)
    abline(fit2)
    par(new=T)
    plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
    text(0.2,0.9,bquote(fit1~R^2~"="~.(round(summary(fit1)$r.squared,2))))
    text(0.2,0.8,bquote(fit2~R^2~"="~.(round(summary(fit2)$r.squared,2))))
    text(0.2,0.7,riverName)
    
    dailyDischarge[river==riverName&logWbDischarge< -4,
                   logDepth:=predict(fit1,
                              data.frame(logWbDischarge=logWbDischarge))] 
    dailyDischarge[river==riverName&logWbDischarge> -4,
               logDepth:=predict(fit2,
                          data.frame(logWbDischarge=logWbDischarge))]
  } else {
  
  fit<-lm(logDepth~logWbDischarge,data=data)
  plot(logDepth~logWbDischarge,data=data)
  par(new=T)
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  text(0.2,0.8,bquote(R^2~"="~.(round(summary(fit)$r.squared,2))))
  text(0.2,0.7,riverName)
  
  if(riverName=='obear'){
    dailyDischarge[river==riverName,
               logDepth:=predict(fit,
                          data.frame(logWbDischarge=logWbDischarge))]
  } else{
  dailyDischarge[is.na(logDepth)&river==riverName,
                 logDepth:=predict(fit,
                            data.frame(logWbDischarge=logWbDischarge))]
  }
  }
}

for(r in c("jimmy","mitchell","obear")){
  fillDepth(r)
}

predictDischarge<-function(riverName){
 q<-measuredDischarge[river==riverName]
 setkey(q,date)
 
 depth<-dailyDischarge[river==riverName,list(date,logDepth)]
 setkey(depth,date)

 data<-depth[q]
 
 fit<-lm(log(discharge)~logDepth,data=data)
 plot(log(discharge)~logDepth,data=data)
 abline(fit)
 par(new=T)
 plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
 text(0.2,0.8,bquote(R^2~"="~.(round(summary(fit)$r.squared,2))))
 text(0.2,0.7,riverName)

 
 dailyDischarge[river==riverName&!is.na(logDepth),
                discharge:=exp(predict(fit,
                                  data.frame(logDepth=logDepth)))]
}

for(r in c("jimmy","mitchell","obear")){
  predictDischarge(r)
}

dailyDischarge<-dailyDischarge[,list(date,river,discharge)]

dbWriteTable(link$conn,"data_predicted_discharge",dailyDischarge,
             row.names=FALSE, overwrite=TRUE, append=FALSE)
assign("dailyDischarge",dailyDischarge,env=shared_data)