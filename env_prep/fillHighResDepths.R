library(waterData)

daily<-highResEnv[,mean(depth,na.rm=T),by=list(as.Date(dateTime),river)]
allDates<-data.table(date=rep(seq(as.Date("1997-05-27"),
                              as.Date("2015-05-29"),
                              "day"),each=4),
                     river=c("jimmy","mitchell","obear","wb"))
setnames(daily,c("date","river","depth"))
setkey(daily,date,river)
setkey(allDates,date,river)

daily<-daily[allDates]


cleanQImport<-function(riverCode){
  q<-data.table(importDVs(riverCode,code="00060",stat="00003",sdate="1997-05-27",
               edate="2015-05-29"))
  setnames(q,c('val','dates'),c('discharge','date'))
  q[discharge<0,discharge:=NA]
  q<-q[,list(date,discharge)]

  setkey(q,date)

  return(q)
}

riverCodes<-c("01169900",
              "01171500",
              "01161000")

for(code in riverCodes){
  q<-cleanQImport(code)
  q<-q[,list(date,log(discharge))]
  setnames(q,c("date",paste0("q",which(code==riverCodes))))
  setkey(q,date)
  
  daily<-q[daily]
}


k<-c(10,20,-3,20)

badDates<-list(jimmy=
                 c(seq(as.Date("2009-07-30"),as.Date("2009-12-27"),"day"),
                   seq(as.Date("2009-01-15"),as.Date("2009-04-29"),"day")),
               mitchell=
                 c(seq(as.Date("2013-03-01"),as.Date("2015-05-31"),"day"),
                   seq(as.Date("2009-12-12"),as.Date("2010-01-12"),"day"),
                   seq(as.Date("2010-04-02"),as.Date("2010-07-19"),"day")),
               obear=
                 c(seq(as.Date("2009-07-30"),as.Date("2009-12-27"),"day")),
               wb=
                 c(seq(as.Date("2008-02-26"),as.Date("2010-02-26"),"day"))
               )


predictDepth<-function(riverName){
  daily[date %in% badDates[[riverName]] & river==riverName,depth:=NA]
  
  r<-which(unique(daily$river)==riverName)
  daily[river==riverName,logDepth:=log(depth+k[r])]
  
  fit<-lm(logDepth~q2*q1*q3,
          data=daily[river==riverName])
  
#   plot(logDepth~q2,
#           data=daily[river==riverName])
   
  cat(riverName,"R^2 = ",summary(fit)$r.squared,"\n")

  daily[is.na(logDepth) & river==riverName,
        logDepth:=predict(fit,data.frame(q1=q1,q2=q2,q3=q3))]
}
for(riv in unique(daily$river)){
predictDepth(riv)
}

daily[river=="mitchell",logDepth:=log(exp(logDepth)-19)]
daily[river=="obear",logDepth:=log(exp(logDepth)+100)]

assign("dailyDischarge",daily,env=shared_data)
