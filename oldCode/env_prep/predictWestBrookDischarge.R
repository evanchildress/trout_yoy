library(waterData)

#summaryize high res depth measurements by day and river
#and fill with NAs when measurements are absent
daily<-highResEnv[,mean(depth,na.rm=T),by=list(as.Date(dateTime),river)]
allDates<-data.table(date=rep(seq(as.Date("1997-05-27"),
                              as.Date("2015-05-29"),
                              "day"),each=4),
                     river=c("jimmy","mitchell","obear","wb"))
setnames(daily,c("date","river","depth"))
setkey(daily,date,river)
setkey(allDates,date,river)

daily<-daily[allDates]


#load YSI data for west brook discharge
wbFile<-"~/process-data/data_store/original_data/West Brook YSI Data with discharge calculated.xls"
wb<-data.table(read_excel(wbFile,
                           col_types=c("numeric","numeric",
                                       "numeric","numeric","numeric",
                                       "numeric","numeric","numeric",
                                       "numeric","numeric","numeric",
                                       "numeric")))
setnames(wb,tolower(names(wb)))
wb<-wb[,list(date,time,discharge,depth)]
wb[,date:=as.POSIXct(date*60*60*24,origin="1899-12-30",tz="UTC")]
wb<-wb[,list(discharge=mean(discharge)),by=date]
wb<-wb[,list(date,discharge)]
wb[,date:=as.Date(date)]
setkey(wb,date)
allDates<-allDates[river=='wb',list(date)]
wb<-wb[allDates]

#Function imports discharge data from USGS database online
cleanQImport<-function(riverCode){
  q<-data.table(importDVs(riverCode,code="00060",stat="00003",sdate="1997-05-27",
               edate="2015-05-29"))
  setnames(q,c('val','dates'),c('discharge','date'))
  q[discharge<0,discharge:=NA]
  q<-q[,list(date,discharge)]

  setkey(q,date)

  return(q)
}
#import these rivers
riverCodes<-c("01169900",
              "01171500",
              "01161000")
k<-c(10,5,0)
#run the imports and merge with YSI data
for(code in riverCodes){
  q<-cleanQImport(code)
  q<-q[,list(date,log(discharge+k[which(code==riverCodes)]))]
  setnames(q,c("date",paste0("q",which(code==riverCodes))))
  setkey(q,date)
  wb<-q[wb]
}
wb<-wb[!is.na(q1)&!is.na(q2)&!is.na(q3)]

wbLm<-lm(log(discharge)~q1*q2*q3,data=wb)

wb[is.na(discharge),
   discharge:=exp(predict(wbLm,data.frame(q1=q1,q2=q2,q3=q3)))]




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
for(r in unique(daily$river)){
  daily[river==r & 
          date %in% badDates[[r]],depth:=NA]
}

wb<-daily[river=='wb'][wb]

wb[,logDepth:=log(depth)]
depthLm<-lm(log(discharge)~logDepth,data=wb)
wb[!is.na(depth),
   discharge:=exp(predict(depthLm,data.frame(logDepth=logDepth)))]
daily[river=='wb'&date %in% wb$date,discharge:=wb$discharge]



assign("dailyDischarge",daily,env=shared_data)