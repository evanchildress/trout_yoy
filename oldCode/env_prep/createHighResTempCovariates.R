summerT<-highResEnv[month(dateTime) %in% 6:9,
                       length(which(temp>=20))/length(temp),
                       by=list(river,year(dateTime))]

summerT<-acast(melt(summerT,id.vars=c("river","year")),year~river)

meanSummerT<-highResEnv[month(dateTime) %in% 6:8,
                       mean(temp),
                       by=list(river,year(dateTime))]

meanSummerT<-acast(melt(meanSummerT,id.vars=c("river","year")),year~river)

dailyDischarge[,threshold:=quantile(discharge,0.02,na.rm=T),by=river]
lowFlowDays<-dailyDischarge[discharge<threshold,date]
summerInt<-highResEnv[month(dateTime) %in% 6:8,
                      sum(temp>=20 & as.Date(dateTime) %in% lowFlowDays)/length(temp),
                      by=list(river,year(dateTime))]
summerInt<-acast(melt(summerInt,id.vars=c("river","year")),year~river)

assign("summerT",summerT,env=shared_data)
assign("meanSummerT",meanSummerT,env=shared_data)
assign("summerInt",summerInt,env=shared_data)