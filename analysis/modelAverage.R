library(MuMIn)
library(ggplot2)
fallSize<-trout[stage==2&season==3,
                .(meanLength=mean(observedLength,na.rm=T)),
                by=.(river,year)] %>%
  melt(id.vars=c("river","year")) %>%
  acast(year~river)
fallSize[is.na(fallSize)]<-mean(fallSize[,3],na.rm=T)

summAdults<-jagsData$summerAdultDATA %>%
            apply(2,function(x){scale(x)[,1]}) %>%
            melt() %>%
            data.table() %>%
            setnames(c("year","river","summerAdults")) %>%
            .[,year:=year+1999] %>%
            .[,river:=c("west brook","wb jimmy","wb mitchell","wb obear")[river]] %>%
            setkey(river,year)


eggs<-(0.00187*fallSize^2.19)*fallAdults
eggs[,1:3]<-rowSums(eggs[,1:3])
yoy<-jagsData$yDATA/jagsData$p

yoySurv<-(yoy/eggs) %>%
  melt() %>%
  data.table() %>%
  setnames(c("year","river","yoySurv")) %>%
  .[,river:=c("west brook","wb jimmy","wb mitchell","wb obear")[river]] %>%
  .[,logitYoySurv:=qlogis(yoySurv)] %>%
  setkey(river,year)

t<-t[,month:=month(datetime)] %>%
  .[,yearOfEffect:=ifelse(month(datetime) %in% c(10,11,12),year(datetime)+1,year(datetime))] %>%
  .[,season:=as.numeric(month>=1)+as.numeric(month>=4)+as.numeric(month>=7)+as.numeric(month>=10)] %>%
  .[,performance:=predictPerformance(temperature,13,19.5,4)]

temps<-t[,.(meanTemp=mean(temperature,na.rm=T)),.(river,yearOfEffect,season)] %>%
  spread(season,meanTemp) %>%
  setnames(c("1","2","3","4"),paste0("meanTemp",1:4)) %>%
  .[,":="(meanTemp1=scale(meanTemp1)[,1],
          meanTemp2=scale(meanTemp2)[,1],
          meanTemp3=scale(meanTemp3)[,1],
          meanTemp4=scale(meanTemp4)[,1]),
    by=river] %>%
  setkey(river,yearOfEffect)

perf<-t[,.(meanPerformance=mean(performance,na.rm=T)),.(river,yearOfEffect,season)] %>%
  spread(season,meanPerformance) %>%
  setnames(c("1","2","3","4"),paste0("meanPerf",1:4)) %>%
  .[,":="(meanPerf1=scale(meanPerf1)[,1],
          meanPerf2=scale(meanPerf2)[,1],
          meanPerf3=scale(meanPerf3)[,1],
          meanPerf4=scale(meanPerf4)[,1]),
    by=river] %>%
  setkey(river,yearOfEffect)

q<-q[,month:=month(date)] %>%
  .[,yearOfEffect:=ifelse(month(date) %in% c(10,11,12),year(date)+1,year(date))] %>%
  .[,season:=as.numeric(month>=1)+as.numeric(month>=4)+as.numeric(month>=7)+as.numeric(month>=10)]

discharge<-q[,.(meanDischarge=mean(discharge,na.rm=T)),.(river,yearOfEffect,season)] %>%
  spread(season,meanDischarge) %>%
  setnames(c("1","2","3","4"),paste0("meanDischarge",1:4)) %>%
  .[,":="(meanDischarge1=scale(meanDischarge1)[,1],
          meanDischarge2=scale(meanDischarge2)[,1],
          meanDischarge3=scale(meanDischarge3)[,1],
          meanDischarge4=scale(meanDischarge4)[,1]),
    by=river] %>%
  setkey(river,yearOfEffect)

yoySurv<-summAdults[perf[discharge[temps[yoySurv]]]]

a<-lm(logitYoySurv~river + meanTemp1*river + meanTemp2*river + meanTemp3*river + 
  meanDischarge1*river + meanDischarge2*river + meanDischarge3*river + 
    meanDischarge4*river + summerAdults*river,
  data=yoySurv,
  na.action=na.fail)

d<-dredge(a)
m<-model.avg(d[d$delta<4],fit=T)
saveRDS(m,"results/yoySurvivalModel.rds")

a<-lm(logitYoySurv~river + meanPerf3 + meanDischarge1*river + meanDischarge2*river + 
        meanDischarge3*river + meanDischarge4*river + summerAdults*river,data=yoySurv)

ggplot(data=yoySurv,aes(x=get(paste0("meanDischarge",2)),y=yoySurv)) +
  facet_wrap(~river,scales="free_x") +
  geom_point()

