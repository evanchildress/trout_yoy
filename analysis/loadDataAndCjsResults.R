require(tidyr)
require(getWBData)
require(data.table)
require(dplyr)

#get the raw data from the server
trout<-createCoreData(includeUntagged=T) %>%
       addTagProperties() %>%
       filter(species==whichSpecies &
                area %in% c("inside","trib")&
                detectionDate>as.POSIXct("2002-01-01")) %>%
       addSampleProperties() %>%
       addEnvironmental(sampleFlow=T) %>%
       mutate(stage=as.numeric(year+season/4-0.25-cohort>1)+1) %>%
       mutate(riverName=river) %>%
       mutate(river=as.numeric(factor(river,levels=c('west brook','wb jimmy','wb mitchell','wb obear'),ordered=T))) %>%
       data.table()
nRivers<-length(unique(trout$river))

#make it into a count array
countArray<-trout[stage==1&season==3,.N,by=list(species,year,river)] %>%
  melt(id.vars=c("species","year","river")) %>%
  acast(year~river~species)
countArray[is.na(countArray)]<-0


#get detection probabilities ready
p<-trout[stage==1&season==3,.(meanLength=mean(observedLength,na.rm=T)),by=.(river,year,species)] %>% data.frame()
bktStds<-readRDS(paste0("cjsInputs/",whichSpecies,"YoyStandards.rds") )
bktStds$length<-bktStds$length%>%
  mutate(river=as.numeric(factor(river,levels=c('west brook','wb jimmy','wb mitchell','wb obear'),ordered=T))) %>%
  rename(stdLength=meanLength) %>%
  mutate(species=whichSpecies)
p <- left_join(p,bktStds$length,by=c("river","species")) %>%
              mutate(meanLengthStd=(meanLength-stdLength)/sdLength) %>%
              select(river,year,species,meanLengthStd)
p<- p %>% 
   left_join(trout[stage==1&season==3,.(flowForP=mean(flowForP)),by=.(river,year)] %>%
               data.frame() %>%
               mutate(flowForP=(flowForP-bktStds$flowForP$meanFlow)/bktStds$flowForP$sdFlow)
             ,by=c("river","year"))

p<-readRDS(paste0("cjsInputs/",whichSpecies,"YoyP.rds")) %>% 
  select(mean,river,beta) %>% 
  spread(beta,mean) %>%
  mutate(river=as.numeric(river)) %>%
  right_join(p,by='river') %>%
  data.table() %>%
  setnames(c("1","2","3"),paste0("beta",1:3)) %>%
  mutate(p=1/(1+exp(-(beta1+beta2*meanLengthStd+beta3*flowForP)))) %>%
  filter(species==whichSpecies) %>%
  select(year,river,p) %>%
  melt(id.vars=c("year","river")) %>%
  acast(year~river)

adults<-readRDS(paste0("cjsInputs/",whichSpecies,"Alive.rds")) %>%
  .[,":="(year=as.numeric(year),
       river=as.numeric(river),
       season=as.numeric(season),
       stage=as.numeric(stage))] %>%
  .[stage==2]

fallAdults<-adults[season==3&year %in% 2:14,.(mean,year,river)] %>%
              melt(id.vars=c("year","river")) %>%
              acast(year~river)
wb<-fallAdults[3:13,1]
for(r in 2:nRivers){
  fit<-lm(fallAdults[3:13,r]~wb)
  fallAdults[1,r]<-predict(fit,data.frame(wb=fallAdults[1,1]))
}

#it looks like the first year adult estimates in the tribs might not be so good,
# this estimates those from wb
# for(r in 2:4){
#   fit<-lm(fallAdults[3:13,r]~wb)
#   fallAdults[2,r]<-predict(fit,data.frame(wb=fallAdults[2,1]))
# }

#combine stock from connected rivers
fallAdults[,1:3]<-rowSums(fallAdults[,1:3])


summerAdults<-adults[season==2&year %in% 3:15,.(mean,year,river)] %>%
  melt(id.vars=c("year","river")) %>%
  acast(year~river)

wb<-summerAdults[2:13,1]
for(r in 2:nRivers){
  fit<-lm(summerAdults[2:13,r]~wb)
  summerAdults[1,r]<-predict(fit,data.frame(wb=summerAdults[1,1]))
}


years<-as.numeric(dimnames(countArray)[[1]])

#make jagsData
jagsData<-list(yDATA=countArray[,,1],
               adultDATA=fallAdults,
               summerAdultDATA=summerAdults,
               p=p,
               nRivers=dim(countArray)[2],
               nYears=dim(countArray)[1])
