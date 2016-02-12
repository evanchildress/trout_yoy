#get the raw data from the server
trout<-createCoreData(includeUntagged=T) %>%
       addTagProperties() %>%
       filter(species %in% c("bkt","bnt")&
                area %in% c("inside","trib")&
                detectionDate>as.POSIXct("2002-01-01")) %>%
       addSampleProperties() %>%
       addEnvironmental(sampleFlow=T) %>%
       mutate(stage=as.numeric(year+season/4-0.25-cohort>1)+1) %>%
       mutate(riverName=river) %>%
       mutate(river=as.numeric(factor(river,levels=c('west brook','wb jimmy','wb mitchell','wb obear'),ordered=T))) %>%
       data.table()

#make it into a count array
countArray<-trout[!is.na(stage),.N,by=list(species,year,season,river,stage)] %>%
  melt(id.vars=c("species","year","season","river","stage")) %>%
  acast(year~season~river~stage~species)
countArray[is.na(countArray)]<-0


#get detection probabilities ready
p<-trout[stage==1&season==3,.(meanLength=mean(observedLength,na.rm=T)),by=.(river,year,species)] %>% data.frame()
bktStds<-readRDS("cjsInputs/bktYoyStandards.rds") 
bktStds$length<-bktStds$length%>%
  mutate(river=as.numeric(factor(river,levels=c('west brook','wb jimmy','wb mitchell','wb obear'),ordered=T))) %>%
  rename(stdLength=meanLength) %>%
  mutate(species="bkt")
p <- left_join(p,bktStds$length,by=c("river","species")) %>%
              mutate(meanLengthStd=(meanLength-stdLength)/sdLength) %>%
              select(river,year,species,meanLengthStd)
p<- p %>% 
   left_join(trout[stage==1&season==3,.(flowForP=mean(flowForP)),by=.(river,year)] %>%
               data.frame() %>%
               mutate(flowForP=(flowForP-bktStds$flowForP$meanFlow)/bktStds$flowForP$sdFlow)
             ,by=c("river","year"))

p<-readRDS("cjsInputs/bktYoyP.rds") %>% 
  select(mean,river,beta) %>% 
  spread(beta,mean) %>%
  mutate(river=as.numeric(river)) %>%
  right_join(p,by='river') %>%
  data.table() %>%
  setnames(c("1","2","3"),paste0("beta",1:3)) %>%
  mutate(p=1/(1+exp(-(beta1+beta2*meanLengthStd+beta3*flowForP)))) %>%
  filter(species=="bkt") %>%
  select(year,river,p) %>%
  melt(id.vars=c("year","river")) %>%
  acast(year~river)

#make jagsData
jagsData<-list(y=countArray[,,,1,1],
               #adults=abundanceArray[,,,,2],
               p=p)

