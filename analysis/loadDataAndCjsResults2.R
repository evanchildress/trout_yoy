
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
meanLength<-trout[stage==1&season==3,.(meanLength=mean(observedLength,na.rm=T)),by=.(river,year,species)] %>%
    setkey(year) %>%
    .[,.SD[data.table(year=2002:2015,key="year")],by=.(river,species)] %>%
    .[,riverMean:=mean(meanLength,na.rm=T),by=river] %>%
    .[is.na(meanLength),meanLength:=riverMean] %>%
    .[,":="(riverMean=NULL,
            species=NULL)] %>%
    data.table() %>%
    setkey(river)

bktStds<-readRDS(paste0("cjsInputs/",whichSpecies,"YoyStandards.rds") )

meanLength<-bktStds$length %>%
            data.table() %>%
            setnames("meanLength","stdLength") %>%
            .[,river:=c(2,3,4,1)] %>%
            setkey(river) %>%
            .[meanLength] %>%
            .[,meanLength:=(meanLength-stdLength)/sdLength] %>%
            .[,":="(stdLength=NULL,
                    sdLength=NULL)] %>%
            melt(id.vars=c("river","year")) %>%
            acast(year~river)
flowForP<-trout[stage==1&season==3,.(flowForP=mean(flowForP)),by=.(river,year)] %>%
               .[,flowForP:=(flowForP-bktStds$flowForP$meanFlow)/bktStds$flowForP$sdFlow] %>%
                melt(id.vars=c("river","year")) %>%
                acast(year~river) 

pPriorMean<-readRDS(paste0("cjsInputs/",whichSpecies,"YoyPPosterior.rds")) %>% 
  .[,.(mean,river,beta)] %>%
  melt(id.vars=c("river","beta")) %>%
  acast(river~beta)

pPriorSd<-readRDS(paste0("cjsInputs/",whichSpecies,"YoyPPosterior.rds")) %>% 
  .[,.(sd,river,beta)] %>%
  melt(id.vars=c("river","beta")) %>%
  acast(river~beta)

adults<-readRDS(paste0("cjsInputs/",whichSpecies,"Alive.rds")) %>%
  .[,":="(year=as.numeric(year),
       river=as.numeric(river),
       season=as.numeric(season),
       stage=as.numeric(stage))] %>%
  .[stage==2]

fallAdults<-adults[season==3&year %in% 2:15,.(mean,year,river)] %>%
              melt(id.vars=c("year","river")) %>%
              acast(year~river)
wb<-fallAdults[3:13,1]
for(r in 2:nRivers){
  fit<-lm(fallAdults[3:13,r]~wb)
  fallAdults[1,r]<-predict(fit,data.frame(wb=fallAdults[1,1]))
}

# it looks like the first year adult estimates in the tribs might not be so good,
# this estimates those from wb
for(r in 2:nRivers){
  fit<-lm(fallAdults[3:13,r]~wb)
  fallAdults[2,r]<-predict(fit,data.frame(wb=fallAdults[2,1]))
}

# #combine stock from connected rivers
# fallAdults[,1:3]<-rowSums(fallAdults[,1:3])


summerAdults<-adults[season==2&year %in% 3:16,.(mean,year,river)] %>%
  melt(id.vars=c("year","river")) %>%
  acast(year~river)

wb<-summerAdults[2:14,1]
for(r in 2:nRivers){
  fit<-lm(summerAdults[2:14,r]~wb)
  summerAdults[1,r]<-predict(fit,data.frame(wb=summerAdults[1,1]))
}

otherSpecies<-c("bnt","bkt")[which(whichSpecies!=c("bnt","bkt"))]
otherSpeciesAdults<-readRDS(paste0("cjsInputs/",otherSpecies,"Alive.rds")) %>%
  .[,":="(year=as.numeric(year),
          river=as.numeric(river),
          season=as.numeric(season),
          stage=as.numeric(stage))] %>%
  .[stage==2]

#it looks like the first year adult estimates in the tribs might not be so good,
# this estimates those from wb
# for(r in 2:4){
#   fit<-lm(fallAdults[3:13,r]~wb)
#   fallAdults[2,r]<-predict(fit,data.frame(wb=fallAdults[2,1]))
# }

otherSpeciesSummerAdults<-otherSpeciesAdults[season==2&year %in% 3:16,.(mean,year,river)] %>%
  melt(id.vars=c("year","river")) %>%
  acast(year~river)

wb<-otherSpeciesSummerAdults[2:14,1]
for(r in 2:3){
  fit<-lm(otherSpeciesSummerAdults[2:14,r]~wb)
  otherSpeciesSummerAdults[1,r]<-predict(fit,data.frame(wb=otherSpeciesSummerAdults[1,1]))
}



years<-as.numeric(dimnames(countArray)[[1]])

#make jagsData
jagsData<-list(yDATA=countArray[,,1],
               adultDATA=fallAdults,
               summerAdultDATA=summerAdults,
               pPriorMean=pPriorMean,
               pPriorSd=pPriorSd,
               meanLength=meanLength,
               flowForP=flowForP,
               nRivers=dim(countArray)[2],
               nYears=dim(countArray)[1],
               otherSpeciesSummerAdultDATA=otherSpeciesSummerAdults)
