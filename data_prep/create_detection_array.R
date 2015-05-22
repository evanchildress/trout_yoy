# load("~/westbrookJS/outMSRiver.rDATA") #results of simpleCJS model to provide detection probabilities
# rm(d)
# detection<-apply(out$BUGSoutput$sims.list$pBeta,
#                  c(2,3,4),
#                  FUN=function(x){
#                    mean(inv.logit(x))})
# detection<-detection[,2:14,1:4]
# 
# saveRDS(detection,file.path(data_root,"detectionFromCJS.rds"))

bktDetect<-readRDS(file.path("~/westbrookJS/results/pYoy.rds"))
bktDetect<-acast(melt(bktDetect[season==3,list(year,river,mean)],
                      id.vars=c("river","year")),
                 year~river)
bktDetect<-bktDetect[2:14,]

bntDetect<-readRDS(file.path("~/westbrookJS/resultsBnt/pYoy.rds"))
bntDetect<-acast(melt(bntDetect[season==3,list(year,river,mean)],
                      id.vars=c("river","year")),
                 year~river)
bntDetect<-bntDetect[2:14,]

detection<-array(c(bktDetect,bntDetect),
                 c(nrow(bktDetect),4,2))

propTagged<-propTagged[3:15,,]

detection<-detection/propTagged #correct for capture of untagged fish (that could be tagged later)
detection[13,,]<-apply(detection[1:12,,],c(2,3),mean)
detection[detection>1]<-1#bnt 2006 is 1.01 after adjustment, so this sets it to 1

#saveRDS(detection,file.path(data_root,'detectionFromJS.rds'))
assign('detection',detection,env=shared_data)
