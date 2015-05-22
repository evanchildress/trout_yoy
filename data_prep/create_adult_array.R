bktAdult<-readRDS(file.path("~/westbrookJS/results/nAdult.rds"))

bktAdult<-acast(melt(bktAdult[season==3,list(year,river,mean)],
                      id.vars=c("river","year")),
                 year~river)
bktAdult[1,1:3]<-apply(bktAdult[2:nrow(bktAdult),1:3],2,mean)

bntAdult<-readRDS(file.path("~/westbrookJS/resultsBnt/nAdult.rds"))

bntAdult<-acast(melt(bntAdult[season==3,list(year,river,mean)],
                      id.vars=c("river","year")),
                 year~river)
bntAdult[1,1:3]<-apply(bntAdult[2:nrow(bntAdult),1:3],2,mean)

A<-array(c(bktAdult,bntAdult),dim=c(dim(bktAdult),2))
