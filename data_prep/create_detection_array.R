load("~/simpleCJS/outMSRiver.rDATA") #results of simpleCJS model to provide detection probabilities
rm(d)
detection<-apply(out$BUGSoutput$sims.list$pBeta,
                 c(2,3,4),
                 FUN=function(x){
                   mean(inv.logit(x))})

saveRDS(detection,file.path(data_root,"detectionFromCJS.rds"))