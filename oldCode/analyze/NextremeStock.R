rivers<-c("West Brook","Jimmy","Mitchell","Obear")

ylimit<-c(1500,800,800,800)
plotN<-function(modelFiles,error=T,
                colors=c('black',rgb(1,0,0),rgb(0,1,0),rgb(0,0,1))){
    #load model output
  getSummary<-function(x){
    mean<-mean(x)
    quants<-quantile(x,c(0.025,0.975))
    results<-c(mean,quants)
    names(results)<-c("mean","lower","upper")
    return(results)
  }
  
  for(m in modelFiles){
    assign(paste0("n",which(m==modelFiles)),
           readRDS(paste0("results/modelOutput/",m))$BUGSoutput$sims.list$nPredict %>%
             apply(c(2,3),getSummary)
           )
  }
  for(r in 1:4){
    plot(NA,ylim=c(0,ylimit[r]),xlim=c(2001,2014),
         ylab="YOY",xlab=NA,main=rivers[r])
    for(m in 1:length(modelFiles)){
      points(get(paste0("n",m))[1,,r]~I(2002:2014),type='l',col=colors[m])
      error.bar(x=2002:2014,
                y=get(paste0("n",m))[1,,r],
                upper.y=get(paste0("n",m))[3,,r],
                lower.y=get(paste0("n",m))[2,,r],
                interval.type="absolute")
    }
  }
}

models<-c("randomYearBkt.rds","StockRecruitBkt.rds","meanBkt.rds","meanStockRecruitBkt.rds")
  
plotN(models)

tiff.par("~/trout_yoy/results/figures/manuscript/N.tif",
         mfcol=c(4,2),mgp=c(2.5,0.5,0),mar=c(1.5,3.5,1,0))
plotN(models)
dev.off()
#         legend(0.35,0.65,c("Random Year","Extreme + Stock Recruit","Extreme","Stock Recruit"),lty=1,
#                col=palette()[1:4],bty='n')
