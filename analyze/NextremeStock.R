fileDir<-'~/trout_yoy/results/modelOutput/speciesRiverSeparate'
setwd(fileDir)

riverNames<-c("Jimmy","Mitchell","Obear","West Brook")
rivers<-riverNames<-c("Jimmy","Mitchell","Obear","WestBrook")

load('~/trout_yoy/abundanceData.RData')
maxAdult<-apply(A,c(2,3),max)


#plots predicted stock recruit relationship with gray lines as
#iterations to indicate uncertainty
plotN<-function(species,river,error=F){
    #load model output
  r<-which(river==c("Jimmy","Mitchell","Obear","WestBrook"))
  
  colors<-c('black',rgb(1,0,0),rgb(0,1,0),rgb(0,0,1))
  
  getSummary<-function(modelType){
    data<-readRDS(paste0(modelType,species,river,".rds"))$BUGSoutput$sims.list
    
    summarize<-function(x){
      meanX<-mean(x)
      quantX<-quantile(x,c(0.025,0.975))
      
      return(c(meanX,quantX))
    }
    
    summary<-data.frame(t(apply(data$nPredict,2,summarize)))
    names(summary)<-c("n","lower","upper")
    
    return(summary)
  }
  ylimit<-c(0,c(320,400,600,1500)[r])
  
  plot(NA,ylim=ylimit,xlim=c(2001,2014),
       ylab="YOY",xlab=NA,main=riverNames[r])
  
  lineType<-1
  for(model in c("randomYear","extremeStockRecruit",
                 "extreme","stockRecruit")){
    data<-getSummary(model)
    points(n~c(2002:2014),type='l',col=palette()[lineType],
           lty=1,data=data,lwd=1)
        points(n~c(2002:2014),col=palette()[lineType],
           data=data)
    with(data,error.bar(c(2002:2014),n,upper.y=upper,lower.y=lower,
                        interval.type='absolute'))
    lineType<-lineType+1
  }  

}

tiff.par("~/trout_yoy/results/figures/manuscript/NextremeStock.tif",
         mfcol=c(4,2),mgp=c(2.5,0.5,0),mar=c(1.5,3.5,1,0))
  for(sp in c("Bkt","Bnt")){
    for(r in rivers){
      if(sp=='Bnt'&r=='Mitchell'){
        plot(NA,xlim=c(0,1),ylim=c(0,1),
             axes=F,xlab="",ylab="")
        legend(0.35,0.65,c("Random Year","Extreme + Stock Recruit","Extreme","Stock Recruit"),lty=1,
               col=palette()[1:4],bty='n')
        next
      }
      if(sp=='Bnt' & r=='Obear'){
             plot(NA,xlim=c(0,1),ylim=c(0,1),
             axes=F,xlab="",ylab="")
        next}
      plotN(species=sp,river=r)
    }
  }
dev.off()