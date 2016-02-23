plotCovariate<-function(covariate,fileName){
  #coefficient names, which are the options for covariate
  coefNames<-c("fallQ","winterQ","springQ","summerQ","summerT")
  betaNum<-which(covariate==coefNames)
  
  #load model output
  filePath<-file.path("~/trout_yoy/results/modelOutput/speciesRiverSeparate",
                      fileName)
  data<-readRDS(filePath)$BUGSoutput$sims.list
  
  #identify the species, river, and model type its from
  sp<-ifelse(grepl('Bkt',fileName),1,2)
  species<-c("Bkt","Bnt")[sp]

  
  river<-strsplit(strsplit(fileName,species)[[1]][2],".rds")[[1]]
  r<-which(river==c("Jimmy","Mitchell","Obear","WestBrook"))
  

  
  if(grep("tockRecruit",fileName)){
      load('~/trout_yoy/abundanceData.RData')
      if(r==3){meanAdult<-mean(A[,3,sp])
      } else {
          meanAdult<-mean(apply(A[,c(1,2,4),sp],1,sum))
      }
      
      lambdaStock<-mean(data$c[,1])+mean(data$c[,2]*meanAdult)
  } else {
    lambdaStock<-0
    meanAdult<-1
  }
 
  #lambdaStock<-0 #if stock is centered then 0 is the mean, otherwise, it should be calculated above
  
  modelType<-strsplit(fileName,species)[[1]][1]
    
  x<-seq(-2,2,0.01)
  
  ylimits<-list(Bkt=c(1000,600,1500,1500),
                Bnt=c(500,10,10,1500))
  
  plot(NA,xlim=c(-2,2),ylim=c(0,ylimits[[species]][r]),
       xlab="",ylab="Recruits",main=covariate)

    numIter<-nrow(data$beta)
    numIter<-300 #to speed it up on test runs
    
    colors<-c(rgb(0,0,1,0.1),rgb(1,0,0,0.1),rgb(0,1,0,0.1),rgb(1,0,0,0.1))
    #colors<-c(rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
    #colors<-c('blue','green','red')
    
    if(covariate %in% c("fallQ","winterQ")){  
      coefChain<-cbind(data$beta[,betaNum],data$beta[,betaNum+6])
      for(n in 1:numIter){
        y1<-exp(x*coefChain[n,1]+x^2*coefChain[n,2]+lambdaStock)*meanAdult
        points(y1~x,type='l',col=colors[4])
        }
    } else if(covariate %in% c("summerT","summerQ")){
        coefChain<-cbind(data$beta[,betaNum],data$beta[,6])
        x2<-c(-1,-0.5,0,1)
        if(river=='WestBrook'){sampleFrom<-1:2
        } else {sampleFrom<-c(1,3,4)}
        
          for(n in 1:numIter){
            i<-sample(sampleFrom,1)
            y1<-exp(x*coefChain[n,1]
                    +coefChain[n,2]*x*x2[i]
                    +lambdaStock)*meanAdult
            points(y1~x,type='l',col=colors[i])
          }
    } else {
        coefChain<-data$beta[,betaNum]
        for(n in 1:numIter){
          y1<-exp(x*coefChain[n]+lambdaStock)*meanAdult
          points(y1~x,type='l',col=colors[4])
          }
    }
}

rivers<-c("Jimmy","Mitchell","Obear","WestBrook")
coefNames<-c("fallQ","winterQ","springQ","summerQ","summerT")

for(sp in c("Bkt","Bnt")){
tiff.par(paste0("~/trout_yoy/results/figures/manuscript/",
                "envCovariates",sp,".tif"),
         mfrow=c(length(coefNames),length(rivers)),mar=c(1.5,3.5,1,0),
         cex.main=0.7)
  for(cov in coefNames){
    for(r in rivers){
        if(r %in% c("Mitchell","Obear") & sp=="Bnt"){
          plot(NA,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F)
          next
        }
      
      fileName<-paste0("extremeStockRecruit",sp,r,".rds")
      plotCovariate(cov,fileName)
    }
  }
  dev.off()
}

