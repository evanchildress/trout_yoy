plotCovariate<-function(file,meanAdult=apply(fallAdults,2,mean)){
#   #coefficient names, which are the options for covariate
#   coefNames<-c("fallQ","winterQ","springQ","summerQ","summerT")
#   betaNum<-which(covariate==coefNames)
#   
#   #load model output
#   filePath<-file.path("~/trout_yoy/results/modelOutput/speciesRiverSeparate",
#                       fileName)
  data<-readRDS(paste0("results/modelOutput/",file))$BUGSoutput$sims.list
  
  getSig<-function(x){
    quants<-quantile(x,c(0.025,0.975))
    sig<-all(quants>0)|all(quants<0)
    return(sig)
  }
  
  sig<-apply(data$beta,c(2,3),getSig)
  
  x<-seq(-2,2,0.1)
  predicted<-array(NA,dim=c(dim(data$beta)[1],length(x),8,4))

  for(i in 1:dim(data$beta)[1]){
    for(r in 1:4){
      predicted[i,,1,r]<-exp(x*data$beta[i,1,r]+x^2*data$beta[i,7,r])
      predicted[i,,2,r]<-exp(x*data$beta[i,2,r])
      predicted[i,,3,r]<-exp(x*data$beta[i,3,r])
      predicted[i,,4,r]<-exp(x*data$beta[i,4,r])
      predicted[i,,5,r]<-exp(x*data$beta[i,5,r])
      predicted[i,,6,r]<-exp(x*data$beta[i,8,r])
      predicted[i,,7,r]<-exp(x*data$beta[i,6,r]*-1+x*data$beta[i,5,r])
      predicted[i,,8,r]<-exp(x*data$beta[i,6,r]*1+x*data$beta[i,5,r])
    }
  }
  betaNum<-c(1,2,3,4,5,8)
  tiff.par("results/figures/environmentalBetas.tif",mfrow=c(6,4),mar=c(0.5,0.5,0,0))
    for(b in 1:6){
      for(r in 1:4){
        plot(NA,ylim=c(0,10),xlim=range(x))
        for(i in 1:dim(predicted)[1]){
          points((predicted[i,,b,r])~x,col=ifelse(sig[betaNum[b],r],'black','gray'),type='l')
        }
      }
    }
  dev.off()
  
#   #identify the species, river, and model type its from
#   sp<-ifelse(grepl('Bkt',fileName),1,2)
#   species<-c("Bkt","Bnt")[sp]
# 
#   
#   river<-strsplit(strsplit(fileName,species)[[1]][2],".rds")[[1]]
#   r<-which(river==c("Jimmy","Mitchell","Obear","WestBrook"))
#   

  
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
  
  ylimits<-list(Bkt=c(10,10,10,10),
                Bnt=c(500,10,10,1500))
  
  plot(NA,xlim=c(-2,2),ylim=c(-1.5,1.5),yaxt='n',
       xlab="",ylab="")
  axis(2,at=seq(-2,2,1),c(0.01,0.1,1,10,100))
    numIter<-nrow(data$beta)
    numIter<-1000 #to speed it up on test runs
    
    colors<-c(rgb(0,0,1,0.1),rgb(1,0,0,0.1),rgb(0,1,0,0.1),rgb(1,0,0,0.1))
    #colors<-c(rgb(0,0,1,0.3),rgb(0,1,0,0.3),rgb(1,0,0,0.3))
    #colors<-c('blue','green','red')
    
    if(covariate %in% c("fallQ","winterQ")){  
      coefChain<-cbind(data$beta[,betaNum],data$beta[,betaNum+6])
      for(n in 1:numIter){
        y1<-exp(x*coefChain[n,1]+x^2*coefChain[n,2]+lambdaStock)*meanAdult
        y1<-y1/(exp(lambdaStock)*meanAdult)
        y1<-log10(y1)
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
            y1<-y1/(exp(lambdaStock)*meanAdult)
            y1<-log10(y1)
            points(y1~x,type='l',col=colors[i])
          }
    } else {
        coefChain<-data$beta[,betaNum]
        for(n in 1:numIter){
          y1<-exp(x*coefChain[n]+lambdaStock)*meanAdult
          y1<-y1/(exp(lambdaStock)*meanAdult)
          y1<-log10(y1)
          points(y1~x,type='l',col=colors[4])
          }
    }
}


rivers<-c("Jimmy","Mitchell","Obear","WestBrook")
coefNames<-c("fallQ","winterQ","springQ","summerQ","summerT")
coefLabels<-c("Fall\nDischarge\n(mean)","Winter\nFloods\n(days > q98)",
              "Spring\nFloods\n(days > q98)",
              "Summer\nDrought\n(days < q02)","Summer\nTemperature\n(% > 20C)")

tiff.par("~/trout_yoy/results/figures/manuscript/envCovariates.tif",
         cex.main=0.7,width=9,height=6.1)

layout(matrix(c(0,1,1,1,1,0,6,6,
                0,2,3,4,5,0,7,8,
                9,10,11,12,13,45,14,15,
                16,17,18,19,20,45,21,22,
                23,24,25,26,27,45,28,29,
                30,31,32,33,34,45,35,36,
                37,38,39,40,41,45,42,43,
                0,44,44,44,44,44,44,44),ncol=8,byrow=T),
       widths=c(1.1,2,2,2,2,0.1,2,2),heights=c(0.3,0.3,2,2,2,2,2,0.3))

  
#   if(sp == "Bkt"){
#     layout(matrix(1:24,ncol=4),widths=rep(2,4),heights=c(0.75,2,2,2,2,2))
#   } else {
#     layout(matrix(1:12,ncol=2),widths=rep(2,2),heights=c(0.75,2,2,2,2,2))
#   }
speciesNames<-c("Brook Trout","Brown Trout")
for(sp in 1:2){
  par(mar=c(0,3.5,0,0))
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  text(0.5,0.5,speciesNames[sp],font=2)

  
  for(r in rivers){
      if(r %in% c("Mitchell","Obear") & sp==2){
        #plot(NA,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F)
        next
      }
      par(mar=c(0,3.5,0,0))
      plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
      text(0.5,0.5,ifelse(r=="WestBrook","West Brook",r))
  }
}



for(cov in coefNames){
  #label for coefficients located in left column
  par(mar=c(0,3.5,0,0),xpd=NA)
  plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
  text(-1.7,0.5,coefLabels[which(cov==coefNames)],adj=c(0,0.5),font=2)
  if(cov=="springQ"){
    text(1.4,0.5,"Recruits (as a proportion of the mean)",srt=90,cex=1.2)
  }
  
  par(mar=c(1.5,3.5,1,0),xpd=F)
  
  for(sp in c("Bkt","Bnt"))
    for(r in rivers){
      if(r %in% c("Mitchell","Obear") & sp=="Bnt"){
        #plot(NA,xlim=c(0,1),ylim=c(0,1),ylab='',xlab='',axes=F)
        next
      }
      fileName<-paste0("extremeStockRecruit",sp,r,".rds")
      plotCovariate(cov,fileName)
    }
}
par(mar=c(0,0,0,0),xpd=NA)
plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
text(0.5,0.5,"Standardized Environmental Covariates",cex=1.2)

plot(NA,xlim=c(0,1),ylim=c(0,1),axes=F,xlab='',ylab='')
abline(v=2)

dev.off()

par(xpd=F)
