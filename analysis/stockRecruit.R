fileDir<-'~/trout_yoy/results/modelOutput/speciesRiverSeparate'
setwd(fileDir)

riverNames<-c("Jimmy","Mitchell","Obear","West Brook")

load('~/trout_yoy/abundanceData.RData')


nBetas<-8
betaNames<-c("Fall Q",
             "Winter Q",
             "Spring Q",
             "Summer Q",
             "Summer T",
             "Summer QxT",
             expression("Fall "~Q^2),
             
             expression("Winter "~Q^2))

coefNames<-c(paste0('beta',1:nBetas),paste0('c',1:3))
for(sp in 1:2){
A[,c(1,2,4),sp]<-apply(A[,c(1,2,4),sp],1,sum)}

maxAdult<-apply(A,c(2,3),max)

#plots predicted stock recruit relationship with gray lines as
#iterations to indicate uncertainty
plotStockRecruit<-function(fileName,fun='plot',ylim=c(0,1100),
                           error=F,perCapita=F,...){
    #load model output
  data<-readRDS(fileName)$BUGSoutput$sims.list
  
  #identify the species, river, and model type its from
  sp<-ifelse(grepl('Bkt',fileName),1,2)
  species<-c("Bkt","Bnt")[sp]
  
  river<-strsplit(strsplit(fileName,species)[[1]][2],".rds")[[1]]
  r<-which(river==c("Jimmy","Mitchell","Obear","WestBrook"))
  
  modelType<-strsplit(fileName,species)[[1]][1]
  
  colors<-c('black',rgb(1,0,0),rgb(0,1,0),rgb(0,0,1))
  
  if(modelType=="extremeStockRecruit"){lty<-1
  } else {lty<-2}
  
  if(species=='Bnt'& river %in% c('Obear','Mitchell')){
      return(list(id=c(river,species,modelType),
              mean=NA,
              lower=NA,
              upper=NA))
  }
  
  coef<-apply(data$c,2,mean)
  x<-seq(0,maxAdult[r,sp]+1,length.out=100)
  xCentered<-x#(x-mean(A[,r,sp]))/sd(A[,r,sp])
  y<-exp(coef[1]+xCentered*coef[2])
  if(!perCapita){y<-y*x}
  if(fun=='plot'){plot(y~x,type='l',ylim=ylim,col=colors[r])
  } else {
    points(y~x,type='l',col=palette()[r],lwd=2,lty=lty)}
  if(error==T){
    errorColors<-c(gray(0,0.03),rgb(1,0,0,0.03),
                   rgb(0,1,0,0.03),rgb(0,0,1,0.03))
    #errorColors<-c(gray(0,0.3),rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3))

    #     plotError<-function(coef1,coef2){
#       y<-exp(coef1+x*coef2)
#       if(!perCapita){y<-y*x}
#       points(y~x,col=errorColors[r],type='l')
#     }
#     
#   coefTable<-data.table(data$c)
#   coefTable[,index:=1:nrow(coefTable)]
#   coefTable[,plotError(V1,V2),by=index]
#   
#   points(y~x,type='l',col=palette()[r],lty=lty)
# 
#   }
  
  for(i in 1:300){#nrow(data$c)){
    y1<-exp(data$c[i,1]+x*data$c[i,2])
    if(!perCapita){y1<-y1*x}
    points(y1~x,type='l',col=errorColors[r],lwd=0.5)
  }
  points(y~x,type='l',col=palette()[r],lty=lty)
  }
}

allFiles<-list.files(fileDir)

bkt<-c(grep('meanStockRecruitBkt',allFiles),
  grep('extremeStockRecruitBkt',allFiles))

bnt<-c(grep('meanStockRecruitBnt',allFiles),
  grep('extremeStockRecruitBnt',allFiles))

bnt<-bnt[c(grep('Jimmy',allFiles[bnt]),grep('West',allFiles[bnt]))]


tiff.par("~/trout_yoy/results/figures/manuscript/stockRecruitExtremeError.tif",
         mfrow=c(2,2),mar=c(2.5,3.5,1,0),height=6.5,width=6.5,lend='butt')
#plot brook trout stock recruit
  plot(NA,xlim=c(0,maxAdult[4,1]),
       ylim=c(0,1000),
       xlab="Spawners (t-1)",ylab="")

  title(ylab="Recruits (t)",line=2)

  for(r in 5:8){
    plotStockRecruit(allFiles[bkt[r]],fun='points',error=T) 
  }
  for(r in 5:8){
    plotStockRecruit(allFiles[bkt[r]],fun='points',error=F) 
  }

#plot brook trout you per capita
  plot(NA,xlim=c(0,maxAdult[4,1]),
       ylim=c(0,4),
       xlab="Spawners (t-1)",ylab="Recruits Per Spawner (t)")
    
  for(r in 5:8){
    plotStockRecruit(allFiles[bkt[r]],fun='points',perCapita=T,error=T)
    
  }
  for(r in 5:8){
    plotStockRecruit(allFiles[bkt[r]],fun='points',perCapita=T,error=F)
  }
  legend(400,3,legend=riverNames,lty=1,lwd=2,
         col=palette()[1:4],bty='n')

#plot brown trout stock recruit  
    plot(NA,xlim=c(0,maxAdult[4,1]),
       ylim=c(0,1000),
       xlab="Spawners (t-1)",ylab="")

  title(ylab="Recruits (t)",line=2)

  for(r in c(2,4)){
    plotStockRecruit(allFiles[bnt[r]],fun='points',error=T) 
  }
  for(r in c(2,4)){
    plotStockRecruit(allFiles[bnt[r]],fun='points',error=F) 
  }

#plot brown trout yoy per capita
  plot(NA,xlim=c(0,maxAdult[4,1]),
       ylim=c(0,4),
       xlab="Spawners (t-1)",ylab="Recruits Per Spawner (t)")
    
  for(r in c(2,4)){
    plotStockRecruit(allFiles[bnt[r]],fun='points',perCapita=T,error=T)
    
  }
  for(r in c(2,4)){
    plotStockRecruit(allFiles[bnt[r]],fun='points',perCapita=T,error=F)
  }
#  legend(400,6.5,legend=riverNames,lty=1,lwd=2,
#         col=palette()[1:4],bty='n')
dev.off()