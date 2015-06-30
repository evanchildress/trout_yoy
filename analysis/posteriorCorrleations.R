postCor<-function(x,vars=c('beta','c'),modelName){
  require(reshape2)
  require(data.table)
  for(v in vars){
    assign(v,
      data.table(melt(x[[v]])))
    get(v)[,Var2:=paste0(v,Var2)]
  }

  #rivers<-c("jimmy","mitchell","obear","west brook")
  
  data<-do.call(rbind,args=as.list(c(quote(beta),quote(c)))) #can this be done using vars?
  
  data<-data.table(dcast(data,Var1~Var2))
  setnames(data,c("Var1"),c("iter"))
  
  corMatrix<-cor(data[,2:ncol(data),with=F])
  
#   par(ask=F)
#   layout(matrix(c(1,2),1,2),widths=c(3,1),heights=c(4,4))
#    for(i in 3:(ncol(data)-1)){
#      for(d in (i+1):ncol(data)){
#        par(mar=c(3,3,0,0),mgp=c(2,0.5,0),las=1)
#        plot(data[[names(data)[i]]],
#             data[[names(data)[d]]],
#             col=as.factor(data$river),
#             bty='l',xlab=names(data)[i],ylab=names(data)[d])
#        par(mar=c(0,0,0,0))
#        plot(NA,xlim=c(0,1),ylim=c(0,5),axes=F,xlab=NA,ylab=NA)
#        for(r in 1:4){
#          a<-lm(data[river==r][[names(data)[d]]]~
#                  data[river==r][[names(data)[i]]])
#          text(0.5,r,
#               bquote(.(rivers[r]) ~ "=" ~ 
#                     .(round(summary(a)$r.squared,2))),
#               col=palette()[r])
#          corByRiver[i-2,d-2,r]<-summary(a)$r.squared
#        }
#      }
#    }
#   
#   meanCor<-array(NA,dim=c(ncol(corByRiver),4))
#   for(r in 1:4){
#     for(i in 1:ncol(corByRiver)){
#       meanCor[i,r]<-mean(c(corByRiver[i,,r],corByRiver[,i,r]),na.rm=T)
#     }
#   }
#   #return(corMatrix)
#   par(ask=F)
#   tiff.par(file.path("~/trout_yoy/results/figures",
#                      paste0("posteriorCor",modelName,".tif")),
#            mfrow=c(6,10),mar=c(0.2,0.2,1,0))
#   for(i in 1:ncol(corByRiver)){
#     for(d in 1:nrow(corByRiver)){
#       if(is.na(corByRiver[i,d,1])){next}
#       barplot(corByRiver[i,d,],col=palette()[1:4],names.arg=NA,axes=F,ylim=c(0,1))
#       title(main=do.call(paste,args=as.list(dimnames(corByRiver)[[1]][c(i,d)])),
#             cex.main=0.8)
#     }
#   }
#   for(i in 1:nrow(meanCor)){
#     barplot(meanCor[i,],col=palette()[1:4],axes=F,xlab=NA,ylab=NA,ylim=c(0,1))
#     title(main=dimnames(corByRiver)[[1]][i])
#   }
#   dev.off()
  return(corMatrix)
}

riverNames<-c("Jimmy","Mitchell","Obear","WestBrook")
modelType<-c("meanStockRecruit","extremeStockRecruit")
for(r in riverNames){
  for(m in modelType){
    modelName<-paste0(m,"Bkt",r,".rds")
    sims<-readRDS(file.path("~/trout_yoy/results/modelOutput/speciesRiverSeparate",
                            modelName))$BUGSoutput$sims.list
    assign(paste0(m,r),postCor(sims,modelName=modelName))
  }
}

# 
# beta1=fall discharge
# beta2=winter discharge
# beta3=spring discharge
# beta4=summer flow
# beta5=summer temp
# beta6 summer flow*temp
# beta7=fall discharge^2
# beta8=winter discharge^2+
