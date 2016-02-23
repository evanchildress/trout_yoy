runModel<-function(data,stockRecruit=F,env='none',randomYear=F,species){
  if(species=='bkt'){sp<-1
  } else {if(species=='bnt'){sp<-2}}
  speciesName<-paste0(toupper(substr(species,1,1)),substr(species,2,3))
  resultsDir<-"~/trout_yoy/results/"
  
  outName<-ifelse(randomYear,paste0('randomYear',speciesName),
                  ifelse(env=='none',paste0('stockRecruit',speciesName),
                         ifelse(stockRecruit,paste0(env,'StockRecruit',speciesName),
                                paste0(env,speciesName))))

  createModel(stockRecruit=stockRecruit,env=env,
              randomYear=randomYear)

  if(env=='pca'){
  covPCA<-array(NA,dim=dim(covariates))
    pca<-prcomp(covariates,center=T,scale=T)
    scores<-pca$x
#     loadings<-pca$rotation
#     loadings[,1]<-loadings[,1]
    covPCA<-scores
    
  covariates<-covPCA
#     tiff.par(file="~/trout_yoy/results/figures/envPca.tif",
#          height=4,width=8,mai=c(.8,.8,0,0),mar=c(3,3,0,0),
#          mfrow=c(1,2))
#     
#       plot(scores[,1], scores[,2], xlab="", ylab="",
#       type="n", xlim=c(-8,8),
#       ylim=c(-8,5),bty='l')
#       arrows(0,0,loadings[,1]*10,loadings[,2]*10, length=0.1,
#       angle=20, col="red")
#       title(xlab="PC1",ylab="PC2",line=2)
#       text(loadings[,1]*10*1.2,loadings[,2]*10*1.2,
#       rownames(loadings), col="red", cex=0.7)
#       points(scores[,1],scores[,2], col="blue",cex=0.7,type="p")
#     
#       plot(scores[,3], scores[,4], xlab="", ylab="",
#       type="n", xlim=c(-8,10),
#       ylim=c(-6,8),bty='l')
#       arrows(0,0,loadings[,3]*10,loadings[,4]*10, length=0.1,
#       angle=20, col="red")
#       title(xlab="PC3",ylab="PC4",line=2)
#       text(loadings[,3]*10*1.2,loadings[,4]*10*1.2,
#       rownames(loadings), col="red", cex=0.7)
#       points(scores[,3],scores[,4], col="blue",cex=0.7,type="p")
#     
#     dev.off()
    
  }
  
  N<-data$y #use naive estimate as initial value
  
#Add orthogonal versions of quadratics for relevant variables
#   covariates<-array(c(covariates,rep(NA,dim(covariates)[1]*4)),
#                     dim=c(13,15))
#   covariates[,12]<-poly(covariates[,1],2)[,2]
#   covariates[,13]<-poly(covariates[,2],2)[,2]
#   covariates[,14]<-poly(covariates[,6],2)[,2]
#   
#     if(env=='extreme'){
#   covariates[,15]<-scale(covariates[,4]*covariates[,5])[,1]
#   } else if(env=='mean'){
#     covariates[,15]<-scale(covariates[,8]*covariates[,9])[,1]
#   }
  
  # Acentered<-A#scale(A)[,1]
  
#   win.data<-list(yDATA=y,
#                  nsamples=length(y),
#                  covariates=covariates,
#                  p=detection,
#                  adultDATAcentered=Acentered,
#                  adultDATA=A)
#                  #,otherSpDATA=otherSp)
  
  inits<-function(){list(N=N)
  }
  
  params<-c("N","c","nPredict","rmse","yExp")
  if(env!='none'){
      params<-c(params,"beta")
    }
  
  ni=10000
  nt=5
  nb=7000
  nc=3
  
  out<-jags(data,inits=inits,params,"~/trout_yoy/abundanceModel.txt",n.chains=nc,n.iter=ni,
            n.thin=nt,n.burnin=nb)
  saveRDS(out,file.path(resultsDir,paste0("modelOutput/",outName,".rds")))
#   
#   sims<-out$BUGSoutput$sims.list
  
################### make posterior predictive check figure ###################
# tiff.par(file.path(resultsDir,paste0("figures/predictiveChecks/bayesP_",outName,".tif")),
#            mfrow=c(1,1),width=4,height=4,mar=c(2.5,2.5,1,0))
#   plot(sims$fitNew[,1]~sims$fit[,1],
#        xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
#        bty='l',main='bkt',ylim=c(min(sims$fitNew),max(sims$fitNew)),
#        xlim=c(min(sims$fit),max(sims$fit)))
#   for(r in 1:4){
#     points(sims$fitNew[,r]~sims$fit[,r],col=palette()[r])
#   }
#   abline(0,1,lwd=2)
#   text(max(sims$fit[,])*0.9,max(sims$fitNew[,])*0.9,
#        round(mean(sims$fitNew[,]>sims$fit[,]),3))
#   
# #   plot(sims$fitNew[,,2]~sims$fit[,,2],
# #        xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
# #        bty='l',main='bnt',pch=NA)
# #   for(r in 1:4){
# #     points(sims$fitNew[,r,2]~sims$fit[,r,2],col=palette()[r])
# #   }
# #   abline(0,1,lwd=2)
# #   text(max(sims$fit[,,2])*0.9,max(sims$fitNew[,,2])*0.9,
# #               round(mean(sims$fitNew[,,2]>sims$fit[,,2]),3))
#   dev.off()
#####################################################################  
}