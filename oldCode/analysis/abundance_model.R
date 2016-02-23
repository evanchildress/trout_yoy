source('~/trout_yoy/createAbundanceModel.R')
load('~/trout_yoy/abundanceData.RData')

sumTempMean<-data.table(covariates[,,9])
sumTempFill<-lm(V4~V3+V2+V1,data=sumTempMean)
covariates[6:9,4,9]<-predict(sumTempFill,data.frame(sumTempMean[6:9,list(V3,V2,V1)]))
rm(sumTempMean)
rm(sumTempFill)
#nsections<-c(14,14,15,47)
#covariates[,,10]<-covariates


y<-y[3:15,,,1]

covariates<-covariates[1:13,,]
detection<-detection[1:13,,]
A<-A[1:14,,]

N<-y #use naive estimate as initial value

win.data<-list(yDATA=y,
               nsamples=dim(y)[1],
               #nsections=nsections,
               covariates=covariates,
               p=detection,
               adultDATA=A)

inits<-function(){list(c=array(c(runif(8,0,5),rep(0,16)),dim=c(4,2,3)),
                       beta=array(rnorm(64,0,0.01),dim=c(4,8,2)),
                       N=N
)}

params<-c("N","c","beta","fit","fitNew","yExp")

ni=15000
nt=5
nb=10000
nc=3

#out<-bugs(win.data,inits,params,"nmixture.txt",n.chains=nc,n.iter=ni,
#          n.thin=nt,n.burnin=nb,debug=T,working.directory=getwd())
out<-jags(win.data,inits,params,"abundanceModel.txt",n.chains=nc,n.iter=ni,
          n.thin=nt,n.burnin=nb,working.directory=getwd())
saveRDS(out,file="~/trout_yoy/results/modelOutputExtreme.rds")

sims<-out$BUGSoutput$sims.list



stockRecruit<-data.table(melt(sims$c))
setnames(stockRecruit,c("sim","river","species","parameter","estimate"))
setkey(stockRecruit,parameter,species,river)
stockRecruit<-stockRecruit[,list(round(mean(estimate),3),
                                 round(quantile(estimate,probs=0.025),3),
                                 round(quantile(estimate,probs=0.975),3)),
             by=list(parameter,species,river)]
setnames(stockRecruit,c("V1","V2","V3"),c("mean","lower","upper"))

popEst<-data.table(melt(sims$N))
setnames(popEst,c("sim","year","river","species","estimate"))
setkey(popEst,species,year)
popEst<-popEst[,list(round(mean(estimate),3),round(quantile(estimate,probs=0.025),3),
                   round(quantile(estimate,probs=0.975),3)),
             by=list(species,river,year)]
setnames(popEst,c("V1","V2","V3"),c("mean","lower","upper"))
popEst[,year:=year+2001]

betas<-data.table(melt(sims$beta))
setnames(betas,c("sim","river","beta","species","estimate"))
setkey(betas,beta,species)
betas<-betas[,list(round(mean(estimate),3),round(quantile(estimate,probs=0.025),3),
                   round(quantile(estimate,probs=0.975),3)),
                   by=list(beta,species,river)]
setnames(betas,c("V1","V2","V3"),c("mean","lower","upper"))



species=c('bkt','bnt')

river=c('wb jimmy','wb mitchell','wb obear','west brook')


betas[,sig:=((lower<0 & upper <0) | (lower>0 & upper>0))]

tiff.par("~/trout_yoy/results/figures/envParameterEstimates.tif",
         width=6.5,height=6.5,mar=c(1.7,3.0,1,0),mgp=c(2.2,0.5,0),
         mfrow=c(2,1))
for(sp in 1:2){

    plot(mean~I(1:length(mean)),col=as.factor(river),pch=19,cex=1.5,
         data=betas[species==unique(species)[sp]],
         xaxt='n',xlab="",ylab='parameter estimate',bty='l',ylim=c(-2,2),
         main=species[sp])
    
    if(sp==1) legend(20,2,c("J","M","O'B","WB"),pt.cex=1.5,pch=19,col=palette()[1:4],bty='n')
    
    with(betas[species==unique(species)[sp]],
         error.bar(1:length(mean),mean,upper.y=upper,lower.y=lower,interval.type='not distance'))
    axis(1,at=seq(2.5,22.5,4),labels=c("fall.q","winter.q","spring.q","sum.temp","sum.q","sum.temp x sum.q"))
    abline(h=0,lty=2)
    abline(v=seq(4.5,24.5,4),col='gray')
  }

dev.off()


tiff.par("~/trout_yoy/results/figures/stockRecruitEstimates.tif",
         width=6.5,height=6.5,mar=c(1.7,3.0,1,0),mgp=c(2.2,0.5,0),
         mfrow=c(2,4))

lims<-list(c(0,5),c(-2,2),c(-0.2,0.2))
for(sp in 1:2){
  for(p in 1:3){

    plot(mean~I(1:length(mean)),col=as.factor(river),pch=19,cex=1.5,
         data=stockRecruit[species==sp & parameter==p],
         xaxt='n',xlab="",ylab='parameter estimate',bty='l',
         ylim=lims[[p]]*ifelse(sp==1 & p==2,0.1,1),
         main=paste0(c('bkt','bnt')[sp]," c",p))
    
    if(sp==1 & p==1) legend(20,2,c("J","M","O'B","WB"),pt.cex=1.5,pch=19,col=palette()[1:4],bty='n')
    
    with(stockRecruit[species==sp & parameter==p],
         error.bar(1:length(mean),mean,upper.y=upper,lower.y=lower,interval.type='not distance'))
    abline(h=0,lty=2)
  }
  
  plot(NA,xlim=c(0,ceiling(max(A[,,sp]))),ylim=c(0,10),bty='l',
      xlab='Spawners',ylab='Per Capita Recruits')
  stockPredict<-data.table(x=seq(0,ceiling(max(A[,,sp])),1))
  for(r in 1:4){
    stockPredict[,y:=as.numeric(NA)]
    stockPredict[x<max(A[,r,sp]),y:=exp(stockRecruit[parameter==1 & species==sp & river==r,mean]+
            stockRecruit[parameter==2 & species==sp & river==r,mean]*x)]
    points(y~x,type='l',col=palette()[r],data=stockPredict)
  }
}
dev.off()


#windows()

  for(sp in 1:2){
tiff.par(paste0("~/trout_yoy/results/figures/population_estimates_",
                'yoy',species[sp],".tif"),
         width=6.5,height=6.5,
         mar=c(1.7,3.1,1,0),mgp=c(2.2,0.5,0),
         mfrow=c(4,1))

    for(r in 1:4){
      ylimit<-c(0,popEst[river   == unique(river)[r]&
                         species == unique(species)[sp]
                         #&age     == unique(age)[g]
                         ,max(upper)+10])
      plot(mean~year,main=get("river",env=.GlobalEnv)[r],
           data=popEst[river   == unique(river)[r]&
                       species == unique(species)[sp]
                       #& age     == unique(age)[g]
                                 ],
           type='b',
           pch=19,
           lty=1,
           col=palette()[r],
           ylim=ylimit)
      
      with(popEst[river   == unique(river)[r]&
                  species == unique(species)[sp]
                  #& age     == unique(age)[g]
                  ],
           error.bar(year,mean,upper.y=upper,lower.y=lower,interval.type="addota"))
    }
      dev.off()
  }


tiff.par("~/trout_yoy/results/figures/posteriorPredictiveCheck.tif",
         mfrow=c(1,2),width=8,height=4,mar=c(2.5,2.5,1,0))
plot(sims$fitNew[,,1]~sims$fit[,,1],
     xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
     bty='l',main='bkt')
for(r in 1:4){
  points(sims$fitNew[,r,1]~sims$fit[,r,1],col=palette()[r])
}
abline(0,1,lwd=2)
text(max(sims$fit[,,1])*0.9,max(sims$fitNew[,,1])*0.9,
     round(mean(sims$fitNew[,,1]>sims$fit[,,1]),3))

plot(sims$fitNew[,,2]~sims$fit[,,2],
     xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
     bty='l',main='bnt',pch=NA)
for(r in 1:4){
  points(sims$fitNew[,r,2]~sims$fit[,r,2],col=palette()[r])
}
abline(0,1,lwd=2)
text(max(sims$fit[,,2])*0.9,max(sims$fitNew[,,2])*0.9,
            round(mean(sims$fitNew[,,2]>sims$fit[,,2]),3))
dev.off()

saveRDS(popEst,"~/trout_yoy/results/popEstEnv.rds")
