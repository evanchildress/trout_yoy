source('createAbundanceModel.R')
data_root<-"~/process-data/data_store/processed_data"

load(file.path(data_root,"abundance_arrays.rDATA"))

# detection<-readRDS(file.path(data_root,"detectionFromCJS.rds"))
# detection<-detection[3,,1:4]

bktDetect<-readRDS(file.path("~/westbrookJS/results/pYoy.rds"))
bktDetect<-acast(melt(bktDetect[season==3,list(year,river,mean)],
                      id.vars=c("river","year")),
                 year~river)
bktDetect<-bktDetect[2:14,]

bntDetect<-readRDS(file.path("~/westbrookJS/resultsBnt/pYoy.rds"))
bntDetect<-acast(melt(bntDetect[season==3,list(year,river,mean)],
                      id.vars=c("river","year")),
                 year~river)
bntDetect<-bntDetect[2:14,]

detection<-array(c(bktDetect,bntDetect),
                 c(nrow(bktDetect),4,2))

propTagged<-propTagged[3:15,,]

detection<-detection/propTagged #correct for capture of untagged fish (that could be tagged later)
detection[13,,]<-apply(detection[1:12,,],c(2,3),mean)
detection[detection>1]<-1#bnt 2006 is 1.01 after adjustment, so this sets it to 1

covariates<-readRDS(file.path(data_root,"covariates.rds")) #load env covariates created in README

sumTempMean<-data.table(covariates[,,9])
sumTempFill<-lm(V4~V3+V2+V1,data=sumTempMean)
covariates[6:9,4,9]<-predict(sumTempFill,data.frame(sumTempMean[6:9,list(V3,V2,V1)]))
rm(sumTempMean)
rm(sumTempFill)
nsections<-c(14,14,15,47)

y<-y[,3:15,,,1]
N<-y #use naive estimate as initial value
# N<-array(10,dim=dim(y))
# for(r in 1:3){
#   N[(nsections[r]+1):47,,r,,]<-NA
# }





win.data<-list(y=y,
               nsamples=dim(y)[2],
               nsections=nsections,
               covariates=covariates,
               p=detection)

inits<-function(){list(mu=array(rnorm(8,0,0.01),dim=c(4,2)),
                       sd.alpha=array(runif(8,0,4),dim=c(4,2)),
                       sd.eps=array(runif(8,0,4),dim=c(4,2)),
                       #mu.p=array(rnorm(8,0,0.01),dim=c(4,2)),sd.p=runif(4,0,3),N=N,
                       beta=array(rnorm(48,0,0.01),dim=c(4,6,2)),
                       N=N)}

params<-c("tot.pop","mu","sd.eps","beta","fit","fitNew",'tau.alpha','sd.alpha')

ni=10000
nt=5
nb=7000
nc=3

#out<-bugs(win.data,inits,params,"nmixture.txt",n.chains=nc,n.iter=ni,
#          n.thin=nt,n.burnin=nb,debug=T,working.directory=getwd())
out<-jags(win.data,inits,params,"abundanceModel.txt",n.chains=nc,n.iter=ni,
          n.thin=nt,n.burnin=nb,working.directory=getwd())
saveRDS(out,file="~/trout_yoy/results/model_output.rds")

sims<-out$BUGSoutput$sims.list




#windows()
#par(mfrow=c(5,6))
#plot(out.mcmc,density=F,auto.layout=T,ask=T)

popEst<-data.table(melt(sims$tot.pop))
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


#get.names<-function(species,age,river){
species=c('bkt','bnt')
#age=c('yoy','adult')
river=c('wb jimmy','wb mitchell','wb obear','west brook')
#return(list())


#set(betas,i=NULL,j=2:4,as.character(betas[,2:4,with=F]))
#betas[,species:=as.character(species)]
#betas[,age:=as.character(age)]
#betas[,river:=as.character(river)]



# for(i in 1:nrow(betas)){
#   set(betas,i=as.integer(i),j=2:3,as.list(c(species[as.numeric(betas[i,species])],
#                         #age[as.numeric(betas[i,age])],
#                         river[as.numeric(betas[i,river])])))
# }

betas[,sig:=((lower<0 & upper <0) | (lower>0 & upper>0))]

tiff.par("~/trout_yoy/results/figures/parameter_estimates.tif",width=6.5,height=6.5,mar=c(1.7,3.1,0.5,0),mgp=c(2.2,0.5,0),
         mfrow=c(2,1))
for(sp in 1:2){

    plot(mean~I(1:length(mean)),col=as.factor(river),pch=19,cex=1.5,
         data=betas[species==unique(species)[sp]],
         xaxt='n',xlab="",ylab='parameter estimate',bty='l',ylim=c(-2,2),
         main=paste(unique(species)[sp]))
    
    if(sp==1) legend(20,2,c("J","M","O'B","WB"),pt.cex=1.5,pch=19,col=palette()[1:4],bty='n')
    
    with(betas[species==unique(species)[sp]],
         error.bar(1:length(mean),mean,upper.y=upper,lower.y=lower,interval.type='not distance'))
    axis(1,at=seq(2.5,22.5,4),labels=c("fall.q","winter.q","spring.q","sum.temp","sum.q","sum.temp x sum.q"))
    abline(h=0,lty=2)
    abline(v=seq(4.5,24.5,4))
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


tiff.par("results/figures/posteriorPredictiveCheck.tif")
plot(sims$fitNew[,,1]~sims$fit[,,1],
     xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
     bty='l')
for(r in 1:4){
  points(sims$fitNew[,r,1]~sims$fit[,r,1],col=palette()[r])
}
abline(0,1,lwd=2)
text(200,300,paste0("YOY BKT Bayesian P = ",
                    round(mean(sims$fitNew[,,1]>sims$fit[,,1]),3)))
dev.off()

#plot(envPred[,1,4]~c(2002:2014),type='b',col='blue',pch=19,lwd=2,ylim=c(-2,2.5))
#points(envPred[,2,4]~c(2002:2014),type='b',col='pink',pch=19,lwd=2)
#points(envPred[,3,4]~c(2002:2014),type='b',col='black',pch=19,lwd=2)
#points(envPred[,4,4]~c(2002:2014),type='b',col='orange',pch=19,lwd=2)
#points(envPred[,5,4]~c(2002:2014),type='b',col='red',pch=19,lwd=2)

