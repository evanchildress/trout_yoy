data_root<-"~/process-data/data_store/processed_data"

load(file.path(data_root,"abundance_arrays.rDATA"))

covariates<-readRDS(file.path(data_root,"covariates.rds"))

nsections<-c(14,14,15,47)

N<-array(0,dim=dim(y)[c(1:4,6)])

N<-N+apply(y,c(1,2,3,4,6),sum,na.rm=T)
for(r in 1:3){
  N[(nsections[r]+1):47,,r,,]<-NA
}

cat("model{
  #Priors
  #State process
  for(g in 1:2){
    for(sp in 1:2){
      for(r in 1:4){
        mu[r,sp,g]~dnorm(0,0.01) #grand mean
  
  #for(f in 1:4){
  #beta[r,f,sp,g]~dnorm(0,0.01)
  }}}


  #Random section effect
  for(g in 1:2){
    for(sp in 1:2){  
      for(r in 1:4){
        for(i in 1:nsections[r]){
          alpha[i,r,sp,g]~dnorm(0,tau.alpha[r,g,sp])
          }
      tau.alpha[r,sp,g]<-pow(sd.alpha[r,g,sp],-2)
      sd.alpha[r,sp,g]~dunif(0,5)
  }}}
  
  #Random year effect
  for(g in 1:2){
    for(sp in 1:2){  
      for(r in 1:4){
        for(j in 1:nsamples){
          eps[j,r,sp,g]~dnorm(0,tau.eps[r,g,sp])
        }
        tau.eps[r,sp,g]<-pow(sd.eps[r,g,sp],-2)
        sd.eps[r,sp,g]~dunif(0,5)
  }}}

  #Observation process
  
  for(r in 1:4){
    mu.p[r]~dnorm(0,0.01)
    for(i in 1:nsections[r]){
      alpha.p[i,r]~dnorm(0,tau.p[r])
    }
    tau.p[r]<-pow(sd.p[r],-2)
    sd.p[r]~dunif(0,5)
  }
  
  #Likelihood for yoy
  for(g in 1:2){
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
          for(i in 1:nsections[r]){
            N[i,j,r,sp,g]~dpois(lambda[i,j,r,sp,g])
              log(lambda[i,j,r,sp,g])<-mu[r,sp,g]

              N2[i,j,r,sp,g]<-N[i,j,r,sp,g]-y[i,j,r,sp,1,g]
              y[i,j,r,sp,1,g]~dbin(p[i,j,r,sp,g],N[i,j,r,sp,g])
              y[i,j,r,sp,2,g]~dbin(p[i,j,r,sp,g],N2[i,j,r,sp,g])
              p[i,j,r,sp,g]<-exp(lp[i,j,r,sp,g])/(1+exp(lp[i,j,r,sp,g]))
              lp[i,j,r,sp,g]<-mu.p[r]+alpha.p[i,r]
               }}}}}
  
  #Derived output
  for(g in 1:2){
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
          tot.pop[j,r,sp,g]<-sum(N[1:nsections[r],j,r,sp,g])
  }}}}
}",file="nmixture.txt")

win.data<-list(y=y,nsamples=dim(y)[2],nsections=nsections)

inits<-function(){list(mu=array(rnorm(16,0,0.01),dim=c(4,2,2)),
                       sd.alpha=array(runif(16,0,4),dim=c(4,2,2)),
                       sd.eps=array(runif(16,0,4),dim=c(4,2,2)),
                       mu.p=rnorm(4,0,0.01),sd.p=runif(4,0,3),N=N)}

params<-c("tot.pop","mu","sd.alpha","sd.eps","mu.p","sd.p")

ni=5000
nt=2
nb=3000
nc=3

#out<-bugs(win.data,inits,params,"nmixture.txt",n.chains=nc,n.iter=ni,
#          n.thin=nt,n.burnin=nb,debug=T,working.directory=getwd())
out<-jags(win.data,inits,params,"nmixture.txt",n.chains=nc,n.iter=ni,
          n.thin=nt,n.burnin=nb,working.directory=getwd())


out.mcmc<-as.mcmc(out)

#windows()
#par(mfrow=c(5,6))
#plot(out.mcmc,density=F,auto.layout=T,ask=T)

popEst<-apply(out$BUGSoutput$sims.list$tot.pop,c(2,3,4,5),mean)
dimnames(popEst)<-list(c(2000:2014),c("jimmy","mitchell","obear","west"),c("bkt","bnt"),c("yoy","adult"))
popEst[,year:=]

#windows()

tiff.par("Figures/Population Estimates.tif",width=6.5,height=6.5,mar=c(1.7,3.1,0,0),mgp=c(2.2,0.5,0))
plot(west~year,data=popEst,ylim=c(0,550),type='b',pch=19,ylab="YOY Population Estimate",xlab="")
points(jimmy~year,data=popEst,type='b',pch=19,col='blue')
points(mitchell~year,data=popEst,type='b',pch=19,col='red')
points(obear~year,data=popEst,type='b',pch=19,col='gray')
legend(2005.5,550,c("WB","J","M","O'B"),lty=1,pch=19,col=c('black',"blue","red","gray"),bty='n')
dev.off()

#plot(envPred[,1,4]~c(2002:2014),type='b',col='blue',pch=19,lwd=2,ylim=c(-2,2.5))
#points(envPred[,2,4]~c(2002:2014),type='b',col='pink',pch=19,lwd=2)
#points(envPred[,3,4]~c(2002:2014),type='b',col='black',pch=19,lwd=2)
#points(envPred[,4,4]~c(2002:2014),type='b',col='orange',pch=19,lwd=2)
#points(envPred[,5,4]~c(2002:2014),type='b',col='red',pch=19,lwd=2)


windows()
par(mfrow=c(4,4))
for(b in 1:4){
for(r in 1:4){
  plot(density(out$BUGSoutput$sims.list$beta[,r,b]),main=names(popEst)[r])
}
}

for(i in 1:4){
  print(cor(out$BUGSoutput$sims.list$beta[,i,]))
}

