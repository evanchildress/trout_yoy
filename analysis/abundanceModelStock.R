data_root<-"~/process-data/data_store/processed_data"

load(file.path(data_root,"abundance_arrays.rDATA"))

detection<-readRDS(file.path(data_root,"detectionFromCJS.rds"))
covariates<-readRDS(file.path(data_root,"covariates.rds")) #load previously 



sumTempMean<-data.table(covariates[,,9])
sumTempFill<-lm(V4~V3+V2+V1,data=sumTempMean)
covariates[6:9,4,9]<-predict(sumTempFill,data.frame(sumTempMean[6:9,list(V3,V2,V1)]))
rm(sumTempMean)
rm(sumTempFill)
nsections<-c(14,14,15,47)

y<-y[,3:15,,,]
N<-y #use naive estimate as initial value
# N<-array(10,dim=dim(y))
# for(r in 1:3){
#   N[(nsections[r]+1):47,,r,,]<-NA
# }



cat("model{
  #Priors
  #State process
  for(g in 1:2){
    for(sp in 1:2){
      for(r in 1:4){
        mu[r,sp,g]~dnorm(0,0.01) #grand mean
  
  for(f in 1:6){
  beta[r,f,sp,g]~dnorm(0,0.01)
  }}}}


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
#   for(g in 1:2){

    for(sp in 1:2){  
      for(r in 1:4){
        for(j in 1:nsamples){
          eps[j,r,sp,2]~dnorm(0,tau.eps[r,2,sp])
        }
        tau.eps[r,sp,2]<-pow(sd.eps[r,2,sp],-2)
        sd.eps[r,sp,2]~dunif(0,5)
  }}

    for(sp in 1:2){  
      for(r in 1:4){
        for(j in 1:nsamples){
          eps[j,r,sp,1]~dnorm(0,tau.eps[r,1,sp])
        }
        tau.eps[r,sp,1]<-pow(sd.eps[r,1,sp],-2)
        sd.eps[r,sp,1]~dunif(0,5)
  }}


  #Initial year for yoy
  for(sp in 1:2){
    for(r in 1:4){
      for(i in 1:nsections[r]){
        logLambdaI[i,r,sp]~dnorm(0,0.01)
        lambdaI[i,r,sp]<-exp(logLambdaI[i,r,sp])
      }
    }
  }
  #Observation process
#   for(sp in 1:2){
#     for(r in 1:4){
#       mu.p[r,sp]~dnorm(0,0.01)
#     }
#   }
 
  #Likelihood for adults
  for(g in 2:2){
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
          for(i in 1:nsections[r]){
            N[i,j,r,sp,g]~dpois(lambda[i,j,r,sp,g])
              log(lambda[i,j,r,sp,g])<-mu[r,sp,g]+ alpha[i,r,sp,g]+eps[j,r,sp,g]

                                       +beta[r,1,sp,g]*covariates[j,r,1]

                                       #extreme covariates (either this or the means one should be commented out)
                                       +beta[r,2,sp,g]*covariates[j,r,2]+beta[r,3,sp,g]*covariates[j,r,3]#winter and spring high flow extremes
                                       +beta[r,4,sp,g]*covariates[j,r,4]+beta[r,5,sp,g]*covariates[j,1,5]#summer low flow and temp
                                       +beta[r,6,sp,g]*covariates[j,r,4]*covariates[j,r,5]#summer low flow/temp interaction
              
                                       #mean covariates
#                                        +beta[r,2,sp,g]*covariates[j,r,6]+beta[r,3,sp,g]*covariates[j,r,7]#winter and spring mean flow
#                                        +beta[r,4,sp,g]*covariates[j,r,8]+beta[r,5,sp,g]*covariates[j,1,9]#summer mean flow and temp
#                                        +beta[r,6,sp,g]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction


              #N2[i,j,r,sp,g]<-N[i,j,r,sp,g]-y[i,j,r,sp,1,g]
              y[i,j,r,sp,g]~dbin(p[j,r],N[i,j,r,sp,g])
              #y[i,j,r,sp,2,g]~dbin(p[j,r],N2[i,j,r,sp,g])
              #p[i,j,r,sp,g]<-exp(lp[i,j,r,sp,g])/(1+exp(lp[i,j,r,sp,g]))
              #lp[i,j,r,sp,g]<-mu.p[r,sp]#+alpha.p[i,r]
               }}}}}

  #Likelihood for first yoy year
    for(sp in 1:2){
      for(r in 1:4){
        for(i in 1:nsections[r]){
             N[i,1,r,sp,1]~dpois(lambdaI[i,r,sp])
             y[i,1,r,sp,1]~dbin(p[1,r],N[i,1,r,sp,1])
        }
      }
    }

  #Likelihood for yoy
  for(g in 1:1){
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 2:nsamples){
          for(i in 1:nsections[r]){
            N[i,j,r,sp,g]~dpois(lambda[i,j,r,sp,g])
              log(lambda[i,j,r,sp,g])<-mu[r,sp,g]+ alpha[i,r,sp,g]+eps[j,r,sp,g]
                                       +c[r,sp,g]*N[i,j-1,r,sp,2]
                                       +beta[r,1,sp,g]*covariates[j,r,1]

                                       #extreme covariates (either this or the means one should be commented out)
                                       +beta[r,2,sp,g]*covariates[j,r,2]+beta[r,3,sp,g]*covariates[j,r,3]#winter and spring high flow extremes
                                       +beta[r,4,sp,g]*covariates[j,r,4]+beta[r,5,sp,g]*covariates[j,1,5]#summer low flow and temp
                                       +beta[r,6,sp,g]*covariates[j,r,4]*covariates[j,r,5]#summer low flow/temp interaction
              
                                       #mean covariates
#                                        +beta[r,2,sp,g]*covariates[j,r,6]+beta[r,3,sp,g]*covariates[j,r,7]#winter and spring mean flow
#                                        +beta[r,4,sp,g]*covariates[j,r,8]+beta[r,5,sp,g]*covariates[j,1,9]#summer mean flow and temp
#                                        +beta[r,6,sp,g]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction


              #N2[i,j,r,sp,g]<-N[i,j,r,sp,g]-y[i,j,r,sp,1,g]
              y[i,j,r,sp,g]~dbin(p[j,r],N[i,j,r,sp,g])
              #y[i,j,r,sp,2,g]~dbin(p[j,r],N2[i,j,r,sp,g])
              #p[i,j,r,sp,g]<-exp(lp[i,j,r,sp,g])/(1+exp(lp[i,j,r,sp,g]))
              #lp[i,j,r,sp,g]<-mu.p[r,sp]#+alpha.p[i,r]
               }}}}}
  
  #Derived output
  for(g in 1:2){
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
          tot.pop[j,r,sp,g]<-sum(N[1:nsections[r],j,r,sp,g])
           

            for(i in 1:nsections[r]){
                 yExp[i,j,r,sp,g]<- N[i,j,r,sp,g]*p[j,r]
                 E[i,j,r,sp,g]<-pow((y[i,j,r,sp,g]-yExp[i,j,r,sp,g]),2)/
                                  (yExp[i,j,r,sp,g]+0.5)

                 yNew[i,j,r,sp,g]~dbin(p[j,r],N[i,j,r,sp,g])
                 ENew[i,j,r,sp,g]<-pow((yNew[i,j,r,sp,g]-yExp[i,j,r,sp,g]),2)/
                                        (yExp[i,j,r,sp,g]+0.5)
            }
  
      }}
    fit[sp,g]<-sum(E[1:nsections[1],,1,sp,g])+sum(E[1:nsections[2],,2,sp,g])+
                 sum(E[1:nsections[3],,3,sp,g])+sum(E[1:nsections[4],,4,sp,g])
    fitNew[sp,g]<-sum(ENew[1:nsections[1],,1,sp,g])+sum(ENew[1:nsections[2],,2,sp,g])+
                 sum(ENew[1:nsections[3],,3,sp,g])+sum(ENew[1:nsections[4],,4,sp,g])
  }}



}",file="nmixture.txt")

win.data<-list(y=y,
               nsamples=dim(y)[2],
               nsections=nsections,
               covariates=covariates,
               p=detection[3,,1:4])

inits<-function(){list(mu=array(rnorm(16,0,0.01),dim=c(4,2,2)),
                       sd.alpha=array(runif(16,0,4),dim=c(4,2,2)),
                       sd.eps=array(runif(16,0,4),dim=c(4,2,2)),
                       #mu.p=array(rnorm(8,0,0.01),dim=c(4,2)),sd.p=runif(4,0,3),N=N,
                       beta=array(rnorm(96,0,0.01),dim=c(4,6,2,2)),
                       N=N)}

params<-c("tot.pop","mu","sd.eps","beta","fit","fitNew")

ni=10000
nt=5
nb=7000
nc=3

#out<-bugs(win.data,inits,params,"nmixture.txt",n.chains=nc,n.iter=ni,
#          n.thin=nt,n.burnin=nb,debug=T,working.directory=getwd())
out<-jags(win.data,inits,params,"nmixture.txt",n.chains=nc,n.iter=ni,
          n.thin=nt,n.burnin=nb,working.directory=getwd())
saveRDS(out,file="~/trout_yoy/results/model_output.rds")

out.mcmc<-as.mcmc(out)




#windows()
#par(mfrow=c(5,6))
#plot(out.mcmc,density=F,auto.layout=T,ask=T)

popEst<-data.table(melt(out$BUGSoutput$sims.list$tot.pop))
setnames(popEst,c("sim","year","river","species","age","estimate"))
setkey(popEst,species,age,year)
popEst<-popEst[,list(round(mean(estimate),3),round(quantile(estimate,probs=0.025),3),
                   round(quantile(estimate,probs=0.975),3)),
             by=list(species,age,river,year)]
setnames(popEst,c("V1","V2","V3"),c("mean","lower","upper"))
popEst[,year:=year+2001]

betas<-data.table(melt(out$BUGSoutput$sims.list$beta))
setnames(betas,c("sim","river","beta","species","age","estimate"))
setkey(betas,beta,species,age)
betas<-betas[,list(round(mean(estimate),3),round(quantile(estimate,probs=0.025),3),
                   round(quantile(estimate,probs=0.975),3)),
                   by=list(beta,species,age,river)]
setnames(betas,c("V1","V2","V3"),c("mean","lower","upper"))


#get.names<-function(species,age,river){
species=c('bkt','bnt')
age=c('yoy','adult')
river=c('wb jimmy','wb mitchell','wb obear','west brook')
#return(list())


#set(betas,i=NULL,j=2:4,as.character(betas[,2:4,with=F]))
betas[,species:=as.character(species)]
betas[,age:=as.character(age)]
betas[,river:=as.character(river)]



for(i in 1:nrow(betas)){
  set(betas,i=as.integer(i),j=2:4,as.list(c(species[as.numeric(betas[i,species])],
                        age[as.numeric(betas[i,age])],
                        river[as.numeric(betas[i,river])])))
  
  betas[,sig:=((lower<0 & upper <0) | (lower>0 & upper>0))]
}

tiff.par("~/trout_yoy/results/figures/parameter_estimates.tif",width=6.5,height=6.5,mar=c(1.7,3.1,0.5,0),mgp=c(2.2,0.5,0),
         mfrow=c(4,1))
for(sp in 1:2){
  for(g in 1:2){
    plot(mean~I(1:length(mean)),col=as.factor(river),pch=19,cex=1.5,
         data=betas[species==unique(species)[sp]&age==unique(age)[g]],
         xaxt='n',xlab="",ylab='parameter estimate',bty='l',ylim=c(-2,2),
         main=paste(unique(species)[sp],unique(age)[g]))
    
    if(g==1 & sp==1) legend(20,2,c("J","M","O'B","WB"),pt.cex=1.5,pch=19,col=palette()[1:4],bty='n')
    
    with(betas[species==unique(species)[sp]&age==unique(age)[g]],
         error.bar(1:length(mean),mean,upper.y=upper,lower.y=lower,interval.type='not distance'))
    axis(1,at=seq(2.5,22.5,4),labels=c("fall.q","winter.q","spring.q","sum.temp","sum.q","sum.temp x sum.q"))
    abline(h=0,lty=2)
    abline(v=seq(4.5,24.5,4))
  }
}
dev.off()


#windows()
for(g in 1:2){
  for(sp in 1:2){
tiff.par(paste0("~/trout_yoy/results/figures/population_estimates_",
                age[g],species[sp],".tif"),
         width=6.5,height=6.5,
         mar=c(1.7,3.1,1,0),mgp=c(2.2,0.5,0),
         mfrow=c(4,1))

    for(r in 1:4){
      ylimit<-c(0,popEst[river   == unique(river)[r]&
                         species == unique(species)[sp]&
                         age     == unique(age)[g],
                         max(upper)+10])
      plot(mean~year,main=get("river",env=.GlobalEnv)[r],
           data=popEst[river   == unique(river)[r]&
                       species == unique(species)[sp]&
                       age     == unique(age)[g]
                                 ],
           type='b',
           pch=19,
           lty=1,
           col=palette()[r],
           ylim=ylimit)
      
      with(popEst[river   == unique(river)[r]&
                  species == unique(species)[sp]&
                  age     == unique(age)[g]
                  ],
           error.bar(year,mean,upper.y=upper,lower.y=lower,interval.type="addota"))
    }
      dev.off()
  }
}


plot(out$BUGSoutput$sims.list$fitNew~out$BUGSoutput$sims.list$fit,
     xlab="Discrepancy Actual Data",ylab="Discrepancy Replicated Data",
     bty='l')
abline(0,1,lwd=2)
mean(out$BUGSoutput$sims.list$fitNew>out$BUGSoutput$sims.list$fit)
mean(out$BUGSoutput$sims.list$fitNew<out$BUGSoutput$sims.list$fit)

#plot(envPred[,1,4]~c(2002:2014),type='b',col='blue',pch=19,lwd=2,ylim=c(-2,2.5))
#points(envPred[,2,4]~c(2002:2014),type='b',col='pink',pch=19,lwd=2)
#points(envPred[,3,4]~c(2002:2014),type='b',col='black',pch=19,lwd=2)
#points(envPred[,4,4]~c(2002:2014),type='b',col='orange',pch=19,lwd=2)
#points(envPred[,5,4]~c(2002:2014),type='b',col='red',pch=19,lwd=2)

