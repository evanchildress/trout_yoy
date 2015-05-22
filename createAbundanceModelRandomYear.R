cat("model{
  #Priors
  
  #Random year effect

    for(sp in 1:2){  
      for(r in 1:4){
        for(j in 1:nsamples){
          eps[j,r,sp]~dnorm(0,tau.eps[r,sp])
        }
        tau.eps[r,sp]<-pow(sd.eps[r,sp],-2)
        sd.eps[r,sp]~dunif(0,10)
      }
    }

  #Likelihood for yoy
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
            N[j,r,sp]~dpois(lambda[j,r,sp])
              log(lambda[j,r,sp])<-eps[j,r,sp]

              yDATA[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])

        }
      }
    }
  
  #Derived output
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
                 yExp[j,r,sp]<-N[j,r,sp]*p[j,r,sp]
                 E[j,r,sp]<-pow((yDATA[j,r,sp]-yExp[j,r,sp]),2)/
                                  (yExp[j,r,sp]+0.5)
                
                 #NNew[j,r,sp]~dpois(lambda[j,r,sp])
                 yNew[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])
                 ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                  (yExp[j,r,sp]+0.5)
        }
        fit[r,sp]<-sum(E[,r,sp])
        fitNew[r,sp]<-sum(ENew[,r,sp])
      }
    }
}",file="abundanceModel.txt")