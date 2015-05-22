cat("model{
  #Priors
  #State process
    for(sp in 1:2){
      for(r in 1:4){

    #stock-recruit
        c[r,sp,1]~dunif(0,50)
        c[r,sp,2]~dnorm(0,0.001)
        c[r,sp,3]~dnorm(0,0.001)
      }
    }
  
  #Likelihood for yoy
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
            N[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
              log(lambda[j,r,sp])<-
                          #stock-recruit
                          c[r,sp,1]+c[r,sp,2]*adultDATA[j,r,sp]+c[r,sp,3]*pow(adultDATA[j,r,sp],2)

              yDATA[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])

        }
      }
    }
  
  #Derived output
    for(sp in 1:2){
      for(r in 1:4){
        for(j in 1:nsamples){
                 yExp[j,r,sp]<- N[j,r,sp]*p[j,r,sp]
                 E[j,r,sp]<-pow((yDATA[j,r,sp]-yExp[j,r,sp]),2)/
                                  (yExp[j,r,sp]+0.5)
                
                 #NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                 yNew[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])
                 ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                        (yExp[j,r,sp]+0.5)
        }
        fit[r,sp]<-sum(E[,r,sp])
        fitNew[r,sp]<-sum(ENew[,r,sp])
      }
    }
}",file="abundanceModel.txt")