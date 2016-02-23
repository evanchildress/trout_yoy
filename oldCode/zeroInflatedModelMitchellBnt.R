  cat("model{
    #Priors
    #State process
      
   #zero year fall Q covariate
    zbeta~dnorm(0,0.01)

          for(f in 1:8){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
  #stock-recruit(c[,,1] and c[,,2])/competition parameters(c[,,3])
          c[1]~dunif(0,50)
          c[2]~dnorm(0,0.01)
        
  
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              z[j]~dbern(omega[j])
              logit(omega[j])<-zbeta*covariates[j,1]

              N[j]~dpois(lambda[j])#*adultDATA[j]*z[j])
                log(lambda[j])<-
                            #stock-recruit
                            c[1]#+c[2]*adultDATA[j]
                            
                            #environmental covariates
                            +beta[1]*covariates[j,1]+beta[7]*covariates[j,1]^2
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[2]*covariates[j,2]#winter flows>0.99
                            +beta[3]*covariates[j,3]#spring flows>0.99
                            +beta[4]*covariates[j,4]#summer flows<XX
                            +beta[5]*covariates[j,5]#summer temp>20C
                            +beta[6]*covariates[j,4]*covariates[j,5]#summer low flow/temp interaction
                            +beta[8]*pow(covariates[j,2],2)
  
                yDATA[j]~dbin(p[j],N[j])
  
          }
        
    
    #Derived output
      
        
          for(j in 1:nsamples){
                   yExp[j]<- N[j]*p[j]
                   E[j]<-pow((yDATA[j]-yExp[j]),2)/
                                    (yExp[j]+0.5)
                  
                   NNew[j]~dpois(lambda[j]*adultDATA[j])
                   yNew[j]~dbin(p[j],NNew[j])
                   ENew[j]<-pow((yNew[j]-yExp[j]),2)/
                                          (yExp[j]+0.5)
  
                   #error to compute rmse
                   Ey[j]<-pow(yDATA[j]-yExp[j],2)
                   pCheck1[j]<-yDATA[j] > yExp[j]
  
          }
          fit<-sum(E[])
          fitNew<-sum(ENew[])
          pCheck<-mean(pCheck1[])
        
        #rmse[sp]<-pow(sum(Ey[,])/(nsamples*4),-2)
      
  }",file="~/trout_yoy/abundanceModel.txt")