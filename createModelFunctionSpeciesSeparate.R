createModel<-function(stockRecruit=T,env='none',randomYear=F){
if((stockRecruit|env!='none') & randomYear){
  stop("not set up to have both random year and covariates")
}
 
if(stockRecruit & env=='mean'){
  cat("model{
    #Priors
    #State process
        for(r in 1:4){
          for(f in 1:8){
  #covariate betas
            beta[r,f]~dnorm(0,0.0001)
          }
  #stock-recruit(c[,,1] and c[,,2])/competition parameters(c[,,3])
          c[r,1]~dunif(0,50)
          c[r,2]~dnorm(0,0.001)
          c[r,3]~dnorm(0,0.001)
        }
  
    #Likelihood for yoy
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-
                            #stock-recruit
                            c[r,1]+c[r,2]*adultDATA[j,r]+c[r,3]*pow(adultDATA[j,r],2)
  
                            #environmental covariates
                            +beta[r,1]*covariates[j,r,1]+beta[r,7]*covariates[j,r,1]^2
  
                            #mean covariates
                          +beta[r,2]*covariates[j,r,6]#winter mean flow
                          +beta[r,3]*covariates[j,r,7]#spring mean flow
                          +beta[r,4]*covariates[j,r,8]#summer
                          +beta[r,5]*covariates[j,1,9]#summer mean flow and temp
                          +beta[r,6]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction
                          +beta[r,8]*pow(covariates[j,r,6],2)
  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }
        #rmse<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='mean'){
  cat("model{
    #Priors
    #State process
      
        for(r in 1:4){
        c[r,1]~dunif(0,50)
          for(f in 1:8){
  #covariate betas
            beta[r,f]~dnorm(0,0.0001)
          }
        }
  
    #Likelihood for yoy
      
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-c[r,1]

                            #environmental covariates
                          +beta[r,1]*covariates[j,r,1]+beta[r,7]*covariates[j,r,1]^2
  
                            #mean covariates
                          +beta[r,2]*covariates[j,r,6]
                          +beta[r,3]*covariates[j,r,7]#winter and spring mean flow
                          +beta[r,4]*covariates[j,r,8]
                          +beta[r,5]*covariates[j,1,9]#summer mean flow and temp
                          +beta[r,6]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction
                          +beta[r,8]*pow(covariates[j,r,6],2)
  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }
        #rmse<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
      
        for(r in 1:4){
          for(f in 1:8){
  #covariate betas
            beta[r,f]~dnorm(0,0.0001)
          }
  #stock-recruit(c[,,1] and c[,,2])/competition parameters(c[,,3])
          c[r,1]~dunif(0,50)
          c[r,2]~dnorm(0,0.0001)
          c[r,3]~dnorm(0,0.0001)
        }
  
    #Likelihood for yoy
      
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-
                            #stock-recruit
                            c[r,1]+c[r,2]*adultDATA[j,r]+c[r,3]*pow(adultDATA[j,r],2)
                            
                            #environmental covariates
                            +beta[r,1]*covariates[j,r,1]+beta[r,7]*covariates[j,r,1]^2
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[r,2]*covariates[j,r,2]#winter flows>0.99
                            +beta[r,3]*covariates[j,r,3]#spring flows>0.99
                            +beta[r,4]*covariates[j,r,4]#summer flows<XX
                            +beta[r,5]*covariates[j,1,5]#summer temp>20C
                            +beta[r,6]*covariates[j,r,4]*covariates[j,r,5]#summer low flow/temp interaction
                            +beta[r,8]*pow(covariates[j,r,2],2)
  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
  
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }
        #rmse[sp]<-pow(sum(Ey[,])/(nsamples*4),-2)
      
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
      
        for(r in 1:4){
        c[r,1]~dunif(0,50)
          for(f in 1:8){
  #covariate betas
            beta[r,f]~dnorm(0,0.0001)
          }
        }
  
    #Likelihood for yoy
      
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-c[r,1]

                            #environmental covariates
                            +beta[r,1]*covariates[j,r,1]+beta[r,7]*covariates[j,r,1]^2
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[r,2]*covariates[j,r,2]#winter flows>0.99
                            +beta[r,3]*covariates[j,r,3]#spring flows>0.99
                            +beta[r,4]*covariates[j,r,4]#summer flows<XX
                            +beta[r,5]*covariates[j,1,5]#summer temp>20C
                            +beta[r,6]*covariates[j,r,4]*covariates[j,r,5]#summer low flow/temp interaction
                            +beta[r,8]*pow(covariates[j,r,2],2)
  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
  
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }
        #rmse[sp]<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='none'){
    cat("model{
    #Priors
    #State process
      
        for(r in 1:4){
  
      #stock-recruit
          c[r,1]~dunif(0,50)
          c[r,2]~dnorm(0,0.001)
          c[r,3]~dnorm(0,0.001)
        }
      
    
    #Likelihood for yoy
      
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-
                            #stock-recruit
                            c[r,1]+c[r,2]*adultDATA[j,r]+c[r,3]*pow(adultDATA[j,r],2)
  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
      
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
  
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }

      
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      
        for(r in 1:4){
          for(f in 1:5){
  #covariate betas
            beta[r,f]~dnorm(0,0.01)
          }
  #stock-recruit
          c[r,1]~dunif(0,50)
          c[r,2]~dnorm(0,0.001)
          c[r,3]~dnorm(0,0.001)
        }

  
    #Likelihood for yoy
      
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-
                            #stock-recruit
                            c[r,1]+c[r,2]*adultDATA[j,r]+c[r,3]*pow(adultDATA[j,r],2)
  
                            #environmental covariates (based on pca in west brook)
                            +beta[r,1]*covariates[j,r,1]
                            +beta[r,2]*covariates[j,r,2]
                            +beta[r,3]*covariates[j,r,3]
                            +beta[r,4]*covariates[j,r,4]
                            +beta[r,5]*covariates[j,r,5]
#                             +beta[r,6]*covariates[j,r,6]
#                             +beta[r,7]*covariates[j,r,7]
#                             +beta[r,8]*covariates[j,r,8]
#                             +beta[r,9]*covariates[j,r,9]

  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }
        #rmse[sp]<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      
        for(r in 1:4){
          c[r,1]~dunif(0,50)
          for(f in 1:5){
  #covariate betas
            beta[r,f]~dnorm(0,0.0001)
          }
        }
  
    #Likelihood for yoy
      
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                log(lambda[j,r])<-c[r,1]

                            #environmental covariates (based on pca in west brook)
                            +beta[r,1]*covariates[j,r,1]
                            +beta[r,2]*covariates[j,r,2]
                            +beta[r,3]*covariates[j,r,3]
                            +beta[r,4]*covariates[j,r,4]
                            +beta[r,5]*covariates[j,r,5]
#                             +beta[r,6]*covariates[j,r,6]
#                             +beta[r,7]*covariates[j,r,7]
#                             +beta[r,8]*covariates[j,r,8]
#                             +beta[r,9]*covariates[j,r,9]

  
                yDATA[j,r]~dbin(p[j,r],N[j,r])
  
          }
        }
    
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }
        #rmse[sp]<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(randomYear){
  cat("model{
  #Priors
  
  #Random year effect

      
      for(r in 1:4){
        c[r,1]~dunif(0,50)
        for(j in 1:nsamples){
          eps[j,r]~dnorm(0,tau.eps[r])
        }
        tau.eps[r]<-pow(sd.eps[r],-2)
        sd.eps[r]~dunif(0,10)
      }

  #Likelihood for yoy
    
      for(r in 1:4){
        for(j in 1:nsamples){
            N[j,r]~dpois(lambda[j,r])
              log(lambda[j,r])<-c[r,1]+eps[j,r]

              yDATA[j,r]~dbin(p[j,r],N[j,r])

        }
      }
  
    #Derived output
      
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r]<- N[j,r]*p[j,r]
                   E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
                                    (yExp[j,r]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
                   yNew[j,r]~dbin(p[j,r],NNew[j,r])
                   ENew[j,r]<-pow((yNew[j,r]-yExp[j,r]),2)/
                                          (yExp[j,r]+0.5)
  
                   #error to compute rmse
                   Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
                   pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
          }
          fit[r]<-sum(E[,r])
          fitNew[r]<-sum(ENew[,r])
          pCheck[r]<-mean(pCheck1[,r])
        }

}",file="~/trout_yoy/abundanceModel.txt")
}

}