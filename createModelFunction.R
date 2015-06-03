createModel<-function(stockRecruit=T,env='none',randomYear=F){
if((stockRecruit|env!='none') & randomYear){
  stop("not set up to have both random year and covariates")
}
 
if(stockRecruit & env=='mean'){
  cat("model{
    #Priors
    #State process
      for(sp in 1:2){
        for(r in 1:4){
          for(f in 1:8){
  #covariate betas
            beta[r,f,sp]~dnorm(0,0.0001)
          }
  #stock-recruit(c[,,1] and c[,,2])/competition parameters(c[,,3])
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
  
                            #environmental covariates
                            +beta[r,1,sp]*covariates[j,r,1]#+beta[r,7,sp]*covariates[j,r,1]^2
  
                            #mean covariates
                          +beta[r,2,sp]*covariates[j,r,6]#winter mean flow
                          +beta[r,3,sp]*covariates[j,r,7]#spring mean flow
                          +beta[r,4,sp]*covariates[j,r,8]#summer
                          +beta[r,5,sp]*covariates[j,1,9]#summer mean flow and temp
                          #+beta[r,6,sp]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction
                          +beta[r,8,sp]*pow(covariates[j,r,6],2)
  
                yDATA[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])
  
          }
        }
      }
    
    #Derived output
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r,sp]<- N[j,r,sp]*p[j,r,sp]
                   E[j,r,sp]<-pow((yDATA[j,r,sp]-yExp[j,r,sp]),2)/
                                    (yExp[j,r,sp]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                   yNew[j,r,sp]~dbin(p[j,r,sp],NNew[j,r,sp])
                   ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                          (yExp[j,r,sp]+0.5)
  
                   #error to compute rmse
                   Ey[j,r,sp]<-pow(yDATA[j,r,sp]-yExp[j,r,sp],2)
          }
          fit[r,sp]<-sum(E[,r,sp])
          fitNew[r,sp]<-sum(ENew[,r,sp])
        }
        rmse[sp]<-pow(sum(Ey[,,sp])/(nsamples*4),-2)
      }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='mean'){
  cat("model{
    #Priors
    #State process
      for(sp in 1:2){
        for(r in 1:4){
        c[r,sp,1]~dunif(0,50)
          for(f in 1:8){
  #covariate betas
            beta[r,f,sp]~dnorm(0,0.0001)
          }
        }
      }
  
    #Likelihood for yoy
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                log(lambda[j,r,sp])<-c[r,sp,1]

                            #environmental covariates
                          +beta[r,1,sp]*covariates[j,r,1]#+beta[r,7,sp]*covariates[j,r,1]^2
  
                            #mean covariates
                          +beta[r,2,sp]*covariates[j,r,6]
                          +beta[r,3,sp]*covariates[j,r,7]#winter and spring mean flow
                          +beta[r,4,sp]*covariates[j,r,8]
                          +beta[r,5,sp]*covariates[j,1,9]#summer mean flow and temp
                          #+beta[r,6,sp]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction
                          +beta[r,8,sp]*pow(covariates[j,r,6],2)
  
                yDATA[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])
  
          }
        }
      }
    
    #Derived output
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r,sp]<- N[j,r,sp]*p[j,r,sp]
                   E[j,r,sp]<-pow((yDATA[j,r,sp]-yExp[j,r,sp]),2)/
                                    (yExp[j,r,sp]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                   yNew[j,r,sp]~dbin(p[j,r,sp],NNew[j,r,sp])
                   ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                          (yExp[j,r,sp]+0.5)
  
                   #error to compute rmse
                   Ey[j,r,sp]<-pow(yDATA[j,r,sp]-yExp[j,r,sp],2)
          }
          fit[r,sp]<-sum(E[,r,sp])
          fitNew[r,sp]<-sum(ENew[,r,sp])
        }
        rmse[sp]<-pow(sum(Ey[,,sp])/(nsamples*4),-2)
      }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
      for(sp in 1:2){
        for(r in 1:4){
          for(f in 1:8){
  #covariate betas
            beta[r,f,sp]~dnorm(0,0.0001)
          }
  #stock-recruit(c[,,1] and c[,,2])/competition parameters(c[,,3])
          c[r,sp,1]~dunif(0,50)
          c[r,sp,2]~dnorm(0,0.0001)
          c[r,sp,3]~dnorm(0,0.0001)
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
                            
                            #environmental covariates
                            +beta[r,1,sp]*covariates[j,r,1]#+beta[r,7,sp]*covariates[j,r,1]^2
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[r,2,sp]*covariates[j,r,2]#winter flows>0.99
                            +beta[r,3,sp]*covariates[j,r,3]#spring flows>0.99
                            +beta[r,4,sp]*covariates[j,r,4]#summer flows<XX
                            +beta[r,5,sp]*covariates[j,1,5]#summer temp>20C
                            #+beta[r,6,sp]*covariates[j,r,4]*covariates[j,r,5]#summer low flow/temp interaction
                            +beta[r,8,sp]*pow(covariates[j,r,2],2)
  
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
                  
                   NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                   yNew[j,r,sp]~dbin(p[j,r,sp],NNew[j,r,sp])
                   ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                          (yExp[j,r,sp]+0.5)
  
                   #error to compute rmse
                   Ey[j,r,sp]<-pow(yDATA[j,r,sp]-yExp[j,r,sp],2)
  
          }
          fit[r,sp]<-sum(E[,r,sp])
          fitNew[r,sp]<-sum(ENew[,r,sp])
        }
        rmse[sp]<-pow(sum(Ey[,,sp])/(nsamples*4),-2)
      }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
      for(sp in 1:2){
        for(r in 1:4){
        c[r,sp,1]~dunif(0,50)
          for(f in 1:8){
  #covariate betas
            beta[r,f,sp]~dnorm(0,0.0001)
          }
        }
      }
  
    #Likelihood for yoy
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                log(lambda[j,r,sp])<-c[r,sp,1]

                            #environmental covariates
                            +beta[r,1,sp]*covariates[j,r,1]#+beta[r,7,sp]*covariates[j,r,1]^2
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[r,2,sp]*covariates[j,r,2]#winter flows>0.99
                            +beta[r,3,sp]*covariates[j,r,3]#spring flows>0.99
                            +beta[r,4,sp]*covariates[j,r,4]#summer flows<XX
                            +beta[r,5,sp]*covariates[j,1,5]#summer temp>20C
                            #+beta[r,6,sp]*covariates[j,r,4]*covariates[j,r,5]#summer low flow/temp interaction
                            +beta[r,8,sp]*pow(covariates[j,r,2],2)
  
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
                  
                   NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                   yNew[j,r,sp]~dbin(p[j,r,sp],NNew[j,r,sp])
                   ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                          (yExp[j,r,sp]+0.5)
  
                   #error to compute rmse
                   Ey[j,r,sp]<-pow(yDATA[j,r,sp]-yExp[j,r,sp],2)
  
          }
          fit[r,sp]<-sum(E[,r,sp])
          fitNew[r,sp]<-sum(ENew[,r,sp])
        }
        rmse[sp]<-pow(sum(Ey[,,sp])/(nsamples*4),-2)
      }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='none'){
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
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      for(sp in 1:2){
        for(r in 1:4){
          for(f in 1:9){
  #covariate betas
            beta[r,f,sp]~dnorm(0,0.01)
          }
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
  
                            #environmental covariates (based on pca in west brook)
                            +beta[r,1,sp]*covariates[j,r,1]
                            +beta[r,2,sp]*covariates[j,r,2]
                            +beta[r,3,sp]*covariates[j,r,3]
                            +beta[r,4,sp]*covariates[j,r,4]
                            +beta[r,5,sp]*covariates[j,r,5]
#                             +beta[r,6,sp]*covariates[j,r,6]
#                             +beta[r,7,sp]*covariates[j,r,7]
#                             +beta[r,8,sp]*covariates[j,r,8]
#                             +beta[r,9,sp]*covariates[j,r,9]

  
                yDATA[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])
  
          }
        }
      }
    
    #Derived output
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r,sp]<- N[j,r,sp]*p[j,r,sp]
                   E[j,r,sp]<-pow((yDATA[j,r,sp]-yExp[j,r,sp]),2)/
                                    (yExp[j,r,sp]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                   yNew[j,r,sp]~dbin(p[j,r,sp],NNew[j,r,sp])
                   ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                          (yExp[j,r,sp]+0.5)
  
                   #error to compute rmse
                   Ey[j,r,sp]<-pow(yDATA[j,r,sp]-yExp[j,r,sp],2)
          }
          fit[r,sp]<-sum(E[,r,sp])
          fitNew[r,sp]<-sum(ENew[,r,sp])
        }
        rmse[sp]<-pow(sum(Ey[,,sp])/(nsamples*4),-2)
      }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      for(sp in 1:2){
        for(r in 1:4){
          c[r,sp,1]~dunif(0,50)
          for(f in 1:9){
  #covariate betas
            beta[r,f,sp]~dnorm(0,0.0001)
          }
        }
      }
  
    #Likelihood for yoy
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
              N[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                log(lambda[j,r,sp])<-c[r,sp,1]

                            #environmental covariates (based on pca in west brook)
                            +beta[r,1,sp]*covariates[j,r,1]
                            +beta[r,2,sp]*covariates[j,r,2]
                            +beta[r,3,sp]*covariates[j,r,3]
                            +beta[r,4,sp]*covariates[j,r,4]
                            +beta[r,5,sp]*covariates[j,r,5]
#                             +beta[r,6,sp]*covariates[j,r,6]
#                             +beta[r,7,sp]*covariates[j,r,7]
#                             +beta[r,8,sp]*covariates[j,r,8]
#                             +beta[r,9,sp]*covariates[j,r,9]

  
                yDATA[j,r,sp]~dbin(p[j,r,sp],N[j,r,sp])
  
          }
        }
      }
    
    #Derived output
      for(sp in 1:2){
        for(r in 1:4){
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j,r,sp]<- N[j,r,sp]*p[j,r,sp]
                   E[j,r,sp]<-pow((yDATA[j,r,sp]-yExp[j,r,sp]),2)/
                                    (yExp[j,r,sp]+0.5)
                
                   #generate new data with estimated parameters  
                   NNew[j,r,sp]~dpois(lambda[j,r,sp]*adultDATA[j,r,sp])
                   yNew[j,r,sp]~dbin(p[j,r,sp],NNew[j,r,sp])
                   ENew[j,r,sp]<-pow((yNew[j,r,sp]-yExp[j,r,sp]),2)/
                                          (yExp[j,r,sp]+0.5)
  
                   #error to compute rmse
                   Ey[j,r,sp]<-pow(yDATA[j,r,sp]-yExp[j,r,sp],2)
          }
          fit[r,sp]<-sum(E[,r,sp])
          fitNew[r,sp]<-sum(ENew[,r,sp])
        }
        rmse[sp]<-pow(sum(Ey[,,sp])/(nsamples*4),-2)
      }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(randomYear){
  cat("model{
  #Priors
  
  #Random year effect

    for(sp in 1:2){  
      for(r in 1:4){
        c[r,sp,1]~dunif(0,50)
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
              log(lambda[j,r,sp])<-c[r,sp,1]+eps[j,r,sp]

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
}",file="~/trout_yoy/abundanceModel.txt")
}

}