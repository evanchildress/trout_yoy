createModel<-function(stockRecruit=T,env='none',randomYear=F){
if((stockRecruit|env!='none') & randomYear){
  stop("not set up to have both random year and covariates")
}
 
if(stockRecruit & env=='mean'){
  cat("model{
    #Priors
    #State process
        
          for(f in 1:6){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
      ",
      ifelse(species=='bkt'&r==4,
        "beta[8]~dnorm(0,0.01)",
        "beta[8]~dnorm(0,0.01)"),
  "
      beta[7]~dnorm(0,0.01)

  #stock-recruit
          c[1]~dunif(0,50)
          c[2]~dnorm(0,0.01)
          #c[3]~dnorm(0,0.01)
  
    #Likelihood for yoy
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-
                            #stock-recruit
                            c[1]+c[2]*adultDATA[j]#+c[3]*pow(adultDATA[j],2)
  
                            #environmental covariates
                            +beta[1]*covariates[j,1]+beta[7]*pow(covariates[j,1],2)
  
                            #mean covariates
                          +beta[2]*(covariates[j,6])#winter mean flow
                          +beta[3]*covariates[j,7]#spring mean flow
                          +beta[4]*covariates[j,8]#summer
                          +beta[5]*covariates[j,9]#summer mean flow and temp
                          +beta[6]*covariates[j,8]*covariates[j,9]#summer mean flow/temp interaction
                          +beta[8]*pow(covariates[j,6],2)
  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
          }
    
    #Derived output
      
        
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j]<- N[j]*p[j]
                   E[j]<-pow((yDATA[j]-yExp[j]),2)/
                                    (yExp[j]+0.5)
                
                   #generate new data with estimated parameters  
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

        #rmse<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='mean'){
  cat("model{
    #Priors
    #State process
      
        
        c[1]~dunif(0,50)

          for(f in 1:6){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
      ",
      ifelse(species=='bkt'&r==4,
        "beta[8]~dnorm(0,0.01)",
        "beta[8]~dnorm(0,0.01)"),
  "
      beta[7]~dnorm(0,0.01)


  
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-c[1]

                            #environmental covariates
                          +beta[1]*covariates[j,1]+beta[7]*pow(covariates[j,1],2)
  
                            #mean covariates
                          +beta[2]*covariates[j,6]
                          +beta[3]*covariates[j,7]#winter and spring mean flow
                          +beta[4]*covariates[j,8]
                          +beta[5]*covariates[j,9]#summer mean flow and temp
                          +beta[6]*covariates[j,8]*covariates[j,9]#summer mean flow/temp interaction
                          +beta[8]*pow(covariates[j,6],2)
  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
          }
        
    
    #Derived output
      
        
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j]<- N[j]*p[j]
                   E[j]<-pow((yDATA[j]-yExp[j]),2)/
                                    (yExp[j]+0.5)
                
                   #generate new data with estimated parameters  
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
        
        #rmse<-pow(sum(Ey[,])/(nsamples*4),-2)
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
      
          for(f in 1:6){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
      ",
      ifelse(species=='bkt'&r==4,
        "beta[8]~dnorm(0,0.01)",
        "beta[8]~dnorm(0,0.01)"),
  "
      beta[7]~dnorm(0,0.01)

  #stock-recruit
          c[1]~dunif(0,50)
          c[2]~dnorm(0,0.01)
          #c[3]~dnorm(0,0.01)
        
  
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-
                            #stock-recruit
                            c[1]+c[2]*adultDATA[j]#+pow(adultDATA[j],2)*c[3]#+c[3]*otherSpDATA[j]
                            
                            #environmental covariates
                            +beta[1]*covariates[j,1]+beta[7]*pow(covariates[j,1],2)
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[2]*(covariates[j,2])#winter flows>0.99
                            +beta[3]*covariates[j,3]#spring flows>0.99
                            +beta[4]*covariates[j,4]#summer flows<XX
                            +beta[5]*covariates[j,5]#summer temp>20C
                            +beta[6]*covariates[j,8]*covariates[j,9]#summer mean flow/temp interaction
                            +beta[8]*pow(covariates[j,2],2)
  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
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
}

if(!stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
      
        
        c[1]~dunif(0,50)

  for(f in 1:6){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
      ",
      ifelse(species=='bkt'&r==4,
        "beta[8]~dnorm(0,0.01)",
        "beta[8]~dnorm(0,0.01)"),
  "
      beta[7]~dnorm(0,0.01)

  
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-c[1]

                            #environmental covariates
                            +beta[1]*covariates[j,1]+beta[7]*pow(covariates[j,1],2)
  
                            #extreme covariates (either this or the means one should be commented out)
                            +beta[2]*covariates[j,2]#winter flows>0.99
                            +beta[3]*covariates[j,3]#spring flows>0.99
                            +beta[4]*covariates[j,4]#summer flows<XX
                            +beta[5]*covariates[j,5]#summer temp>20C
                            +beta[6]*covariates[j,8]*covariates[j,9]#summer mean flow/temp interaction
                            +beta[8]*pow(covariates[j,2],2)
  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
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
}

if(stockRecruit & env=='none'){
    cat("model{
    #Priors
    #State process
      
        
  
      #stock-recruit
          c[1]~dunif(0,50)
          c[2]~dnorm(0,0.01)
          #c[3]~dnorm(0,0.01)
      
    
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-
                            #stock-recruit
                            c[1]+c[2]*adultDATA[j]#+c[3]*pow(adultDATA[j],2)
  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
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
        

      
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      
        
          for(f in 1:5){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
  #stock-recruit
          c[1]~dunif(0,50)
          c[2]~dnorm(0,0.01)
          #c[3]~dnorm(0,0.01)

  
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-
                            #stock-recruit
                            c[1]+c[2]*adultDATA[j]#+c[3]*pow(adultDATA[j],2)
  
                            #environmental covariates (based on pca in west brook)
                            +beta[1]*covariates[j,1]
                            +beta[2]*covariates[j,2]
                            +beta[3]*covariates[j,3]
                            +beta[4]*covariates[j,4]
                            +beta[5]*covariates[j,5]
#                             +beta[6]*covariates[j,6]
#                             +beta[7]*covariates[j,7]
#                             +beta[8]*covariates[j,8]
#                             +beta[9]*covariates[j,9]

  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
          }
         
    
    #Derived output
      
        
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j]<- N[j]*p[j]
                   E[j]<-pow((yDATA[j]-yExp[j]),2)/
                                    (yExp[j]+0.5)
                
                   #generate new data with estimated parameters  
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
}

if(!stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      
        
          c[1]~dunif(0,50)
          for(f in 1:5){
  #covariate betas
            beta[f]~dnorm(0,0.01)
          }
         
  
    #Likelihood for yoy
      
        
          for(j in 1:nsamples){
              N[j]~dpois(lambda[j]*adultDATA[j])
                log(lambda[j])<-c[1]

                            #environmental covariates (based on pca in west brook)
                            +beta[1]*covariates[j,1]
                            +beta[2]*covariates[j,2]
                            +beta[3]*covariates[j,3]
                            +beta[4]*covariates[j,4]
                            +beta[5]*covariates[j,5]
#                             +beta[6]*covariates[j,6]
#                             +beta[7]*covariates[j,7]
#                             +beta[8]*covariates[j,8]
#                             +beta[9]*covariates[j,9]

  
                yDATA[j]~dbin(p[j],N[j])

              nPredict[j]<-lambda[j]*adultDATA[j]
  
          }
        
    
    #Derived output
      
        
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j]<- N[j]*p[j]
                   E[j]<-pow((yDATA[j]-yExp[j]),2)/
                                    (yExp[j]+0.5)
                
                   #generate new data with estimated parameters  
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
}

if(randomYear){
  cat("model{
  #Priors
  
  #Random year effect

      
      
        c[1]~dunif(0,50)
        for(j in 1:nsamples){
          eps[j]~dnorm(0,tau.eps)
        }
        tau.eps<-pow(sd.eps,-2)
        sd.eps~dunif(0,10)
      

  #Likelihood for yoy
    
      
        for(j in 1:nsamples){
            N[j]~dpois(lambda[j])
              log(lambda[j])<-c[1]+eps[j]

              yDATA[j]~dbin(p[j],N[j])

            nPredict[j]<-lambda[j]

        }
      
  
    #Derived output
      
        
          for(j in 1:nsamples){
                   #calculate expected values and chi-squared for observed
                   yExp[j]<- N[j]*p[j]
                   E[j]<-pow((yDATA[j]-yExp[j]),2)/
                                    (yExp[j]+0.5)
                
                   #generate new data with estimated parameters  
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
        

}",file="~/trout_yoy/abundanceModel.txt")
}

}