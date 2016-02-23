#split into term for stock and term for env combined after transformation by addition

createModel<-function(stockRecruit=T,env='none',randomYear=F){
if((stockRecruit|env!='none') & randomYear){
  stop("not set up to have both random year and covariates")
}
 
if(stockRecruit & env=='mean'){
  cat("model{
#Priors
#State process
  for(r in 1:nRivers){
    for(b in 1:8){
      beta[b,r]~dnorm(0,0.001)T(-3,3)
    }
    c[1,r]~dunif(0,3)
    c[2,r]~dnorm(0,0.001)T(-0.1,0.1)
    c[3,r]~dnorm(0,0.01)
  }

#Likelihood for yoy
for(r in 1:nRivers){
  for(j in 1:nYears){
    N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
    log(lambda[j,r])<-
      #stock-recruit
      c[1,r]+c[2,r]*adultDATA[j,r]+c[3,r]*summerAdultDATA[j,r]
  
      #environmental covariates
      +beta[1,r]*covariates[j,r,1]+beta[7,r]*pow(covariates[j,r,1],2)
  
      #mean covariates
      +beta[2,r]*(covariates[j,r,6])#winter mean flow
      +beta[3,r]*covariates[j,r,7]#spring mean flow
      +beta[4,r]*covariates[j,r,8]#summer
      +beta[5,r]*covariates[j,r,9]#summer mean flow and temp
      +beta[6,r]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction
      +beta[8,r]*pow(covariates[j,r,6],2)

    yDATA[j,r]~dbin(p[j,r],N[j,r])

    nPredict[j,r]<-lambda[j,r]*adultDATA[j,r]
  
    }
  }
    
    #Derived output
      
 for(r in 1:nRivers){       
  for(j in 1:nYears){
     #calculate expected values and chi-squared for observed
#      nExp[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
#      yExp[j,r]~dbin(p[j,r],nExp[j,r])
     yExp[j,r]<-lambda[j,r]*adultDATA[j,r]*p[j,r]
#      E[j,r]<-pow((yDATA[j,r]-yExp[j,r]),2)/
#                       (yExp[j,r]+0.5)
#   
#      #generate new data with estimated parameters  
#      NNew[j]~dpois(lambda[j]*adultDATA[j])
#      yNew[j]~dbin(p[j],NNew[j])
#      ENew[j]<-pow((yNew[j]-yExp[j]),2)/
#                             (yExp[j]+0.5)
#   
#      #error to compute rmse
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
#      pCheck1[j,r]<-yDATA[j,r] > yExp[j,r]
#   }
# }
# 
#   fit<-sum(E[])
#   fitNew<-sum(ENew[])
#   pCheck[r]<-mean(pCheck1[])
    }
    rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
  }
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='mean'){
  cat("model{
#Priors
#State process
  for(r in 1:nRivers){
    for(b in 1:8){
      beta[b,r]~dnorm(0,0.001)T(-3,3)
    }
    c[1,r]~dunif(0,5)
  }
      
#Likelihood for yoy
  for(r in 1:nRivers){
    for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r])
      log(lambda[j,r])<-
      #stock-recruit
      c[1,r]
      
      #environmental covariates
      +beta[1,r]*covariates[j,r,1]+beta[7,r]*pow(covariates[j,r,1],2)
      
      #mean covariates
      +beta[2,r]*(covariates[j,r,6])#winter mean flow
      +beta[3,r]*covariates[j,r,7]#spring mean flow
      +beta[4,r]*covariates[j,r,8]#summer
      +beta[5,r]*covariates[j,r,9]#summer mean flow and temp
      +beta[6,r]*covariates[j,r,8]*covariates[j,r,9]#summer mean flow/temp interaction
      +beta[8,r]*pow(covariates[j,r,6],2)
      
      yDATA[j,r]~dbin(p[j,r],N[j,r])
      
      nPredict[j,r]<-lambda[j,r]
      
    }
  }

 for(r in 1:nRivers){       
  for(j in 1:nYears){
      yExp[j,r]<-lambda[j,r]*p[j,r]
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
}
  rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
}

}",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='extreme'){
  cat("model{
    #Priors
    #State process
#Priors
#State process
  for(r in 1:nRivers){
  #environmental betas
    for(b in 1:8){
      beta[b,r]~dnorm(0,0.01)T(-3,3)
    }
  #stock recruit
    c[1,r]~dunif(0,3)
    c[2,r]~dnorm(0,0.01)T(-0.1,0.1)
    c[3,r]~dnorm(0,0.01)
    }
    #Likelihood for yoy
 #Likelihood for yoy
  for(r in 1:nRivers){
      for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
      log(lambda[j,r])<-
      #stock-recruit
      c[1,r]
      +c[2,r]*adultDATA[j,r]
      +c[3,r]*summerAdultDATA[j,r]
      
      #environmental covariates
      #+beta[1,r]*covariates[j,r,1]+beta[7,r]*pow(covariates[j,r,1],2)
      #+beta[1,r]*covariates[j,r,11]+beta[7,r]*covariates[j,r,12]
      
      #mean covariates
      +beta[2,r]*(covariates[j,r,2])#winter high flow
      +beta[3,r]*covariates[j,r,3]#spring floods
      +beta[4,r]*covariates[j,r,4]#summer low flow
      +beta[5,r]*covariates[j,r,5]#summer high temps
      +beta[6,r]*covariates[j,r,5]*covariates[j,r,4]#summer mean flow/temp interaction
      +beta[8,r]*covariates[j,r,10] #winter low flow
      
      yDATA[j,r]~dbin(p[j,r],N[j,r])
      
      nPredict[j,r]<-lambda[j,r]*adultDATA[j,r]
      
    }
  }  

 for(r in 1:nRivers){       
  for(j in 1:nYears){
      yExp[j,r]<-lambda[j,r]*adultDATA[j,r]*p[j,r]
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
}
  rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
} 
     
  }",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='extreme'){
  cat("model{
#Priors
#State process
  for(r in 1:nRivers){
    for(b in 1:8){
      beta[b,r]~dnorm(0,0.001)T(-3,3)
    }
    c[1,r]~dunif(0,5)
  }
      
#Likelihood for yoy
  for(r in 1:nRivers){
    for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r])
      log(lambda[j,r])<-
      #stock-recruit
        c[1,r]
      
      #environmental covariates
      +beta[1,r]*covariates[j,r,1]+beta[7,r]*pow(covariates[j,r,1],2)
      #+beta[1,r]*covariates[j,r,11]+beta[7,r]*covariates[j,r,12]
      
      #mean covariates
      +beta[2,r]*(covariates[j,r,2])#winter high flow
      +beta[3,r]*covariates[j,r,3]#spring floods
      +beta[4,r]*covariates[j,r,4]#summer low flow
      +beta[5,r]*covariates[j,r,5]#summer high temps
      +beta[6,r]*covariates[j,r,5]*covariates[j,r,4]#summer mean flow/temp interaction
      +beta[8,r]*covariates[j,r,10] #winter low flow
      
      yDATA[j,r]~dbin(p[j,r],N[j,r])
      
      nPredict[j,r]<-lambda[j,r]
      
    }
  } 

 for(r in 1:nRivers){       
  for(j in 1:nYears){
      yExp[j,r]<-lambda[j,r]*p[j,r]
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
}
  rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
}
}",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='none'){
    cat("model{
    #Priors
    #State process
      
  for(r in 1:nRivers){      
      #stock-recruit
          c[1,r]~dunif(0,3)
          c[2,r]~dnorm(0,0.001)T(-0.1,0.1)
  }
    
    #Likelihood for yoy
  for(r in 1:nRivers){
    for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
      log(lambda[j,r])<-
      #stock-recruit
        c[1,r]+c[2,r]*adultDATA[j,r]
  
      yDATA[j,r]~dbin(p[j,r],N[j,r])

      nPredict[j,r]<-lambda[j,r]*adultDATA[j,r]
  
    }
  }    

 for(r in 1:nRivers){       
        for(j in 1:nYears){
        yExp[j,r]<-lambda[j,r]*adultDATA[j,r]*p[j,r]
        Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
        }
        rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
}
}",file="~/trout_yoy/abundanceModel.txt")
}

if(stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
      
  for(r in 1:nRivers){      
    for(f in 1:5){
  #covariate betas
      beta[f,r]~dnorm(0,0.01)
    }
  #stock-recruit
    c[1,r]~dunif(0,3)
    c[2,r]~dnorm(0,0.001)T(-0.1,0.1)
  
    #Likelihood for yoy
      
  for(r in 1:nRivers)      
    for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r]*adultDATA[j,r])
      log(lambda[j,r])<-
      #stock-recruit
        c[1,r]+c[2,r]*adultDATA[j,r]

      #environmental covariates (based on pca in west brook)
        +beta[1,r]*covariates[j,r,1]
        +beta[2,r]*covariates[j,r,2]
        +beta[3,r]*covariates[j,r,3]
        +beta[4,r]*covariates[j,r,4]
        +beta[5,r]*covariates[j,r,5]
  
      yDATA[j,r]~dbin(p[j,r],N[j,r])

      nPredict[j,r]<-lambda[j,r]*adultDATA[j,r]
    }
  }

 for(r in 1:nRivers){       
      for(j in 1:nYears){
      yExp[j,r]<-lambda[j,r]*adultDATA[j,r]*p[j,r]
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
      }
      rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
}
}",file="~/trout_yoy/abundanceModel.txt")
}

if(!stockRecruit & env=='pca'){
  cat("model{
    #Priors
    #State process
  for(r in 1:nRivers){  
    c[1,r]~dunif(0,5)
    for(f in 1:5){
      #covariate betas
      beta[f,r]~dnorm(0,0.01)
    }
  }
         
  for(r in 1:nRivers)      
    for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r])
      log(lambda[j,r])<-
      #stock-recruit
        c[1,r]+
      
      #environmental covariates (based on pca in west brook)
        +beta[1,r]*covariates[j,r,1]
        +beta[2,r]*covariates[j,r,2]
        +beta[3,r]*covariates[j,r,3]
        +beta[4,r]*covariates[j,r,4]
        +beta[5,r]*covariates[j,r,5]
      
      yDATA[j,r]~dbin(p[j,r],N[j,r])
      
      nPredict[j,r]<-lambda[j,r]*adultDATA[j,r]
    }
  } 

 for(r in 1:nRivers){       
  for(j in 1:nYears){
      yExp[j,r]<-lambda[j,r]*p[j,r]
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
}
  rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
}
}",file="~/trout_yoy/abundanceModel.txt")
}

if(randomYear){
  cat("model{
  #Priors
  
  #Random year effect
  for(r in 1:nRivers){
    c[1,r]~dunif(0,5)
    for(j in 1:nYears){
      eps[j,r]~dnorm(0,tau.eps[r])
    }
    tau.eps[r]<-pow(sd.eps[r],-2)
    sd.eps[r]~dunif(0,10)
  }
      

  #Likelihood for yoy
  for(r in 1:nRivers){
    for(j in 1:nYears){
      N[j,r]~dpois(lambda[j,r])
      log(lambda[j,r])<-c[1,r]+eps[j,r]

      yDATA[j,r]~dbin(p[j,r],N[j,r])
      nPredict[j,r]<-lambda[j,r]
    }
  }
  
 for(r in 1:nRivers){       
  for(j in 1:nYears){
      yExp[j,r]<-lambda[j,r]*p[j,r]
      Ey[j,r]<-pow(yDATA[j,r]-yExp[j,r],2)
}
  rmse[r]<-pow(sum(Ey[,r])/nYears,0.5)
}
}",file="~/trout_yoy/abundanceModel.txt")
}
}