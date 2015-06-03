#sp = species
#s = site
#p = pass
#ob= observer
cat("model{
#Likelihood
  for(s in 1:nSites){
    N[s,1]~dpois(lambda[s,1])
    log(lambda[s,1])<-beta[1,1]*envDATA[s,1]+beta[2,1]*envDATA[s,2]

    for(sp in 2:nSpecies){
      N[s,sp]~dpois(lambda[s,sp])
      log(lambda[s,sp])<-beta[1,sp]*envDATA[s,1]+beta[2,sp]*envDATA[s,2]
                        +beta[3,sp]*N[s,1]
    }

  for(sp in 1:3){
    for(p in 1:nPass){
      yDATA[s,sp,p]~dbin(pObs,N[s,sp])
    }
  }
}

#Priors
for(sp in 1:nSpecies){
    for(b in 1:3){
      beta[b,sp]~dnorm(0,0.01)
  }
}

pObs~dunif(0,1)

}",file="model.txt")
###is it possible to account for overestimates, which are a result of 
# counting fish outside the quad, by using an inclusion probability
# so it would be two binomials, one with the detection probability, and
# a second given detection the probability that it is in the quadrat
# Might not be able to parse detection from inclusion

# need to see what the exact method for estimating capture 
# probability from spatially replicated samples is before being sure
