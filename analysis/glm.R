library(glmulti)
load('~/trout_yoy/abundanceData.RData')
  sumTempMean<-data.table(covariates[,,9])
  sumTempFill<-lm(V4~V3+V2+V1,data=sumTempMean)
  covariates[6:9,4,9]<-predict(sumTempFill,data.frame(sumTempMean[6:9,list(V3,V2,V1)]))
  rm(sumTempMean)
  rm(sumTempFill)
  
  y<-y[3:15,,,1]
  
  covariates<-covariates[1:13,,]


popEst<-readRDS("~/trout_yoy/results/modelOutput/randomYear.rds")$BUGSoutput$sims.list$N
  popEst<-data.table(melt(popEst))
  setnames(popEst,c("sim","year","river","species","estimate"))
  setkey(popEst,species,year)
  popEst<-popEst[,list(round(mean(estimate),3),round(quantile(estimate,probs=0.025),3),
                   round(quantile(estimate,probs=0.975),3)),
             by=list(species,river,year)]
  setnames(popEst,c("V1","V2","V3"),c("mean","lower","upper"))
  popEst[,year:=year+2001]

env<-dcast(melt(covariates),Var1+Var2~Var3)
  env<-data.table(env)

adult<-data.table(melt(A))
  setnames(adult,c("year","river","species","adults"))
  adult[,year:=year+1]
  adult[,species:=as.integer(ifelse(species=='bkt',1,2))]
  adult<-adult[year<15]
  setnames(env,c("Var1","Var2"),c("year","river"))
  setkey(env,river,year)
  setkey(popEst,river,year,species)

popEnv<-env[popEst]
  setkey(popEst,river,year,species)
  setkey(adult,river,year,species)
  popEnv[,adults:=adult$adults]
  popEnv[,riverFactor:=as.factor(river)]

popEnv[,fall2:=fall_discharge^2]
  popEnv[,winter_mean2:=winter_flow_mean^2]
  popEnv[,winter_floods2:=winter_floods^2]
  popEnv[,adults2:=adults^2]
  popEnv[,mean:=round(mean)]

predictors<-c("spring_floods","fall_discharge",
        "winter_floods","summer_low","summer_temp",
        "fall2","winter_floods2",
        "adults","adults2")

mPredictors<-c("winter_flow_mean","spring_flow_mean","summer_flow_mean",
        "summer_temp_mean","winter_mean2","adults","adults2")

# for(r in 1:4){
#   for(sp in 1:2){
#     assign(paste0("glmulti",sp,r),
#            glmulti("mean",predictors,
#                    data=popEnv[river==r & species==sp],
#                    family=poisson,maxsize=12,level=1))
#   }
# }


covPCA<-array(NA,dim=dim(covariates))

for(r in 1:4){
    assign(paste0("pca",r),prcomp(covariates[,r,],center=T,scale=T))
    scores<-get(paste0("pca",r))$x
    covPCA[,r,]<-scores
}

covPCA<-data.table(dcast(melt(covPCA),Var1+Var2~Var3))
setnames(covPCA,c("year","river",paste0(rep("pc",9),1:9)))
covPCA[,year:=year+2001]
setkey(covPCA,river,year)
popEnv<-covPCA[popEnv]

glmAic<-array(NA,dim=c(7,4,2),
              dimnames=list(c("x","m","pca","sr","xSr","mSr","pcaSr"),
                            c("jimmy","mitchell","obear","wb"),
                            c("bkt","bnt")))


glmFun<-function(r,sp){  
    assign(paste0("x",sp,r),glm(round(mean)~spring_floods+
               fall_discharge+winter_floods+summer_low*summer_temp
               +fall2+winter_floods2+offset(adults),
             data=popEnv[river==r&species==sp],family=poisson))
    
    assign(paste0("m",sp,r),glm(round(mean)~spring_flow_mean+
               fall_discharge+winter_flow_mean+
               summer_flow_mean*summer_temp_mean+fall2+winter_mean2,
             data=popEnv[river==r&species==sp],family=poisson))

    assign(paste0("pca",sp,r),glm(round(mean)~pc1+pc2+pc3+pc4+pc5,
             data=popEnv[river==r&species==sp],family=poisson))

    assign(paste0("sr",sp,r),glm(round(mean)~adults+adults2,
             data=popEnv[river==r&species==sp],family=poisson))
    
    assign(paste0("xSr",sp,r),glm(round(mean)~spring_floods+
               fall_discharge+winter_floods+summer_low*summer_temp+
               fall2+winter_floods2+adults+adults2,
             data=popEnv[river==r&species==sp],family=poisson))
    
    assign(paste0("mSr",sp,r),glm(round(mean)~spring_flow_mean+
               fall_discharge+winter_flow_mean+
               summer_flow_mean*summer_temp_mean+fall2+
               winter_mean2+adults+adults2,
             data=popEnv[river==r&species==sp],family=poisson))
    
    assign(paste0("pcaSr",sp,r),glm(round(mean)~pc1+pc2+pc3+pc4+pc5+
                                      adults+adults2,
         data=popEnv[river==r&species==sp],family=poisson))

    a<-c(aicc(get(paste0("x",sp,r))),
         aicc(get(paste0("m",sp,r))),
         aicc(get(paste0("pca",sp,r))),
         aicc(get(paste0("sr",sp,r))),
         aicc(get(paste0("xSr",sp,r))),
         aicc(get(paste0("mSr",sp,r))),
         aicc(get(paste0("pcaSr",sp,r)))
         )
    return(a)
}

for(r in 1:4){
  for(sp in 1:2){
    glmAic[,r,sp]<-glmFun(r,sp)
  }
}

    assign(paste0("xSr",sp,r),glm(round(mean)~spring_floods+
               fall_discharge+winter_floods+
               fall2+winter_floods2+adults+adults2+adults*fall_discharge+adults*fall2,
             data=popEnv[river==r&species==sp],family=poisson))


  
poissonmfx(round(mean)~spring_floods+
               fall_discharge+winter_floods+summer_low*summer_temp+
               fall2+winter_floods2+adults+adults2,
           data=popEnv[species==1&river==1])

