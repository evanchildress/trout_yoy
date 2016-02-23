library(glmulti)
if(!exists("jagsData")){
  jagsData<-readRDS("cjsInputs/jagsData.rds")
}
riverNames<-data.table(number=1:4,name=c("west brook","wb jimmy","wb mitchell","wb obear"))

popEst<-readRDS("~/trout_yoy/results/modelOutput/randomYearBnt.rds")$BUGSoutput$sims.list$N %>%
  apply(c(2,3),mean) %>%
  melt() %>%
  data.table() %>%
  setnames(c("year","river","popEst")) %>%
  .[,year:=year+2001] %>%
  .[,river:=riverNames[match(river,number),name]] %>%
  setkey(river,year) %>%
  .[,logLambda:=log(popEst)]

env<-melt(jagsData$covariates) %>%
  dcast(Var1+Var2~Var3) %>%
  data.table() %>%
  setnames(c("Var1","Var2"),c("year","river")) %>%
  setkey(river,year) %>%
  .[,":="(fall2=fall_discharge^2,
          winter2=winter_flow_mean^2,
          summerIntExtreme=summer_low*summer_temp,
          summerIntMean=summer_flow_mean*summer_temp_mean)]

adult<-jagsData$adultDATA %>%
  melt() %>%
  data.table() %>%
  setnames(c("year","river","fallAdults")) %>%
  .[,year:=year+2000] %>%
  .[,river:=riverNames[match(river,number),name]] %>%
  setkey(river,year)

summerAdult<-jagsData$summerAdultDATA %>%
  melt() %>%
  data.table() %>%
  setnames(c("year","river","summerAdults")) %>%
  .[,year:=year+1999] %>%
  .[,river:=riverNames[match(river,number),name]] %>%
  setkey(river,year)

data<-popEst[env] %>%
      setkey(river,year) %>%
      .[adult] %>%
      setkey(river,year) %>%
      .[summerAdult] %>%
      .[,logPerCapita:=log(popEst/fallAdults)]

predictors<-c("spring_floods","fall_discharge",
        "winter_floods","summer_low","summer_temp",
        "fall2","fall_floods","fall_low","winter_low",
        "fallAdults","summerAdults","winter_flow_mean","spring_flow_mean","summer_flow_mean",
        "summer_temp_mean","winter2")

# for(r in 1:4){
#   for(sp in 1:2){
#     assign(paste0("glmulti",sp,r),
#            glmulti("mean",predictors,
#                    data=popEnv[river==r & species==sp],
#                    family=poisson,maxsize=12,level=1))
#   }
# }

# 
# covPCA<-array(NA,dim=dim(covariates))
# 
# for(r in 1:4){
#     assign(paste0("pca",r),prcomp(covariates[,r,],center=T,scale=T))
#     scores<-get(paste0("pca",r))$x
#     covPCA[,r,]<-scores
# }
# 
# covPCA<-data.table(dcast(melt(covPCA),Var1+Var2~Var3))
# setnames(covPCA,c("year","river",paste0(rep("pc",9),1:9)))
# covPCA[,year:=year+2001]
# setkey(covPCA,river,year)
# popEnv<-covPCA[popEnv]
# 
# glmAic<-array(NA,dim=c(7,4,2),
#               dimnames=list(c("x","m","pca","sr","xSr","mSr","pcaSr"),
#                             c("jimmy","mitchell","obear","wb"),
#                             c("bkt","bnt")))

# fit<-glm(log(popEst)~spring_floods+fall_discharge+winter_floods+
#             summer_low+summer_temp+fall2+fall_floods+
#             fall_low+winter_low+fallAdults+
#             summerAdults+winter_flow_mean+spring_flow_mean+
#             summer_flow_mean+summer_temp_mean+winter2,
#     data=data[river=="west brook"],family="gaussian")
# leaps<-regsubsets(logLambda~spring_floods+fall_discharge+winter_floods+
#                                 summer_low+summer_temp+fall2+fall_floods+
#                                 fall_low+winter_low+fallAdults+
#                                 summerAdults+winter_flow_mean+spring_flow_mean+
#                                 summer_flow_mean+summer_temp_mean+winter2,
                  
                  
predictorSets<-list(fall=list(extreme=c("fall_floods","fall_low"),
                              mean=c("fall_discharge","fall2"),
                              null=NULL),
                    winter=list(extreme=c("winter_floods","winter_low"),
                                mean=c("winter_flow_mean","winter2"),
                                null=NULL),
                    spring=list(extreme="spring_floods",
                                mean="spring_flow_mean",
                                null=NULL),
                    summer=list(extreme=c("summer_low","summer_temp"),
                                mean=c("summer_flow_mean","summer_temp_mean"),
                                extremeInteraction=c("summer_low","summer_temp","summerIntExtreme"),
                                meanInteraction=c("summer_flow_mean","summer_temp_mean","summerIntMean"),
                                null=NULL),
                    stock=list(fall="fallAdults",
                               summer="summerAdults",
                               both=c("fallAdults","summerAdults"),
                               null=NULL))
results<-array(NA,dim=c(540,8),dimnames=list(NULL,c(names(predictorSets),"aic","rmse","rsq")))
results<-list("west brook"=results,
              "wb jimmy"=results,
              "wb mitchell"=results,
              "wb obear"=results)
rmse<-function(model,standardize=T){
  result<-sqrt(sum(model$residuals^2)/length(model$residuals))
  if(standardize) result<-result/abs(mean(model$fitted.values+model$residuals))
  return(result)
}
row<-1
for(b1 in names(predictorSets$fall)){
  for(b2 in names(predictorSets$winter)){
    for(b3 in names(predictorSets$spring)){
      for(b4 in names(predictorSets$summer)){
        for(b5 in names(predictorSets$stock)){
          modelString<-paste(c(predictorSets$fall[[b1]],
                               predictorSets$winter[[b2]],
                               predictorSets$spring[[b3]],
                               predictorSets$summer[[b4]],
                               predictorSets$stock[[b5]]),
                             collapse="+")
          if(modelString=="") modelString<-"1"
          model<-formula(paste("logPerCapita~",modelString))
          for(r in unique(data$river)){
           fit<-lm(model,data=data[river==r])
           results[[r]][row,]<-c(b1,b2,b3,b4,b5,AIC(fit),rmse(fit),summary(fit)$r.squared)
          }
          row<-row+1
        }
      }
    }
  }
}


bestModels<-NULL
for(r in names(results)){
  bestModels<-data.table(results[[r]]) %>%
  .[,river:=r] %>%
  rbind(bestModels)
}

bestModels<-bestModels[,":="(aic=as.numeric(aic),
                             rmse=as.numeric(rmse),
                             rsq=as.numeric(rsq))] %>%
            .[,":="(minAic=min(aic),
                    minRmse=min(rmse),
                    maxRsq=max(rsq)),
              by=river] %>%
            .[aic-minAic<=2|rmse==minRmse|rsq==maxRsq] %>%
            .[,":="(minAic=NULL,
                    minRmse=NULL,
                    maxRsq=NULL)] %>%
            setkey(river,aic)

bestFits<-list()
for(r in unique(bestModels$river)){
  row<-bestModels[river==r][aic==min(aic)]
  modelString<-paste(c(predictorSets$fall[[row$fall]],
                       predictorSets$winter[[row$winter]],
                       predictorSets$spring[[row$spring]],
                       predictorSets$summer[[row$summer]],
                       predictorSets$stock[[row$stock]]),
                     collapse="+")
  model<-formula(paste("logPerCapita~",modelString))
  fit<-lm(model,data=data[river==r])
  bestFits[[r]]<-fit
}
extremeFits<-list()
meanFits<-list()
tiff.par("results/figures/nFromFrequentist.tif",mfrow=c(4,1))
for(r in unique(data$river)){
  plot(popEst~year,data=data[river==r],type='l')
  predictedN<-exp(bestFits[[r]]$fitted.values)*data[river==r,fallAdults]
  points(predictedN~I(2002:2014),type='l',col='blue')
  
  modelString<-paste(c(predictorSets$fall[['extreme']],
                       predictorSets$winter[['extreme']],
                       predictorSets$spring[['extreme']],
                       predictorSets$summer[['extremeInteraction']],
                       predictorSets$stock[['both']]),
                     collapse="+")
  model<-formula(paste("logPerCapita~",modelString))
  extremeFit<-lm(model,data=data[river==r])
  extremeFits[[r]]<-extremeFit
  predictedN<-exp(extremeFit$fitted.values)*data[river==r,fallAdults]
  points(predictedN~I(2002:2014),type='l',col='red')
  bestModels<-rbind(bestModels,data.table(fall="extreme",
                                          winter="extreme",
                                          spring="extreme",
                                          summer="extremeInteraction",
                                          stock="both",
                                          aic=AIC(extremeFit),
                                          rmse=rmse(extremeFit),
                                          river=r,
                                          rsq=summary(extremeFit)$r.squared))
  
  modelString<-paste(c(predictorSets$fall[['mean']],
                       predictorSets$winter[['mean']],
                       predictorSets$spring[['mean']],
                       predictorSets$summer[['meanInteraction']],
                       predictorSets$stock[['mean']]),
                     collapse="+")
  model<-formula(paste("logPerCapita~",modelString))
  meanFit<-lm(model,data=data[river==r])
  meanFits[[r]]<-meanFit
  predictedN<-exp(meanFit$fitted.values)*data[river==r,fallAdults]
  points(predictedN~I(2002:2014),type='l',col='orange')
  
  bestModels<-rbind(bestModels,data.table(fall="mean",
                                          winter="mean",
                                          spring="mean",
                                          summer="meanInteraction",
                                          stock="both",
                                          aic=AIC(meanFit),
                                          rmse=rmse(meanFit),
                                          river=r,
                                          rsq=summary(meanFit)$r.squared))
}
dev.off()
setkey(bestModels,river,aic)
