env<-dcast(melt(covariates),Var1+Var2~Var3)
env<-data.table(env)
adult<-data.table(melt(A))
setnames(adult,c("year","river","species","adults"))
adult[,year:=year+1]
adult[,species:=as.integer(ifelse(species=='bkt',1,2))]
adult<-adult[year<2015]
setnames(env,c("Var1","Var2"),c("year","river"))
setkey(env,river,year)
setkey(popEst,river,year,species)
popEnv<-env[popEst]
setkey(popEst,river,year,species)
setkey(adult,river,year,species)
popEnv[,adults:=adult$adults]
popEnv[,riverFactor:=as.factor(river)]

envVariables<-c("adults","spring_floods","fall_discharge",
        "winter_floods","summer_low","summer_temp")

for(r in 1:4){
assign(paste0('fit',r),randomForest(I(mean/adults)~adults+
        spring_floods+fall_discharge+
        winter_floods+summer_low+summer_temp,
        data=popEnv[species==1&river==r],
        ntree=2000,
        mtry=4))

tiff.par(paste0("~/trout_yoy/results/figures/randomForestPartial",
                r,".tif"),
         mfrow=c(2,5),width=8,height=4)
for(i in 1:length(envVariables)){
  partialPlot(get(paste0('fit',r)),popEnv,envVariables[[i]],
              main=NA,xlab=envVariables[[i]])
}
#varImpPlot(get(paste0('fit',r)))
dev.off()
}
