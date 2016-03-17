meanX<-rnorm(100,100,10)

x<-array(NA,dim=c(100,100))
for(t in 1:100){
  x[,t]<-rnorm(100,meanX[t],10)
}

y<-apply(x,2,mean)*-0.5+rnorm(100,0,3) #mean
#x70<-apply(x,2,function(d){sum(d>=quantile(x,0.1))})
#y<-x70*-0.5+rnorm(100,0,2)

thresholds<-seq(0.01,0.99,0.01)
results<-data.table(threshold=thresholds,rsq=as.numeric(NA),aic=as.numeric(NA))
for(thresh in thresholds){
  threshold<-quantile(x,thresh)
  xThresh<-apply(x,2,function(d){sum(d>threshold)})
  a<-lm(y~xThresh)
  results[threshold==thresh,":="(rsq=summary(a)$r.squared,
                                aic=AICc(a))]
}

a<-lm(y~apply(x,2,mean))
