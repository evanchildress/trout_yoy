load("~/process-data/data_store/processed_data/abundance_arrays.rDATA")

plotBySection<-function(river,species,stage){
  data<-y[,,river,species,stage]
  data<-data[!is.na(data[,1]),]
  cols<-colorRampPalette(c('yellow','black'))(nrow(data))
  
  plot(NA,xlim=c(2000,2014),ylim=c(0,max(data)))
  for(r in 1:nrow(data)){
    points(data[r,]~I(2000:2014),type='l',col=cols[r])
  }
}

plotByYear<-function(river,species,stage,standardize=F){
  data<-y[,,river,species,stage]
  data<-data[!is.na(data[,1]),]
  cols<-colorRampPalette(c('green','black'))(ncol(data))
  
  if(standardize){
    data<-(data-colMeans(data))/apply(data,2,sd)
  }
  
  plot(NA,xlim=c(1,nrow(data)),ylim=c(-10,max(data)))
  
  for(y in 1:ncol(data)){
    points(data[,y]~I(1:nrow(data)),type='l',col=cols[y],lwd=2)
  }
}

variogram<-function(river,species,stage,year){
  data<-y[,,river,species,stage]
  data<-data[!is.na(data[,1]),]
  data<-data[,year]
  
  dists<-dist(cbind(1:length(data),rep(1,length(data))))
  variog(coords=cbind(1:length(data),rep(1,length(data))),data=data,breaks=seq(0,47,1))
  
}
for(d in 1:15){
  bla<-variogram('west brook','bkt','yoy',d)
  #plot(NA,xlim=c(0,50),ylim=c(0,10))
  plot(bla,main=d)
}


plot(apply(y[,,'west brook','bkt','adult'],c(1),function(x){sd(x)/mean(x)}))
points(apply(y[,,'west brook','bnt','adult'],c(1),function(x){sd(x)/mean(x)}),pch=19)



plot(apply(y[,,'west brook','bnt','adult'],c(1),mean)~
       apply(y[,,'west brook','bnt','yoy'],c(1),mean),ylab='')
bla<-apply(y[,,'west brook','bkt','adult'],c(1),function(x){sd(x)/mean(x)})-
       apply(y[,,'west brook','bnt','adult'],c(1),function(x){sd(x)/mean(x)})
