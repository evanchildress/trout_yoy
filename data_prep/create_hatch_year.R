data_root<-"~/process-data/data_store/processed_data"

yoy_bins<-readRDS(file.path(data_root,"yoy_bins.rds"))

hatch_year<-readRDS(file.path(data_root,"hatch_year.rds"))
setkey(hatch_year,tag)

trout[!is.na(tag),first_capture:=min(date),by=tag]
trout[is.na(tag),first_capture:=date]
setkey(trout,tag)

trout<-hatch_year[trout]

hatch.year<-function(Length,Sample,Species,River){
  if(is.na(Length)){return(as.numeric(NA))} else{
    upper<-yoy_bins[sample==Sample&species==Species&river==River,max_length]
    lower<-yoy_bins[sample==Sample&species==Species&river==River,min_length]
    hatch_year<-yoy_bins[sample==Sample&species==Species&river==River,hatch_year][intersect(which(Length<upper),which(Length>lower))]
  if(length(hatch_year)>0) {return(as.numeric(hatch_year))} else{
    return(as.numeric(NA))
  }
  }
}

trout[is.na(hatch_year),hatch_year:=hatch.year(Length=measured_length,Sample=sample_name,Species=species,River=river),by=index]
yoy<-trout[hatch_year==year(date)]
yoy[,sample_name:=as.numeric(sample_name)]

assign('yoy',yoy,env=shared_data)
assign('trout',trout,env=shared_data)
dbWriteTable(conn=link$conn, name='data_yoy',value=yoy, overwrite=TRUE, row.names=FALSE)


