nsections<-c(14,14,15,47)

sample_dates<-trout[,median(date),by=sample_name]
setnames(sample_dates,"V1","median_sample_date")

season_breaks<-dbGetQuery(link$conn,"SELECT * FROM season_data")

season<-function(date){ 
            length(which(strftime(date, format="%j") >= season_breaks$end_julian_day))+1
}
for(i in 1:nrow(sample_dates)){
sample_dates[i,season:=season(median_sample_date)]
}



abund.array<-function(data,season_set=3,age_class){
  if(age_class=="adult"){                                     
    data_subset<-data[cohort!=year(date) & sample_name %in% sample_dates[season==season_set,sample_name]]
  }
  
  if(age_class=="yoy"){                                     
    data_subset<-data[cohort==year(date) & sample_name %in% sample_dates[season==season_set,sample_name]]
  }
  
  abund_array<-acast(melt(data_subset[,length(tag),by=c("species","sample_name","section","river")],
                          id=c("species","sample_name","section","river")),as.numeric(section)~sample_name~river~species)
  
  abund_array[is.na(abund_array)]<-0
  rivers<-unique(data$river)[order(unique(data$river))]
  
  nsections<-c(14,14,15,47)
  
  for(r in 1:3){
    abund_array[(nsections[r]+1):dim(abund_array)[1],,rivers[r],]<-NA
  }
  
  abund_array[,c("30","36"),rivers[1:3],]<-NA
  
  return(abund_array)
}

yoy_abund_fall<-abund.array(trout,season_set=3,age_class='yoy')
adult_abund_fall<-abund.array(trout,season_set=3,age_class='adult')

y<-array(c(yoy_abund_fall,adult_abund_fall),dim=c(dim(yoy_abund_fall),2))
dimnames(y)<-c(dimnames(yoy_abund_fall),list(c("yoy","adult")))

save(yoy_abund_fall,adult_abund_fall,y,file=file.path(data_root,"abundance_arrays.rDATA"))

#yoy_fall_array[,,c("wb jimmy","wb mitchell","wb obear"),,2]<-NA
  
#This section uses the no catch indicators to replace NAs with 0s, which is overriden by the above
#bad<-NULL
#bad.array<-NULL
#for(i in 1:nrow(no_catch[sample_name %in% sample_dates[season==3,sample_name]])){
#  with(no_catch[sample_name %in% sample_dates[season==3,sample_name]],
#       expr={
#       pass.i<<-pass[i]
#       species.i<<-species[i]
#       river.i<<-river[i]
#       sample_name.i<<-sample_name[i]
#       section.i<<-section[i]
#       }
#       )
  
#  if(!is.na(yoy_fall_array[section.i,sample_name.i,river.i,species.i,pass.i])){
#    bad<-c(bad,i)
#    bad.array<-c(bad.array,yoy_fall_array[section.i,sample_name.i,river.i,species.i,pass.i])
#    cat("warning:",paste(section.i,sample_name.i,river.i,species.i,pass.i,sep=","),"is",yoy_fall_array[section.i,sample_name.i,river.i,species.i,pass.i],"not NA","\n")
#    } else{
#      yoy_fall_array[section.i,sample_name.i,river.i,species.i,pass.i]<-0
#    }
#}


