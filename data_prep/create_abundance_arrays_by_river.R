nsections<-c(14,14,15,47)
trout<-trout[species %in% c('bkt','bnt')]
# sample_dates<-trout[,median(date),by=sample_name]
# setnames(sample_dates,"V1","median_sample_date")
# 
# season_breaks<-dbGetQuery(link$conn,"SELECT * FROM season_data")
# 
# season<-function(date){ 
#             length(which(strftime(date, format="%j") >= season_breaks$end_julian_day))+1
# }
# for(i in 1:nrow(sample_dates)){
# sample_dates[i,season:=season(median_sample_date)]
# }

sampleNames<-readRDS(file.path(data_root,"sampleNames.rds"))
sampleNames[,sample_name:=as.character(sample_name)]
sampleNames<-sampleNames[drainage=='west',list(sample_name,season)]
sampleNames[sample_name=="89",season:="PostSmolt"]



setkey(trout,sample_name)
setkey(sampleNames,sample_name)

trout<-sampleNames[trout]
trout<-trout[season!="Summer"]

trout<-trout[,season:=as.numeric(factor(trout$season,
                                 levels=c("PreSmolt","PostSmolt","Fall","PreWinter"),
                                 ordered=T))]
trout<-trout[,year:=year(date)]

abund.array<-function(data,season_set=3,age_class){
  if(age_class=="adult"){                                     
    data_subset<-data[cohort!=year & season==season_set]
  }
  
  if(age_class=="yoy"){                                     
    data_subset<-data[cohort==year & season==season_set]
  }
  
  abund_array<-acast(melt(data_subset[,length(tag),by=c("species","year","river")],
                          id=c("species","year","river")),year~river~species)
  
  abund_array[is.na(abund_array)]<-0
  
  return(abund_array)
}

yoy_abund_fall<-abund.array(trout,season_set=3,age_class='yoy')
adult_abund_fall<-abund.array(trout,season_set=3,age_class='adult')

y<-array(c(yoy_abund_fall,adult_abund_fall),dim=c(dim(yoy_abund_fall),2))
dimnames(y)<-c(dimnames(yoy_abund_fall),list(c("yoy","adult")))


propTagged<-trout[year==cohort & season==3,
                  length(unique(tag))/(sum(is.na(tag))+length(unique(tag))),
                  by=list(species,river,year)]

propTagged<-acast(melt(propTagged,id.vars=c("species","river","year")),
                  year~river~species)
propTagged[is.na(propTagged)]<-1

# save(yoy_abund_fall,adult_abund_fall,y,propTagged,
#      file=file.path(data_root,"abundanceArraysByRiver.rDATA"))

assign('propTagged',propTagged,env=shared_data)
assign('y',y,env=shared_data)
assign('trout',trout,env=shared_data)
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


