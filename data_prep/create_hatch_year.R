data_root<-"~/process-data/data_store/processed_data"

yoy_bins<-readRDS(file.path(data_root,"yoy_bins.rds"))

#hatch_year<-readRDS(file.path(data_root,"hatch_year.rds"))
#setkey(hatch_year,tag)

trout[!is.na(tag),first_capture:=min(date),by=tag]
trout[is.na(tag),first_capture:=date]
setkey(trout,tag)
trout[,cohort := as.numeric(cohort)]
trout[!is.na(tag),cohort := ifelse(any(!is.na(cohort)),
                                   min(cohort[which(!is.na(cohort))]),as.numeric(NA)),
      by=tag]

#trout<-hatch_year[trout]

hatch.year.no.tag<-function(Length,Sample,Species,River){
  if(all(is.na(Length))){return(as.numeric(NA))} else{
    
    #A few samples/species/rivers don't have cohorts assigned, and these adjustments use the best substitute to fill the gaps
    if(Sample == 41.8){Sample <- 41}
    if(Species== 'bnt' & Sample %in% c(40:54,56:60) & River == 'wb jimmy'){River <- 'west brook'} #No assignments have been made for jimmy pre sample 60, but river specific breaks could be determined
    if(Species== 'bnt' & Sample ==55 & River == 'wb jimmy'){River <- 'wb mitchell'} #No assignments have been made for jimmy pre sample 60, but river specific breaks could be determined
    if(Species== 'bnt' & Sample == 40 & River == 'wb mitchell'){River <- 'west brook'} #no bins assinged for this sample in mitchell
    if(Species== 'bnt' & Sample == 47 & River == 'wb mitchell' & Length>200){River <- 'west brook'} #these guys are bigger than the largest bin that was assigned
    if(Species== 'bnt' & Sample == 82 & River == 'wb jimmy' & Length>210){River <- 'west brook'} #bigger than assigned bins
    if(Species== 'bnt' & Sample == 91 & River == 'wb mitchell' & Length == 90){return(2014)}
    
    upper<-yoy_bins[sample==Sample&species==Species&river==River,max_length]
    lower<-yoy_bins[sample==Sample&species==Species&river==River,min_length]
    hatch_year<-yoy_bins[sample==Sample&species==Species&river==River,hatch_year][
      intersect(which(Length<=upper),which(Length>=lower))]
  if(length(hatch_year)>0) {return(as.numeric(hatch_year))} else{
    return(as.numeric(NA))
  }
  }
}

hatch.year.tag<-function(Length,Sample,Species,River) {
  if(all(is.na(Length))){return(as.numeric(NA))} else{
  firstObs<-which(Length==min(Length,na.rm=T))[1]
  return(hatch.year.no.tag(Length[firstObs],Sample[firstObs],Species[firstObs],River[firstObs]))
  }}

trout[,cohort_estimated:= is.na(cohort)]
trout[is.na(cohort) & is.na(tag),cohort:=
        hatch.year.no.tag(Length=measured_length,
        Sample=sample_name,Species=species,River=river),
      by=index]

trout[is.na(cohort) & !is.na(tag),cohort:=
        hatch.year.tag(Length=measured_length,
        Sample=sample_name,Species=species,River=river),
      by=tag]

trout[is.na(cohort) & grepl('too small to tag',comments),cohort:= year(date)]

if(nrow(trout[is.na(cohort) & !is.na(measured_weight)])>0){
  cohort.from.weight<-function(Weight,Sample,Species,River){
    weight_lm<-lm(log(measured_length)~log(measured_weight),data=trout[species==Species])
    predicted_length<-exp(predict(weight_lm,data.frame(measured_weight=Weight)))  
    cohort<-hatch.year.no.tag(predicted_length,Sample,Species,River)
    return(cohort)
  }
  
  trout[is.na(cohort) & !is.na(measured_weight),cohort:=
          cohort.from.weight(Weight=measured_weight,
          Sample=sample_name,Species=species,River=river),
        by=index]
}

#################################################################################
##This bit is useful for seeing whether the fish with unassigned cohort are actually fish.
##As of 24-Mar-2015, none were real fish, just no catch indicators that didn't have a comment,
##this is indicated by no other fish being caught and the one line having no data
# for(i in 1:nrow(trout[is.na(cohort)])){
#   a<-trout[is.na(cohort)][i]
#   a<-trout[river == a$river &
#               sample_name == a$sample_name &
#               section == a$section &
#               species ==a$species &
#               pass == a$pass &
#               braid == a$braid]
#   if(nrow(a)>1){print(a[])
#   print("***********************************************************")}
# }
#################################################################################

##This removes all of the fish that weren't assigned a cohort 
##because most of these a supposed to be indicators that no fish were caught,
##but it produces a warning if there are any that seem like they might be real fish
for(i in 1:nrow(trout[is.na(cohort)])){
  a<-trout[is.na(cohort)][i]
  ab<-trout[river == a$river &
              sample_name == a$sample_name &
              section == a$section &
              species ==a$species &
              pass == a$pass &
              braid == a$braid]
  if(nrow(ab)>1){warning(paste0("id ",a$id," removed b/c cohort couldn't be assigned, but other fish were caught for that species/braid/pass/section/sample/river"))}
}
trout<-trout[!is.na(cohort)]


# yoy<-trout[cohort==year(date)]
# yoy[,sample_name:=as.numeric(sample_name)]



# dbWriteTable(conn=link$conn, name='data_yoy',value=yoy, overwrite=TRUE, row.names=FALSE)
dbWriteTable(conn=link$conn, name='data_trout',value=trout,overwrite=TRUE,row.names=FALSE)

# yoy<-yoy[area %in% c("inside","trib")]
trout<-trout[area %in% c("inside","trib")]
# assign('yoy',yoy,env=shared_data)
assign('trout',trout,env=shared_data)
 
