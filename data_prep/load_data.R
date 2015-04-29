tribs<-data.table(dbGetQuery(link$conn,
       "SELECT * FROM tags_tribs_wb WHERE species = 'bkt' OR species = 'bnt' OR species = 'ats';"))
tribs[,line_number:=NULL]
trout<-data.table(dbGetQuery(link$conn,
        "SELECT * FROM tags_trout_wb WHERE species = 'bkt' OR species = 'bnt' OR species = 'ats';"))
setnames(trout,tolower(names(trout)))
trout<-rbind(trout,tribs)
rm(tribs)

trout<-trout[survey=="shock"]
trout[,measured_length:=as.numeric(measured_length)]
trout[tag=="0.9999",tag:=as.character(NA)]
trout[,date:=as.POSIXct(parse_date_time(date,orders='mdy'))]
trout[,measured_weight:=as.numeric(measured_weight)]
trout[measured_weight==0,measured_weight:=NA]
trout[measured_length==0,measured_length:=NA]

trout$index<-1:nrow(trout)
trout<-trout[!section %in% c(-1,0)]

assign('trout',trout,env=shared_data)
