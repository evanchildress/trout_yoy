require(readxl)

originalDir<-"~/process-data/data_store/original_data"
solinstFile<-file.path(originalDir,"West Drainage Solinst logger data.xlsx")
tempLogFile<-file.path(originalDir,"finished west brook temperature.xlsx")

sheetNames<-excel_sheets(solinstFile)

objectNames<-NULL
colTypes<-c("text","text","numeric","numeric","numeric","numeric","text","text")
for(n in sheetNames){
  name<-paste0(gsub(" ","",n),"Solinst")
  assign(name,
         data.table(read_excel(solinstFile,n,
                               col_types=colTypes)))
  setnames(get(name),
           c("river","section","date","time","depth","temp",
             "source","comment"))
  assign(name,get(name)[,list(river,section,date,time,temp,depth,source)])
  get(name)[,dateTime:=as.POSIXct(date*24*60*60+time*24*60*60,
                                  origin="1899-12-30",tz="GMT")]
  
  objectNames<-c(objectNames,as.name(name))

}

sheetNamesLog<-c("WB Jimmy","WB Mitchell","WB Obear","section 45")

for(n in 1:4){
  name<-paste0(gsub(" ","",sheetNames[n]),"Log")
  if(grepl("jimmy",name)|grepl("obear",name)){
    colTypes<-c("text","numeric","numeric","numeric","text","text")
  } else if(grepl("westbrook",name)){
    colTypes<-c("text","text","numeric","numeric","numeric","text","text")
  } else {
    colTypes<-c("text","numeric","numeric","numeric","text","text","text")
  }
  assign(name,
         data.table(read_excel(tempLogFile,sheetNamesLog[n],
                               col_types=colTypes)))
  if(grepl("westbrook",name)){
      setnames(get(name),
           c("river","section","date","time","temp",
             "source","notes"))
  } else if(grepl("mitch",name)){
  setnames(get(name),
           c("river","date","time","temp","section",
             "source","notes"))
  } else {
      setnames(get(name),
           c("river","date","time","temp","section",
             "source"))
  }
  
  get(name)[,depth:=NA]
  assign(name,get(name)[,list(river,section,date,time,temp,depth,source)])
    get(name)[,dateTime:=as.POSIXct(date*24*60*60+time*24*60*60,
                                  origin="1899-12-30",tz="GMT")]

  objectNames<-c(objectNames,as.name(name))
}

data<-do.call(rbind,as.list(objectNames))
data[,river:=tolower(river)]
setkey(data,river,dateTime)

dbWriteTable(link$conn, 'data_high_resolution_env', data, row.names=FALSE,
						 overwrite=TRUE, append=FALSE)

assign("highResEnv",data,env=shared_data)