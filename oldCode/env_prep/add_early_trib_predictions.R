
bla<-dcast.data.table(melt(ed[,list(river,date,discharge)],id=c("river","date")),date~river)
setnames(bla,c("date","jimmy","mitchell","obear","westbrook"))
bla[,c("jimmy.ln","mitchell.ln","obear.ln","westbrook.ln"):=
      list(log(jimmy),log(mitchell),log(obear),log(westbrook))]

jimmy.lm<-lm(jimmy.ln~westbrook.ln,data=bla)
mitchell.lm<-lm(mitchell.ln~westbrook.ln,data=bla)
obear.lm<-lm(obear.ln~westbrook.ln,data=bla)

bla[,jimmy:=exp(predict(jimmy.lm,data.frame(westbrook.ln)))]
bla[,mitchell:=exp(predict(mitchell.lm,data.frame(westbrook.ln)))]
bla[,obear:=exp(predict(obear.lm,data.frame(westbrook.ln)))]
bla<-bla[,list(date,jimmy,mitchell,obear,westbrook)]
setnames(bla,c("date","wb jimmy","wb mitchell","wb obear","west brook"))
bla<-data.table(melt(bla,id=c("date")))
setnames(bla,c("date","river","discharge"))
ed[,discharge:=NULL]
setkey(bla,river,date)
setkey(ed,river,date)
ed<-ed[bla]


bla<-dcast.data.table(melt(ed[,list(river,date,temperature)],id=c("river","date")),date~river)
setnames(bla,c("date","jimmy","mitchell","obear","westbrook"))

jimmy.lm<-lm(jimmy~westbrook,data=bla)
mitchell.lm<-lm(mitchell~westbrook,data=bla)
obear.lm<-lm(obear~westbrook,data=bla)

bla[is.na(jimmy),jimmy:=predict(jimmy.lm,data.frame(westbrook))]
bla[is.na(mitchell),mitchell:=predict(mitchell.lm,data.frame(westbrook))]
bla[is.na(obear),obear:=predict(obear.lm,data.frame(westbrook))]
setnames(bla,c("date","wb jimmy","wb mitchell","wb obear","west brook"))
bla<-data.table(melt(bla,id=c("date")))
setnames(bla,c("date","river","temperature"))
ed[,temperature:=NULL]
setkey(bla,river,date)
setkey(ed,river,date)
ed<-ed[bla]

for(i in ed[river=="west brook"&is.na(temperature),date]){
  jimmyTemp<-ed[date==i&river=="wb jimmy",temperature]
  if(length(jimmyTemp)>0){
    ed[date==i&river=="west brook",temp:=jimmyTemp]
  }}

ed<-ed[,list(date,river,discharge,temperature)]
assign('ed',ed,env=shared_data)

## Save cleaned daily data.
dbWriteTable(conn=link$conn, name='data_environmental', value=ed,
             row.names=FALSE, overwrite=TRUE, append=FALSE)
