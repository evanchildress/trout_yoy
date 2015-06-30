if(!exists("highResEnv")){
  highResEnv<-data.table(dbGetQuery(link$conn,
        "SELECT * FROM data_high_res_env_filled;"))

}

summerT<-highResEnv[month(dateTime) %in% 6:9,
                       length(which(temp>=20))/length(temp),
                       by=list(river,year(dateTime))]

summerT<-acast(melt(summerT,id.vars=c("river","year")),year~river)

meanSummerT<-highResEnv[month(dateTime) %in% 6:8,
                       mean(temp),
                       by=list(river,year(dateTime))]

meanSummerT<-acast(melt(meanSummerT,id.vars=c("river","year")),year~river)

assign("summerT",summerT,env=shared_data)
assign("meanSummerT",meanSummerT,env=shared_data)