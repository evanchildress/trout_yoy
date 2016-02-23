if(!exists("highResEnv")){
#   highResEnv<-data.table(dbGetQuery(link$conn,
#                                     "SELECT * FROM data_high_res_env_filled;"))
  highResEnv<-tbl(conn,"data_high_res_env_filled") %>%
              collect() %>%
              data.table()
  
}
assign("highResEnv",highResEnv,env=shared_data)