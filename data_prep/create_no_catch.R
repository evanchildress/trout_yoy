#create table of indicator rows
no_catch<-trout[grepl('no',comments)&grepl('pass',comments)]



#remove rows that serve as no catch indicators and update the shared_data version
trout<-trout[!(grepl('no',comments)&grepl('pass',comments))]
assign('trout',trout,env=shared_data)
assign('no_catch',no_catch,env=shared_data)
