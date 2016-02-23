#create vector of rows that are just place holders indicating that no fish were caught
no_catch_rows<-with(trout, which((grepl('no',comments) & (grepl('pass',comments) | grepl('p2',comments)) | 
                      grepl('no fish',comments) | grepl('no bkt',comments) | grepl('dry',comments))))

#create a table with only those rows
no_catch<-trout[no_catch_rows]



#remove rows that serve as no catch indicators and update the shared_data version
trout<-trout[-no_catch_rows]

#make objects available in shared_data environment
assign('trout',trout,env=shared_data)
assign('no_catch',no_catch,env=shared_data)
