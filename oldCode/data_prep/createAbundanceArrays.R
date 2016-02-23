abundanceArray<-trout[!is.na(stage),.N,by=list(species,year,season,river,stage)] %>%
                melt(id.vars=c("species","year","season","river","stage")) %>%
                acast(year~season~river~stage~species)

propTagged<-trout[stage==1 & season==3,
                  length(unique(tag))/.N,
                  by=list(species,river,year)] %>%
            melt(id.vars=c("species","river","year")) %>%
            acast(year~river~species)
propTagged[is.na(propTagged)]<-1

assign('propTagged',propTagged,env=shared_data)
assign('y',y,env=shared_data)
assign('trout',trout,env=shared_data)


