root <- '~/trout_yoy'

start_year <- 2002
end_year <- 2014

do <- list(
	data_prep = c(
		'load_data.R',
    'create_no_catch.R',
    'create_hatch_year.R',
    'create_abundance_arrays_by_section.R',
    'create_abundance_arrays_by_river.R',
    'create_detection_array.R',
    'create_adult_array.R'
	),
  env_prep = c(
    'load_data.R',
    'add_early_trib_predictions.R',
    
    #need to run these only if there are changes to the input files
     #'importHighResTemp.R', 
     #'fillHighResTemps.R',
    'loadHighResEnv.R',#loads highResEnv from db if it isn't created in the above two scripts
    
    'predictWestBrookDischarge.R',
    'predictTribDischargeFromWestBrook.R',
    
    'fillHighResDepths.R',
    'convertDepthToDischarge.R',
    
    'createHighResTempCovariates.R',
    
    'covariate_inputs.R',
    'create_covariates.R',
    #'covariate_inputs_bnt.R',
    #'create_covariates_bnt.R',
    #'create_covariates_old.R'
    #,'create_loess_covariates.R' #not functional and abandoned
  ),
  save = 'saveData.R',
  
  analyze = c('runModelSpeciesRiversSeparate.R',
              'N.R',
              'NextremeStock.R',
              'stockRecruit.R',
              'envCovariatesProportionalChange.R',
              'compareRmseSpeciesRiversSeparate.R',
              'posteriorCorrelations.R'
              )
)

source(
	file=file.path(root,'shared_data.R'), 
	echo=TRUE, verbose=TRUE
)
for(stage in names(do)){
	for (script in do[[stage]]) {
		temp <- new.env(parent=shared_data)
		temp[['shared_data']] <- shared_data
		with(
			data=temp,
			expr= {
				s <- file.path(root,stage,script)
				cat(s,"\n")
				source(file=s, local=TRUE, echo=TRUE, verbose=TRUE)
			}
		)
		rm(envir=temp, list='shared_data')
		rm(temp)
	}	
}

rm(shared_data)

