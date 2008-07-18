`getNuIDMappingInfo` <-
function(nuID, lib.mapping) {
	if (missing(nuID) || missing(lib.mapping)) stop('Please provide all input parameters: nuID and lib.mapping!')
	if(!require(lib.mapping, character=TRUE)) stop(paste(lib.mapping, 'is required!'))
	dbconn <- sub("\\.db", "_dbconn", lib.mapping)
	conn <- do.call(dbconn, list())
    
	allTableNames <- dbListTables(conn)	
	fieldName <- dbListFields(conn, 'nuID_MappingInfo')
	nuIDMappingInfo <- dbReadTable(conn, 'nuID_MappingInfo')
	rownames(nuIDMappingInfo) <- nuIDMappingInfo[,'nuID']
	nuIDMappingInfo <- nuIDMappingInfo[,-1]
	mappingInfo <- nuIDMappingInfo[nuID,]
	return(mappingInfo)
}

