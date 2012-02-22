`getNuIDMappingInfo` <-
function(nuID=NULL, lib.mapping) {

	if (missing(lib.mapping)) stop('Please specify lib.mapping library!')
	if(!require(lib.mapping, character.only=TRUE)) stop(paste(lib.mapping, 'is required!'))
	# dbconn <- sub("\\.db", "_dbconn", lib.mapping)
	dbconn <- paste(lib.mapping, "_dbconn", sep='')
	conn <- do.call(dbconn, list())
    
	allTableNames <- dbListTables(conn)	
	fieldName <- dbListFields(conn, 'nuID_MappingInfo')
	nuIDMappingInfo <- dbReadTable(conn, 'nuID_MappingInfo')
	allNuID <- nuIDMappingInfo[,'nuID']
	rownames(nuIDMappingInfo) <- allNuID
	nuIDMappingInfo <- nuIDMappingInfo[,-1]
	if (!is.null(nuID)) {
		# check unID
		if (!all(sapply(nuID, is.nuID))) stop('Some inputted nuIDs are not real nuIDs!\n')
		mappingInfo <- matrix(rep(NA, ncol(nuIDMappingInfo)*length(nuID)), ncol=ncol(nuIDMappingInfo))
		rownames(mappingInfo) <- nuID
		colnames(mappingInfo) <- colnames(nuIDMappingInfo)
		selNuID <- nuID[nuID %in% allNuID]
		if (length(selNuID) == 0) {
			warning('No matches were found!\n')
		} else {
			if (length(selNuID) < length(nuID)) warning('Some input IDs can not be matched!\n')
			mappingInfo[selNuID, ] <- as.matrix(nuIDMappingInfo[selNuID,])
		}
	} else {
		mappingInfo <- nuIDMappingInfo
	}
	return(mappingInfo)
}

