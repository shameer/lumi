`nuID2RefSeqID` <-
function(nuID, lib.mapping, filterTh=c(Strength1=95, Uniqueness=95), returnAllInfo=FALSE) {
	if (missing(nuID) || missing(lib.mapping)) stop('Please provide all input parameters: nuID and lib.mapping!')
	# 
	mappingInfo <- getNuIDMappingInfo(nuID, lib.mapping=lib.mapping)
	if (!is.null(filterTh)) {
		if (!all(names(filterTh) %in% names(mappingInfo))) {
			cat('The provided names of filterTh does not match the field names of nuID_MappingInfo table!\n No filtering will be performed!\n')
		} else {
			filterNames <- names(filterTh)
			selInd <- TRUE
			for (i in seq(filterTh)) {
				selInd <- selInd & as.numeric(mappingInfo[, filterNames[i]]) > filterTh[i]
			}
			mappingInfo <- mappingInfo[selInd,]
		}
	}
	if (!returnAllInfo) mappingInfo <- mappingInfo[,'Accession']
	return(mappingInfo)	
}

