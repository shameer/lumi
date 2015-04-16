`nuID2RefSeqID` <-
function(nuID=NULL, lib.mapping, filterTh=c(Strength1=95, Uniqueness=95), returnAllInfo=FALSE) {
	if (missing(lib.mapping)) stop('Please specify lib.mapping library!')
	# 
	mappingInfo <- getNuIDMappingInfo(nuID, lib.mapping=lib.mapping)

	if (!('Refseq_old' %in% colnames(mappingInfo)) && !is.null(filterTh)) {
		if (!all(names(filterTh) %in% names(mappingInfo))) {
			cat('The provided names of filterTh does not match the field names of nuID_MappingInfo table.\n No filtering will be performed.\n')
		} else {
			filterNames <- names(filterTh)
			selInd <- TRUE
			for (i in seq(filterTh)) {
				selInd <- selInd & as.numeric(mappingInfo[, filterNames[i]]) > filterTh[i]
			}
			mappingInfo <- mappingInfo[selInd,]
		}
	}
	if (!returnAllInfo) {
		nuID <- rownames(mappingInfo)
		if ('Accession' %in% colnames(mappingInfo)) {
			mappingInfo <- mappingInfo[,'Accession']			
		} else if ('Refseq' %in% colnames(mappingInfo)) {
			mappingInfo <- mappingInfo[,'Refseq']			
		}
		names(mappingInfo) <- nuID
	}
	return(mappingInfo)	
}

