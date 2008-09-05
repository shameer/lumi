`nuID2IlluminaID` <-
function(nuID, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), idType=c('All', 'Probe', 'Gene', 'Accession', 'Search_key', 'Symbol'), chipVersion=NULL, ...) {
	
	## retrieve Illumina IDs from Id Mapping 
	retrieveIlluminaID <- function(idMapping, idType) {		
		mapName <- colnames(idMapping)
		illuminaID <- switch(idType,
			'All'=idMapping,
			'Probe'= {
				ind <- grep("probe_id", mapName, ignore.case=TRUE)
				if (length(ind) == 0) ind <- grep("probe", mapName, ignore.case=TRUE)
				idMapping[,ind]
				},
			'Gene'= idMapping[,grep('target|gene', mapName, ignore.case=TRUE)],
			'Accession'= idMapping[,grep('accession', mapName, ignore.case=TRUE)],
			'Search_key'= idMapping[,grep('search_key', mapName, ignore.case=TRUE)],
			'Symbol'= idMapping[,grep('symbol', mapName, ignore.case=TRUE)])
		names(illuminaID) <- rownames(idMapping)
		return(illuminaID)
	}
		
	species <- match.arg(species)
	idType <- match.arg(idType)
	chipInfo <- getChipInfo(nuID, lib.mapping=lib.mapping, species=species, idMapping=TRUE, chipVersion=chipVersion, ...)
	idMapping <- chipInfo$idMapping
	if (length(nuID) > 0) {
		if (chipInfo$IDType[1] != 'nuID') {
			cat('The input ID is not nuID. So no ID converstion will be made!\n')
			return(nuID)
		}
	} 
	if (is(idMapping, 'list')) {
		illuminaID <- lapply(idMapping, function(x) retrieveIlluminaID(x, idType))
	} else {
		illuminaID <- retrieveIlluminaID(idMapping, idType)
	}
	return(illuminaID)		
}

