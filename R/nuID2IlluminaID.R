`nuID2IlluminaID` <-
function(nuID, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), idType=c('All', 'Probe', 'Gene', 'Accession', 'Search_key', 'Symbol'), chipVersion=NULL, ...) {
	
	## retrieve Illumina IDs from Id Mapping 
	retrieveIlluminaID <- function(idMapping, idType) {		
		mapName <- colnames(idMapping)
		illuminaID <- switch(idType,
			'All'=idMapping,
			'Probe'= idMapping[,grep('probe', mapName, ignore.case=TRUE)],
			'Gene'= idMapping[,grep('target|gene', mapName, ignore.case=TRUE)],
			'Accession'= idMapping[,grep('accession', mapName, ignore.case=TRUE)],
			'Search_key'= idMapping[,grep('search_key', mapName, ignore.case=TRUE)],
			'Symbol'= idMapping[,grep('symbol', mapName, ignore.case=TRUE)])
		return(illuminaID)
	}
		
	species <- match.arg(species)
	idType <- match.arg(idType)
	if ((is.null(chipVersion) || species == 'Unknown')) {
		chipInfo <- getChipInfo(nuID, lib.mapping=lib.mapping, species=species, idMapping=TRUE, ...)
		if (chipInfo$IDType[1] != 'nuID') {
			cat('The input ID is not nuID. So no ID converstion will be made!\n')
			return(nuID)
		}
		idMapping <- chipInfo$idMapping
		if (is(idMapping, 'list')) {
			illuminaID <- lapply(idMapping, function(x) retrieveIlluminaID(x, idType))
		} else {
			illuminaID <- retrieveIlluminaID(idMapping, idType)
		}
		return(illuminaID)
	} else {
		## directly get the table and return the mapping
		if (is.null(lib.mapping)) {
			lib.mapping <- switch(species,
				'Rat'='lumiRatIDMapping.db',
				'Human'='lumiHumanIDMapping.db',
				'Mouse'='lumiMouseIDMapping.db')			
		}
		if(!require(lib.mapping, character=TRUE)) stop(paste(lib.mapping, 'is required!'))
		dbconn <- sub("\\.db", "_dbconn", lib.mapping)
		conn <- do.call(dbconn, list())

		allTableNames <- dbListTables(conn)	
		if (!(chipVersion %in% allTableNames)) {
			cat('The provided chipVersion cannot be found in the database!\n We will search over entire database for ID mapping!\n')
			return(nuID2IlluminaID(nuID, lib.mapping, chipVersion=NULL, ...))
		}
	}
}

