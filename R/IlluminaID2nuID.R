`IlluminaID2nuID` <-
function(IlluminaID, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), chipVersion=NULL, ...) {
	species <- match.arg(species)
	if ((is.null(chipVersion) || species == 'Unknown')) {
		chipInfo <- getChipInfo(IlluminaID, lib.mapping=lib.mapping, species=species, idMapping=TRUE, ...)
		if (is.null(chipInfo$IDType)) {
			nuID <- rep(NA, length(IlluminaID))
			names(nuID) <- IlluminaID
			return(nuID)
		} else if (chipInfo$IDType[1] == 'nuID') {
			cat('The input ID is nuID. So no ID converstion will be made!\n')
			return(IlluminaID)
		}
		idMapping <- chipInfo$idMapping
		return(idMapping)
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
			return(IlluminaID2nuID(IlluminaID, lib.mapping, chipVersion=NULL, ...))
		}
	}
}

