`IlluminaID2nuID` <-
function(IlluminaID, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), chipVersion=NULL, ...) {
	species <- match.arg(species)

	chipInfo <- getChipInfo(IlluminaID, lib.mapping=lib.mapping, species=species, idMapping=TRUE, chipVersion=chipVersion, ...)
	
	if (is.null(chipInfo$IDType)) {
		nuID <- rep(NA, length(IlluminaID))
		names(nuID) <- IlluminaID
		return(nuID)
	} else if (!is.na(chipInfo$IDType)) {
		if (chipInfo$IDType[1] == 'nuID') {
			cat('The input ID is nuID. So no ID converstion will be made!\n')
			return(IlluminaID)
		}
	}		
	idMapping <- chipInfo$idMapping
	return(idMapping)
}

