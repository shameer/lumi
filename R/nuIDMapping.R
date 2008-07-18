nuID2probeID <- function(nuID, lib="lumiHumanV1") {
	cat('This function is obsoleted!\nPlease use "nuID2IlluminaID" instead.\n')
	if (length(nuID) == 0) return(NULL)

	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (require(lib, character.only=TRUE)) {
		env <- get(paste(lib, 'PROBEID2NUID', sep = ""), mode = "environment")
		probe2nuID <- unlist(as.list(env))	
		allProbe <- names(probe2nuID)
		keepInd <- !is.na(probe2nuID)
		allProbe <- allProbe[keepInd]
		names(probe2nuID) <- allProbe
		probe <- lapply(nuID, function(x) allProbe[probe2nuID == x])
		names(probe) <- nuID
		return(probe)
	} else {
		cat(paste(lib, ' annotation library is required!\n', sep=''))
	}
}

nuID2targetID <- function(nuID, lib="lumiHumanV1") {
	cat('This function is obsoleted!\nPlease use "nuID2IlluminaID" instead.\n')
	if (length(nuID) == 0) return(NULL)
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (require(lib, character.only=TRUE)) {
		env <- get(paste(lib, 'TARGETID2NUID', sep = ""), mode = "environment")
		target2nuID <- unlist(as.list(env))	
		allTarget <- names(target2nuID)
		keepInd <- !is.na(target2nuID)
		allTarget <- allTarget[keepInd]
		names(target2nuID) <- allTarget
		target <- lapply(nuID, function(x) allTarget[target2nuID == x])
		names(target) <- nuID
		return(target)
	} else {
		cat(paste(lib, ' annotation library is required!\n', sep=''))
	}
}

probeID2nuID <- function(probeID, lib="lumiHumanV1") {
	cat('This function is obsoleted!\nPlease use "IlluminaID2nuID" instead.\n')
	if (length(probeID) == 0) return(NULL)
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (!require(annotate)) cat('Please install "annotate" library!\n')
	if (require(lib, character.only=TRUE)) {
		nuID <- unlist(lookUp(probeID, lib, 'PROBEID2NUID'))
		return(nuID)
	} else {
		cat(paste(lib, ' annotation library is required!\n', sep=''))
	}
}

targetID2nuID <- function(targetID, lib="lumiHumanV1") {
	cat('This function is obsoleted!\nPlease use "IlluminaID2nuID" instead.\n')
	if (length(targetID) == 0) return(NULL)
	if (!is.null(lib)) {
		if (length(grep('\\.db', lib)) > 0) {
			warning(paste(lib, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (!require(annotate)) cat('Please install "annotate" library!\n')
	if (require(lib, character.only=TRUE)) {
		nuID <- unlist(lookUp(targetID, lib, 'TARGETID2NUID'))
		return(nuID)
	} else {
		cat(paste(lib, ' annotation library is required!\n', sep=''))
	}
}
