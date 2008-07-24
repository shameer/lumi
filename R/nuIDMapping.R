nuID2probeID <- function(nuID, lib.mapping="lumiHumanIDMapping.db", ...) {

	if (length(nuID) == 0) return(NULL)

	if (!is.null(lib.mapping)) {
		if (length(grep('\\.db', lib.mapping)) > 0) {
			warning(paste(lib.mapping, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (require(lib.mapping, character.only=TRUE)) {
		if (length(grep('IDMapping.db', lib.mapping)) > 0) {
			probe <- nuID2IlluminaID(nuID, lib.mapping=lib.mapping, type='Probe', ...)
		} else {
			env <- get(paste(lib.mapping, 'PROBEID2NUID', sep = ""), mode = "environment")
			probe2nuID <- unlist(as.list(env))	
			allProbe <- names(probe2nuID)
			keepInd <- !is.na(probe2nuID)
			allProbe <- allProbe[keepInd]
			names(probe2nuID) <- allProbe
			probe <- lapply(nuID, function(x) allProbe[probe2nuID == x])
			names(probe) <- nuID			
		}
		return(probe)
	} else {
		cat(paste(lib.mapping, ' ID mapping library is required!\n', sep=''))
	}
}

nuID2targetID <- function(nuID, lib.mapping="lumiHumanIDMapping.db", ...) {

	if (length(nuID) == 0) return(NULL)
	if (!is.null(lib.mapping)) {
		if (length(grep('\\.db', lib.mapping)) > 0) {
			warning(paste(lib.mapping, 'does not include nuID conversion information!'))
			return(nuID)
		}
	}
	if (require(lib.mapping, character.only=TRUE)) {
		if (length(grep('IDMapping.db', lib.mapping)) > 0) {
			target <- nuID2IlluminaID(nuID, lib.mapping=lib.mapping, type='Gene', ...)
		} else {
			env <- get(paste(lib.mapping, 'TARGETID2NUID', sep = ""), mode = "environment")
			target2nuID <- unlist(as.list(env))	
			allTarget <- names(target2nuID)
			keepInd <- !is.na(target2nuID)
			allTarget <- allTarget[keepInd]
			names(target2nuID) <- allTarget
			target <- lapply(nuID, function(x) allTarget[target2nuID == x])
			names(target) <- nuID
		}
		return(target)
	} else {
		cat(paste(lib.mapping, ' ID mapping library is required!\n', sep=''))
	}
}

probeID2nuID <- function(probeID, lib.mapping="lumiHumanIDMapping.db", ...) {

	if (length(probeID) == 0) return(NULL)
	if (!is.null(lib.mapping)) {
		if (length(grep('\\.db', lib.mapping)) > 0) {
			warning(paste(lib.mapping, 'does not include nuID conversion information!'))
			return(probeID)
		}
	}
	if (!require(annotate)) cat('Please install "annotate" library!\n')
	if (require(lib.mapping, character.only=TRUE)) {
		if (length(grep('IDMapping.db', lib.mapping)) > 0) {
			nuID <- IlluminaID2nuID(probeID, lib.mapping=lib.mapping, ...)
		} else {
			nuID <- unlist(lookUp(probeID, lib.mapping, 'PROBEID2NUID'))
		}
		return(nuID)
	} else {
		cat(paste(lib.mapping, ' ID mapping library is required!\n', sep=''))
	}
}

targetID2nuID <- function(targetID, lib.mapping="lumiHumanIDMapping.db", ...) {

	if (length(targetID) == 0) return(NULL)
	if (!is.null(lib.mapping)) {
		if (length(grep('\\.db', lib.mapping)) > 0) {
			warning(paste(lib.mapping, 'does not include nuID conversion information!'))
			return(targetID)
		}
	}
	if (!require(annotate)) cat('Please install "annotate" library!\n')
	if (require(lib.mapping, character.only=TRUE)) {
		if (length(grep('IDMapping.db', lib.mapping)) > 0) {
			nuID <- IlluminaID2nuID(targetID, lib.mapping=lib.mapping, ...)
		} else {
			nuID <- unlist(lookUp(targetID, lib.mapping, 'TARGETID2NUID'))
		}
		return(nuID)
	} else {
		cat(paste(lib.mapping, ' ID mapping library is required!\n', sep=''))
	}
}
