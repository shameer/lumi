`getChipInfo` <-
function(x, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), chipVersion=NULL, idMapping=FALSE, returnAllMatches=FALSE, verbose=TRUE) {
	
	## Function to make sure the output mapping matches the inputID.bak
	matchInputID <- function(mapping, inputID.bak) {
		## We assume mapping IDs are a subset of the names of inputID.bak
		nc <- ncol(mapping)
		if (is.null(nc)) {
			len <- length(mapping)
			if (len == length(inputID.bak)) return(mapping)
			mapping.new <- rep(NA, length(inputID.bak))
			names(mapping.new) <- inputID.bak
			mappingID <- names(mapping)
			mapping.new[mappingID] <- mapping
			type <- 'vector'
		} else {
			nr <- nrow(mapping)
			if (nr == length(inputID.bak)) return(mapping)
			mapping.new <- matrix(NA, nrow=length(inputID.bak), ncol=nc)
			rownames(mapping.new) <- inputID.bak
			colnames(mapping.new) <- colnames(mapping)
			mappingID <- rownames(mapping)
			mapping.new[mappingID,] <- mapping[mappingID,]
			type <- 'matrix'
		}
		selID <- inputID.bak[inputID.bak %in% mappingID]
		dupInd <- which(duplicated(selID))
		if (length(dupInd) > 0) {
			dupID <- unique(selID[dupInd])
			for (dupId.i in dupID) {
				ind.i <- which(inputID.bak == dupId.i)
				if (type == 'matrix') {
					mapping.new[ind.i, ] <- mapping[dupId.i, ]
				} else {
					mapping.new[ind.i] <- mapping[dupId.i]					
				}
			}
		}
		return(mapping.new)
	}

	
	if (is(x, 'ExpressionSet')) {
		inputID <- featureNames(x)
	} else if (is(x, 'matrix') || is(x, 'data.frame')) {
		inputID <- rownames(x)		
	} else {
		inputID <- as.character(x)
	}
	inputID.bak <- inputID
	inputID <- unique(inputID)
	species <- match.arg(species)
	if (species == 'Unknown' && is.null(lib.mapping)) {
		human.match <- getChipInfo(inputID, species='Human', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		mouse.match <- getChipInfo(inputID, species='Mouse', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		rat.match <- getChipInfo(inputID, species='Rat', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		bestChip <- which.max(c(length(human.match$matchedProbeNumber), length(mouse.match$matchedProbeNumber), length(rat.match$matchedProbeNumber)))
		best.match <- switch(bestChip, human.match, mouse.match, rat.match)
		return(best.match)
	} else {
		if (is.null(lib.mapping)) {
			lib.mapping <- switch(species,
				'Rat'='lumiRatIDMapping.db',
				'Human'='lumiHumanIDMapping.db',
				'Mouse'='lumiMouseIDMapping.db')			
		}		
		if(!require(lib.mapping, character=TRUE)) stop(paste(lib.mapping, 'is required!'))
		dbconn <- sub("\\.db", "_dbconn", lib.mapping)
		conn <- do.call(dbconn, list())
		
		metaInfo <- dbReadTable(conn, 'metadata')
		species <- metaInfo$value[metaInfo[,'name'] == "SPECIES"]
		
		allTableNames <- dbListTables(conn)
		allTableNames <- allTableNames[!(allTableNames %in% c('nuID_MappingInfo', 'metadata'))]
		if (!is.null(chipVersion)) {
			if (chipVersion %in% allTableNames) {
				allTableNames <- chipVersion
			} else {
				warning('The provided chipVersion does not exist in the library!\n')
			}
		}
		lenID <- length(inputID)
		if (lenID == 0 && !is.null(chipVersion)) {
			## return the entire table based on chipVersion
			mapping <- dbReadTable(conn, chipVersion)
			returnList <- list(chipVersion=chipVersion, species=species, IDType=NA, chipProbeNumber=nrow(mapping), 
				inputProbeNumber=0, matchedProbeNumber=NA)
		} else {
			matchLen <- NULL
			matchField <- NULL
			tableLen <- NULL
			cut0.probeId <- FALSE  # determine whether to cut the 0s in front of the ProbeId. E.g. 0004760445 --> 4760445
			for (i in seq(allTableNames)) {
				tableName.i <- allTableNames[i]
				table.i <- dbReadTable(conn, tableName.i)
				fieldNames.i <- names(table.i)
				## Not include the Symbol column during matching the IDs.
				fieldNames.i <- fieldNames.i[!(fieldNames.i %in% c('Symbol'))]
				len.i <- NULL
				for (fieldNames.ij in fieldNames.i) {
					field.ij <- as.character(table.i[,fieldNames.ij])
					len.ij <- length(which(inputID %in% field.ij))
					if (fieldNames.ij == 'ProbeId') {
						len.ij.2 <- length(which(inputID %in% as.character(as.integer(field.ij))))
						if (len.ij.2 > len.ij) {
							cut0.probeId <- TRUE
							len.ij <- len.ij.2
						}
					}
					len.i <- c(len.i, len.ij)
				}
				max.ind.i <- which.max(len.i)
				matchLen <- c(matchLen, len.i[max.ind.i])
				matchField <- c(matchField, fieldNames.i[max.ind.i])
				tableLen <- c(tableLen, nrow(table.i))
			}
			bestMatchLen <- max(matchLen)
			if (bestMatchLen == 0) {
				if (verbose) warning('No matches were found!\n')
				return(list(chipVersion=NULL, species=species, IDType=NULL, chipProbeNumber=NULL, 
					matchedProbeNumber=bestMatchLen), idMapping=NULL)
			} else if (bestMatchLen < lenID && verbose) {
				warning('Some input IDs can not be matched!\n')
			}			

			if (!returnAllMatches) {
				bestInd <- which(matchLen == bestMatchLen)
				if (length(bestInd) > 1) {
					bestInd <- bestInd[tableLen[bestInd] == min(tableLen[bestInd])]
				}
				matchLen.match <- bestMatchLen
				tableName.match <- allTableNames[bestInd]
				fieldName.match <- matchField[bestInd]
				bestTable <- dbReadTable(conn, tableName.match[1])
				probeNumber <- nrow(bestTable)			
				if (idMapping) {
					## remove duplicated IDs becasue rownames should be unique
					dupInd <- duplicated(bestTable[,fieldName.match[1]])
					bestTable <- bestTable[!dupInd,]
					nuID <- bestTable[,'nuID']
					if (fieldName.match[1] == 'ProbeId' && cut0.probeId) {
						allID <- as.character(as.integer(bestTable[,fieldName.match[1]]))
					} else {
						allID <- bestTable[,fieldName.match[1]]					
					}
					selInputID <- inputID[inputID %in% allID]
					if (fieldName.match[1] == 'nuID') {
						bestTable <- bestTable[, names(bestTable) != 'nuID']
						rownames(bestTable) <- nuID
						mapping <- as.matrix(bestTable[selInputID, ])
					} else {
						names(nuID) <- allID
						mapping <- nuID[selInputID]
					}
					## match the original input IDs
					mapping <- matchInputID(mapping, inputID.bak)
				}
			} else {
				# if (matchField[which.max(matchLen)] != 'nuID') IDType <- 'nuID'
				matchInd <- which(matchLen > 0)
				tableName.match <- allTableNames[matchInd]
				fieldName.match <- matchField[matchInd]
				matchLen.match <- matchLen[matchInd]
				probeNumber <- NULL
				if (idMapping) mapping <- NULL
				for (i in seq(matchInd)) {
					table.i <- dbReadTable(conn, tableName.match[i])
					if (idMapping) {
						## remove duplicated IDs becasue rownames should be unique
						dupInd.i <- duplicated(table.i[,fieldName.match[i]])
						table.i <- table.i[!dupInd.i,]
						nuID.i <- table.i[,'nuID']
						if (fieldName.match[1] == 'ProbeId' && cut0.probeId) {
							allID.i <- as.character(as.integer(table.i[,fieldName.match[i]]))
						} else {
							allID.i <- table.i[,fieldName.match[i]]				
						}
						selInputID.i <- inputID[inputID %in% allID.i]
						if (fieldName.match[which.max(matchLen.match)] == 'nuID') {
							table.i <- table.i[, names(table.i) != 'nuID']
							rownames(table.i) <- nuID.i
							mapping.i <- as.matrix(table.i[selInputID.i, ])
							mapping.i <- matchInputID(mapping.i, inputID.bak)
							mapping <- c(mapping, list(mapping.i))
						} else {
							names(nuID.i) <- allID.i
							mapping.i <- nuID[selInputID.i]
							mapping.i <- matchInputID(mapping.i, inputID.bak)
							mapping <- cbind(mapping, mapping.i)
						}
					}
					probeNumber <- c(probeNumber, nrow(table.i))
				}
				rank1 <- rank(matchLen.match)
				rank2 <- rank(-probeNumber)
				ord <- order(rank1 + rank2/length(matchInd), decreasing=TRUE)
				tableName.match <- tableName.match[ord]
				fieldName.match <- fieldName.match[ord]
				probeNumber <- probeNumber[ord]
				matchLen.match <- matchLen.match[ord]
				if (idMapping) {
					if (is(mapping, 'list')) {
						mapping <- mapping[ord]
						names(mapping) <- tableName.match	
					} else if (is(mapping, 'matrix')) {
						rownames(mapping) <- inputID
						mapping <- mapping[,ord]
						colnames(mapping) <- tableName.match		
					}
				}
			}
			returnList <- list(chipVersion=tableName.match, species=species, IDType=fieldName.match, chipProbeNumber=probeNumber, 
				inputProbeNumber=length(inputID), matchedProbeNumber=matchLen.match)
		}
		if (idMapping)  returnList <- c(returnList, list(idMapping=mapping))		

		return(returnList)
	}
}

