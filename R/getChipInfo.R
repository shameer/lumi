`getChipInfo` <-
function(x, lib.mapping=NULL, species=c('Human', 'Mouse', 'Rat', 'Unknown'), idMapping=FALSE, returnAllMatches=FALSE, verbose=TRUE) {
	if (is(x, 'ExpressionSet')) {
		inputID <- featureNames(x)
	} else if (is(x, 'matrix') || is(x, 'data.frame')) {
		inputID <- rownames(x)		
	} else {
		inputID <- x
	}
	inputID <- unique(inputID)
	species <- match.arg(species)
	if (species == 'Unknown' && is.null(lib.mapping)) {
		human.match <- getChipInfo(inputID, species='Human', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		mouse.match <- getChipInfo(inputID, species='Mouse', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		rat.match <- getChipInfo(inputID, species='Rat', idMapping=idMapping, returnAllMatches=returnAllMatches, verbose=FALSE)
		bestChip <- which.max(c(length(human.match$matchProbeNumber), length(mouse.match$matchProbeNumber), length(rat.match$matchProbeNumber)))
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

		allTableNames <- dbListTables(conn)
		allTableNames <- allTableNames[!(allTableNames %in% c('nuID_MappingInfo', 'metadata'))]
		lenID <- length(inputID)
		matchLen <- NULL
		matchField <- NULL
		tableLen <- NULL
		for (i in seq(allTableNames)) {
			tableName.i <- allTableNames[i]
			table.i <- dbReadTable(conn, tableName.i)
			fieldNames.i <- names(table.i)
			## Not include the Symbol column during matching the IDs.
			fieldNames.i <- fieldNames.i[!(fieldNames.i %in% c('Symbol'))]
			len.i <- NULL
			for (j in seq(fieldNames.i)) {
				field.ij <- as.character(table.i[,j])
				len.ij <- length(which(inputID %in% field.ij))
				len.i <- c(len.i, len.ij)
			}
			max.ind.i <- which.max(len.i)
			matchLen <- c(matchLen, len.i[max.ind.i])
			matchField <- c(matchField, fieldNames.i[max.ind.i])
			tableLen <- c(tableLen, nrow(table.i))
		}
		bestMatchLen <- max(matchLen)
		if (bestMatchLen < lenID && verbose) warning('Some input IDs can not be matched!\n')
		if (bestMatchLen == 0) {
			return(list(chipVersion=NULL, species=species, IDType=NULL, chipProbeNumber=NULL, 
				matchProbeNumber=bestMatchLen))
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
				if (fieldName.match == 'nuID') {
					bestTable <- bestTable[, names(bestTable) != 'nuID']
					rownames(bestTable) <- nuID
					mapping <- bestTable[inputID, ]
				} else {
					names(nuID) <- bestTable[,fieldName.match[1]]
					mapping <- nuID[inputID]
				}
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
					if (fieldName.match[which.max(matchLen.match)] == 'nuID') {
						table.i <- table.i[, names(bestTable) != 'nuID']
						rownames(table.i) <- nuID.i
						mapping.i <- table.i[inputID, ]
						mapping <- c(mapping, list(mapping.i))
					} else {
						names(nuID.i) <- table.i[,fieldName.match[i]]
						mapping.i <- nuID[inputID]
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
			matchProbeNumber <- matchLen.match[ord]
			if (idMapping) {
				rownames(mapping) <- inputID
				if (is(mapping, 'list')) {
					mapping <- mapping[ord]
					names(mapping) <- tableName.match	
				} else if (is(mapping, 'matrix')) {
					mapping <- mapping[,ord]
					colnames(mapping) <- tableName.match		
				}
			}
		}
		returnList <- list(chipVersion=tableName.match, species=species, IDType=fieldName.match, chipProbeNumber=probeNumber, 
			matchProbeNumber=matchLen.match)
		if (idMapping)  returnList <- c(returnList, list(idMapping=mapping))		

		return(returnList)
	}
}

