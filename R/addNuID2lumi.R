`addNuID2lumi` <-
function(x.lumi, annotationFile=NULL, sep=NULL, lib.mapping=NULL, annotationColName=c(sequence='Probe_Sequence', target='ILMN_Gene', probe='Probe_Id'), verbose=TRUE) {

    history.submitted <- as.character(Sys.time())

	## check whether the object is nuID annotated.
	exprs <- exprs(x.lumi)
	id <- rownames(exprs)
	if(all(sapply(id[1:50], is.nuID))) {
		if (verbose) cat('The lumiBatch object is already nuID annotated!\n')
		return(x.lumi)
	}
	if (!is.null(lib.mapping)) {
		if (length(grep('lumi.*\\.db', lib.mapping)) == 0 && length(grep('lumi.*IDMapping', lib.mapping)) == 0) {
			warning(paste(lib.mapping, 'does not include nuID conversion information!\n'))
			lib.mapping <- NULL
		}
	}

	newId <- id

	## ---------------------------------------
	## First check whether probe sequence information is available 
	annotation <- pData(featureData(x.lumi))
	names(annotation) <- toupper(names(annotation))
	if (toupper(annotationColName['sequence']) %in% names(annotation)) {
		sequence <- annotation[, toupper(annotationColName['sequence'])]
		cat('Directly converting probe sequence to nuIDs ...\n')
		newId <- sapply(sequence, seq2id)
		names(newId) <- id				
	} else if (!is.null(annotationFile)) {
		## identify the Metadata lines 
		info <- readLines(annotationFile, n=10)    # take the first 10 lines to have a taste

		## Use annotationColName[1] as an indicator of Where the metaData stops
		##   intelligently find nMetaDataLines  
		nMetaDataLines <- grep(annotationColName[1], info, ignore.case=TRUE) - 1

		if (is.null(sep)) {
		    ## Find out the separator (sep) by taking the first two line of data, and comparing them.
		    ##  we assume it is either "," or "\t".
	    	titleLine <- info[nMetaDataLines + 1]
			sepNum <- regexpr('\t', titleLine)
			if (sepNum >= 2) {
				sep <- '\t'
			} else {
				sepNum <- regexpr(',', titleLine)
				if (sepNum >= 2) {
					sep <- ','
				} else {
					stop('The seperator is not Tab or comma!\n Please sepecify the seperator used in the file!\n')
				}
			}
		}

		dataLine1 <- strsplit(info[nMetaDataLines + 2], sep)[[1]]
		quoteCount1 <- gregexpr('"', dataLine1[1])[[1]]
		quoteCount2 <- gregexpr('\'', dataLine1[1])[[1]]
		quote <- ''
		if (sep == ',') quote <- '"'
		if (length(quoteCount1) == 2) {
			quote <- '"'
		} else if (length(quoteCount2) == 2) {
			quote <- '\''
		}

		## Read in annotation data
		annotation <- read.table(annotationFile, sep=sep, colClasses="character", header=TRUE, skip=nMetaDataLines,
		 	blank.lines.skip=TRUE, row.names=NULL, check.names=FALSE, quote=quote, comment.char="", strip.white=TRUE, fill=TRUE)
		
		allColName <- toupper(colnames(annotation))
		colnames(annotation) <- allColName
		## Create unique Id based on 50mer sequence
		if (toupper(annotationColName['sequence']) %in% allColName) {
			nuID <- sapply(annotation[, toupper(annotationColName['sequence'])], seq2id)			
		} else {
			stop('The "sequence" column cannot be found!\nPlease check the "annotationColName" of "sequence"!\n')
		}
		## check the TargetID first
		comm_target <- NULL
		if (toupper(annotationColName['target']) %in% allColName) {
			ann_target <- annotation[, toupper(annotationColName['target'])]
			comm_target <- id[id %in% ann_target]
		}
		if (length(comm_target) == 0) {
			## check the ProbeID if id does not match the TargetID
			if (toupper(annotationColName['probe']) %in% allColName) {
				ann_target <- annotation[, toupper(annotationColName['probe'])]
				comm_target <- id[id %in% ann_target]
				if (length(comm_target) == 0) {
					width <- nchar(ann_target[1])
					id <- formatC(as.numeric(id), width=width, flag='0', format='d')
					comm_target <- id[id %in% ann_target]
					if (length(comm_target) == 0) stop('The annotation file does not match the data!\n')
				}
			}
		} 
		if (length(comm_target) < length(id)) {
			diffId.len <- length(id) - length(comm_target)
			warning(paste('The annotation file does not match the data.',  diffId.len, 'ids cannot be replaced!\n'))
		}
		names(nuID) <- ann_target

		newId <- id
		newId[id %in% ann_target] <- nuID[comm_target]
	} else if (!is.null(lib.mapping)) {
		if (!require(lib.mapping, character.only=TRUE)) stop(paste('Annotation library', lib.mapping, 'is not installed!\n'))
		controlId <- c('lysA','pheA','thrB','trpF', 'bla1','bla2','cat1','cat2','cre1','cre2','e1a1',
		'e1a2','gfp1','gfp2','gst1','gst2','gus1','gus2','lux1','lux2','lysA','neo1',
		'neo2','pheA','thrB','trpF')
		if (length(grep('IDMapping', lib.mapping)) == 0) {
			## check the TargetID first
			newId <- mget(id, get(paste(lib.mapping, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
			newId <- unlist(newId)
			if (length(which(!is.na(newId))) == 0) {
				usingTargetID <- FALSE
				## check the ProbeID if id does not match the TargetID
				newId <- mget(id, get(paste(lib.mapping, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
				if (length(which(!is.na(newId))) == 0) {
					width <- nchar(ls(envir=get(paste(lib.mapping, 'PROBEID2NUID', sep=''), mode='environment'))[1])
					id <- formatC(as.numeric(id), width=width, flag='0', format='d')
					newId <- mget(id, get(paste(lib.mapping, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
					if (length(which(!is.na(newId))) == 0) {
						targetID <- pData(featureData(x.lumi))$TargetID
						newId <- mget(targetID, get(paste(lib.mapping, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
						if (length(which(!is.na(newId))) == 0) stop('The library does not match the data!\n')
					}
				} 
			} else {
				usingTargetID <- TRUE
			}
			## Check for the targetIDs cannot be found in the lib.mapping.
			## Some known control genes will not be checked.
			naInd <- which(is.na(newId))
			if (!usingTargetID) {
				TargetID <- featureData(x.lumi)$TargetID
				if (is.null(TargetID)) {
					if (!all(TargetID[naInd] %in% controlId)) {
						if (length(naInd) < 10) {
							warning(paste('Identifiers:', paste(TargetID[naInd], collapse=','), ' cannot be found in the ', lib.mapping, '!\n', sep=''))
						} else {
							warning(paste('Some identifiers cannot be found in the ', lib.mapping, '!\n', sep=''))
						}
					}
				}
			} else if (!all(id[naInd] %in% controlId)) {
				if (length(naInd) < 10) {
					warning(paste('Identifiers:', paste(id[naInd], collapse=','), ' cannot be found in the ', lib.mapping, '!\n', sep=''))
				} else {
					warning(paste('Some identifiers cannot be found in the ', lib.mapping, '!\n', sep=''))
				}
			}
			newId[naInd] <- id[naInd]
		} else {
			chipInfo <- getChipInfo(featureNames(x.lumi), lib.mapping=lib.mapping, idMapping=TRUE, verbose=FALSE)
			newId <- chipInfo$idMapping[,'nuID']
			naInd <- which(is.na(newId))
			if (length(naInd) > 0) {
				newId[naInd] <- id[naInd]
				if (!all(id[naInd] %in% controlId))	{
					if (length(naInd) > 500) {
						warning(paste('More than 500 identifiers cannot be found in the library:', lib.mapping, 
								'!\n The provided library might be wrong or outdated!\n', sep=''))
					} else if (length(naInd) < 10) {
						warning(paste('Identifiers:', paste(id[naInd], collapse=','), ' cannot be found in the ', lib.mapping, '!\n', sep=''))
					} else {
						warning(paste('Some identifiers cannot be found in the ', lib.mapping, '!\n', sep=''))
					}
				}
			}
			selInfoName <- names(chipInfo)
			selInfoName <- selInfoName[selInfoName != 'idMapping']
			chipInfo.print <- paste(selInfoName, unlist(chipInfo[selInfoName]), sep=': ')
			notes(x.lumi) <- c(notes(x.lumi), list('Chip Information'=chipInfo.print))

			annotation(x.lumi) <- switch(chipInfo$species,
					'Human'='lumiHumanAll.db',
					'Mouse'='lumiMouseAll.db',
					'Rat'='lumiRatAll.db')
		}
	} else {
		cat('Please provide Illumina ID Mapping library!\n')
	}
	if (all(newId == id)) {
		conversion <- FALSE
	} else {
		conversion <- TRUE
	}

	if (any(duplicated(newId)))  {
		if (verbose) cat('Duplicated IDs found and were merged!\n')
		dupId <- unique(newId[duplicated(newId)])
		## determine whether the detection p-value close to 0 or 1 is significant
		if (!is.null(detection(x.lumi))) {
			detect.low <- exprs[which.max(detection(x.lumi)[,1]), 1]
			detect.high <- exprs[which.min(detection(x.lumi)[,1]), 1]
		}
		
		rmIndex <- NULL
		for (dupId.i in dupId) {
			dupIndex <- which(newId == dupId.i)
			ave.exp <- colMeans(exprs(x.lumi)[dupIndex,, drop=FALSE])
			exprs(x.lumi)[dupIndex[1],] <- ave.exp
			if (is(x.lumi, 'LumiBatch') && !is.null(beadNum(x.lumi)) && !is.null(detection(x.lumi))) {
				totalBeadNum <- colSums(beadNum(x.lumi)[dupIndex, , drop=FALSE])
				if (detect.low < detect.high) {
					maxDetection <- apply(detection(x.lumi), 2, min)
				} else {
					maxDetection <- apply(detection(x.lumi), 2, max)
				}

				temp <- colSums(se.exprs(x.lumi)[dupIndex, , drop=FALSE]^2 * (beadNum(x.lumi)[dupIndex,, drop=FALSE] - 1))
				temp <- temp / (totalBeadNum - length(dupIndex))
				se.exprs(x.lumi)[dupIndex[1],] <- sqrt(temp * (colSums(1/beadNum(x.lumi)[dupIndex,, drop=FALSE])))
				detection(x.lumi)[dupIndex[1],] <- maxDetection
				beadNum(x.lumi)[dupIndex[1],] <- totalBeadNum
			}
			rmIndex <- c(rmIndex, dupIndex[-1])
		}
		## remove duplicated
		x.lumi <- x.lumi[-rmIndex, ]
		newId <- newId[-rmIndex]
	}

	## update the feature names (probe ids)
	if (conversion) {
		names(newId) <- NULL
		featureNames(x.lumi) <- newId
		## update the feautre data
		featureData <- featureData(x.lumi)
		rownames(pData(featureData)) <- newId
		if (!is.null(pData(featureData)$'PROBE_SEQUENCE')) pData(featureData)$'PROBE_SEQUENCE' <- NULL
		featureData(x.lumi) <- featureData

		## Add history tracking
		if (is(x.lumi, 'LumiBatch')) {
			history.finished <- as.character(Sys.time())
			# history.command <- capture.output(print(match.call(addNuID2lumi)))
			history.command <- paste(deparse(match.call(addNuID2lumi)), collapse='') 
			if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
			lumiVersion <- packageDescription('lumi')$Version
			x.lumi@history<- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
					finished=history.finished, command=history.command, lumiVersion=lumiVersion))
		}
	}

	return(x.lumi)
}

`addNuId2lumi` <-
function(x.lumi, annotationFile=NULL, sep=NULL, lib.mapping=NULL, annotationColName=c(sequence='Probe_Sequence', target='Target', probe='Probe_Id'), verbose=TRUE) {
	cat('Function addNuId2lumi is deprecated!\n Please use addNuID2lumi instead!\n')
	addNuID2lumi(x.lumi, annotationFile, sep, lib.mapping, annotationColName, verbose)
}