`addNuId2lumi` <-
function(x.lumi, annotationFile=NULL, sep=NULL, lib=NULL, annotationColName=c(sequence='Probe_Sequence', target='Target', probe='ProbeId')) {

    history.submitted <- as.character(Sys.time())

	## check whether the object is nuID annotated.
	exprs <- exprs(x.lumi)
	id <- rownames(exprs)
	if(is.nuID(id[1]) & is.nuID(id[2])) {
		print('The lumiBatch object is already nuID annotated!')
		return(x.lumi)
	}

	## ---------------------------------------
	## identify the Metadata lines 
	if (!is.null(annotationFile)) {
		info <- readLines(file(annotationFile), n=10)    # take the first 10 lines to have a taste

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
					stop('The seperator is not Tab or comma!\n Please sepecify the seperator used in the file!')
				}
			}
		}

		## Read in annotation data
		annotation <- read.table(annotationFile, sep=sep, colClasses="character", header=TRUE, skip=nMetaDataLines,
		 	blank.lines.skip=TRUE, check.names=FALSE, fill=TRUE)

		colnames(annotation) <- toupper(colnames(annotation))
		## Create unique Id based on 50mer sequence
		nuID <- sapply(annotation[, toupper(annotationColName['sequence'])], seq2id)
		## check the TargetID first
		ann_target <- annotation[, toupper(annotationColName['target'])]
		comm_target <- id[id %in% ann_target]
		if (length(comm_target) == 0) {
			## check the ProbeID if id does not match the TargetID
			ann_target <- annotation[, toupper(annotationColName['probe'])]
			comm_target <- id[id %in% ann_target]
			if (length(comm_target) == 0)
				stop('The annotation file does not match the data!')
		} 
		if (length(comm_target) != length(id)) {
			warning('The annotation file does not match the data. Partial ids cannot be replaced!')
		}
		names(nuID) <- ann_target

		newId <- id
		newId[id %in% ann_target] <- nuID[comm_target]
	} else if (!is.null(lib)) {
		if (require(lib, character.only=TRUE)) {
			## check the TargetID first
			newId <- mget(id, get(paste(lib, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
			newId <- unlist(newId)
			if (length(which(!is.na(newId))) == 0) {
				## check the ProbeID if id does not match the TargetID
				newId <- mget(id, get(paste(lib, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
				if (length(which(!is.na(newId))) == 0) {
					targetID <- pData(featureData(x.lumi))$TargetID
					newId <- mget(targetID, get(paste(lib, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
					if (length(which(!is.na(newId))) == 0) stop('The library does not match the data!')
				}
			}
			## Check for the targetIDs cannot be found in the lib.
			## Some known control genes will not be checked.
			naInd <- is.na(newId)
			controlId <- c('lysA','pheA','thrB','trpF', 'bla1','bla2','cat1','cat2','cre1','cre2','e1a1',
			'e1a2','gfp1','gfp2','gst1','gst2','gus1','gus2','lux1','lux2','lysA','neo1',
			'neo2','pheA','thrB','trpF')		
			if (!all(id[naInd] %in% controlId)) {
				if (length(which(naInd)) < 10) {
					warning(paste('Identifiers:', paste(id[naInd], collapse=','), ' cannot be found in the ', lib, '!', sep=''))
				} else {
					warning(paste('Some identifiers cannot be found in the ', lib, '!', sep=''))
				}
			}
			newId[naInd] <- id[naInd]
		} else {
			stop(paste('Annotation library', lib, 'is not installed!'))
		}
	} else {
		stop('Please provide the annotation file or lumi annotation library!')
	}

	if (any(duplicated(newId)))  {
		warning('Duplicated IDs found and were merged!')
		dupId <- unique(newId[duplicated(newId)])
		rmIndex <- NULL
		for (dupId.i in dupId) {
			dupIndex <- which(newId == dupId.i)
			ave.exp <- colMeans(exprs(x.lumi)[dupIndex, ])
			totalBeadNum <- colSums(beadNum(x.lumi)[dupIndex, ])
			maxDetection <- apply(detection(x.lumi), 2, max)

			temp <- colSums(se.exprs(x.lumi)[dupIndex,]^2 * (beadNum(x.lumi)[dupIndex,] - 1))
			temp <- temp / (totalBeadNum - length(dupIndex))
			se.exprs(x.lumi)[dupIndex[1],] <- sqrt(temp * (colSums(1/beadNum(x.lumi)[dupIndex,])))
			exprs(x.lumi)[dupIndex[1],] <- ave.exp
			detection(x.lumi)[dupIndex[1],] <- maxDetection
			beadNum(x.lumi)[dupIndex[1],] <- totalBeadNum
			rmIndex <- c(rmIndex, dupIndex[-1])
		}

		## remove duplicated
		x.lumi <- x.lumi[-rmIndex, ]
		newId <- newId[-rmIndex]
	}

	## update the feature names (probe ids)
	names(newId) <- NULL
	featureNames(x.lumi) <- newId
	## update the feautre data
	featureData <- featureData(x.lumi)
	rownames(pData(featureData)) <- newId
	featureData(x.lumi) <- featureData

	## Add history tracking
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(addNuId2lumi)))
	if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
	lumiVersion <- packageDescription('lumi')$Version
	x.lumi@history<- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))

	return(x.lumi)
}