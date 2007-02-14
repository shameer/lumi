`addNuId2lumi` <-
function(x.lumi, annotationFile=NULL, sep=NULL, lib=NULL, annotationColName=c(sequence='Probe_Sequence', target='Target')) {

    history.submitted <- as.character(Sys.time())

	## check whether the object is nuID annotated.
	exprs <- exprs(x.lumi)
	targetID <- rownames(exprs)
	if(is.nuID(targetID[1]) & is.nuID(targetID[2])) {
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
		id <- sapply(annotation[, toupper(annotationColName['sequence'])], seq2id)
		ann_target <- annotation[, toupper(annotationColName['target'])]
		names(id) <- ann_target
		comm_target <- targetID[targetID %in% ann_target]
		if (length(comm_target) == 0) {
			stop('The annotation file does not match the data!')
		} else if (length(comm_target) != length(targetID)) {
			warning('The annotation file does not match the data. Partial ids cannot be replaced!')
		}

		newId <- targetID
		newId[targetID %in% ann_target] <- id[comm_target]
	} else if (!is.null(lib)) {
		if (require(lib, character.only=TRUE)) {
			newId <- mget(targetID, get(paste(lib, 'TARGETID2NUID', sep=''), mode='environment'), ifnotfound=NA)
			newId <- unlist(newId)
			if (length(which(!is.na(newId))) == 0) {
				newId <- mget(targetID, get(paste(lib, 'PROBEID2NUID', sep=''), mode='environment'), ifnotfound=NA)
				if (length(which(!is.na(newId))) == 0) stop('The library does not match the data!')
			}
			## Check for the targetIDs cannot be found in the lib.
			## Some known control genes will not be checked.
			naInd <- is.na(newId)
			controlId <- c('lysA','pheA','thrB','trpF', 'bla1','bla2','cat1','cat2','cre1','cre2','e1a1',
			'e1a2','gfp1','gfp2','gst1','gst2','gus1','gus2','lux1','lux2','lysA','neo1',
			'neo2','pheA','thrB','trpF')		
			if (!all(targetID[naInd] %in% controlId)) {
				if (length(which(naInd)) < 10) {
					warning(paste('TargetIDs:', paste(targetID[naInd], collapse=','), ' cannot be found in the ', lib, '!', sep=''))					
				} else {
					warning(paste('Some TargetIDs cannot be found in the ', lib, '!', sep=''))					
				}
			}
			newId[naInd] <- targetID[naInd]
		} else {
			stop(paste('Annotation library', lib, 'is not installed!'))
		}
	} else {
		stop('Please provide the annotation file or lumi annotation library!')
	}

	if (any(duplicated(newId)))  {
		warning('Duplicated IDs found!!!')
		dupId <- newId[duplicated(newId)]
		rmIndex <- NULL
		for (dupId.i in dupId) {
			dupIndex <- which(newId == dupId.i)
			ave.exp <- colMeans(exprs(x.lumi)[dupIndex, ])
			beadNum <- colSums(beadNum(x.lumi)[dupIndex, ])
			detection <- apply(detection(x.lumi), 2, max)
			# varExp <- (se.exprs(x.lumi)[dupIndex, ])^2 * beadNum(x.lumi)[dupIndex, ] / 
			#	matrix(rep(beadNum, each=length(dupIndex)), nrow=length(dupIndex))
			# stdExp <- apply(varExp, 2, function(x) sqrt(sum(x)))
			# se.exprs(x.lumi)[dupIndex[1],] <- stdExp
			temp <- colSums(se.exprs[dupIndex,]^2 * (beadNum[dupIndex,] - 1))
			temp <- temp / (beadNum - length(dupIndex))
			se.exprs(x.lumi)[dupIndex[1],] <- sqrt(temp * (colSums(1/beadNum[dupIndex,])))
			exprs(x.lumi)[dupIndex[1],] <- ave.exp				
			detection(x.lumi)[dupIndex[1],] <- detection
			beadNum(x.lumi)[dupIndex[1],] <- beadNum
			rmIndex <- c(rmIndex, dupIndex[-1])
		}

		## remove duplicated
		x.lumi <- x.lumi[-rmIndex, ]
		newId <- newId[-rmIndex]
	}		

	## update the feature names (probe ids)
	featureNames(x.lumi) <- newId
	## update the feautre data
	featureData <- featureData(x.lumi)
	rownames(pData(featureData)) <- newId
	featureData(x.lumi) <- featureData

	## Add history tracking
    history.finished <- as.character(Sys.time())
    history.command <- capture.output(print(match.call(addNuId2lumi)))  
    x.lumi@history<- rbind(x.lumi@history,
                       data.frame(submitted=history.submitted, finished=history.finished, command=history.command))

	return(x.lumi)
}