`lumiR` <-
function(fileName, sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = NULL) 
{
	## the patterns used to grep columns in the BeadStudio output text file 
	## Advanced users can modify this based on the data file
	columnGrepPattern <- c(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='Detection', beadNum='Avg_NBEADS')
	history.submitted <- as.character(Sys.time())
	## set "stringsAsFactors" as FALSE
	oldSetting <- options("stringsAsFactors")[[1]]
	options(stringsAsFactors = FALSE)
	version <- 2

	## ---------------------------------------
	## identify the Metadata lines 
	info <- readLines(file(fileName), n=20)    # take the first 20 lines to have a taste

	## Use "AVG_SIGNAL" as an indicator of Where the metaData stops
	##   intelligently find nMetaDataLines  
	nMetaDataLines <- grep(columnGrepPattern['exprs'], info, ignore.case=TRUE) - 1
    
	if (is.null(sep)) {
	    ## Find out the separator (sep) by taking the first two line of data, and comparing them.
	    ##  we assume it is either "," or "\t".
    	titleLine <- info[nMetaDataLines + 1]
	    dataLine1 <- info[nMetaDataLines + 2]
		dataLine2 <- info[nMetaDataLines + 3]
		sepNum1 <- gregexpr('\t', dataLine1)[[1]]
		sepNum2 <- gregexpr('\t', dataLine2)[[1]]
		if (sepNum1[1] > 0 & length(sepNum1) == length(sepNum2)) {
			sep <- '\t'
		} else {
			sepNum1 <- gregexpr(',', dataLine1)[[1]]
			sepNum2 <- gregexpr(',', dataLine2)[[1]]
			if (sepNum1[1] > 0 & length(sepNum1) == length(sepNum2)) {
				sep <- ','
			} else {
				stop('The seperator is not Tab or comma!\n Please sepecify the seperator used in the file!')
			}
		}
	}
    
	## ---------------------------------------
    # get data info
	if (nMetaDataLines > 0) {
		info <- readLines(file(fileName), n=nMetaDataLines)
		## check the version of the beadStudio output
		markerInd <- grep('^\\[.*\\]', info, ignore.case=TRUE)
		if (length(markerInd) > 0) {
			if (length(grep('^\\[Header\\]', info[markerInd[1]], ignore.case=TRUE)) == 0) 
				warning('The data file may not be in the Illumia BeadStudio output format!')
			if (length(markerInd) > 1) {
				if (length(grep('^\\[.*\\Profile]', info[markerInd[2]], ignore.case=TRUE)) == 0) 
					warning('The data file may not be in the Illumia BeadStudio output format!')
			}
			version <- 3
			info <- info[-markerInd]
		}
		## remove the blanks
		info <- sub("[[:blank:]]+$", "", info)
		info <- gsub(sep, "", info)
		## check the meta info of the file
		if (version == 2) {
			ind <- grep("Illumina Inc. BeadStudio version", info, ignore.case=TRUE)
		} else {
			ind <- grep("BSGX Version", info, ignore.case=TRUE)
		}
		if (length(ind) == 0) 
		    warning("The data file is not in the Illumia BeadStudio output format.")

		## should not be normalized in BeadStudio
		ind <- grep("Normalization", info, ignore.case=TRUE)  # find where is the row index
		if (version == 2) {
			normalization <- strsplit(info, split='=')[[ind]][2]
			normalization <- gsub(pattern=" |,", replace="", normalization) # remove space or ","
		} else {
			normalization <- strsplit(info, split='\t')[[ind]][2]
		}
		if (length(grep("none", normalization, ignore.case=TRUE)) == 0) {
		    warning("The raw data should not be normalized in BeadStudio.")
		}
	} else {
		info <- NULL
	}
    
	allData <- read.table(file=fileName, header=TRUE, sep=sep, skip=nMetaDataLines, row.names=NULL,
		quote='', as.is=TRUE, check.names=FALSE, strip.white=TRUE, comment.char="", fill=TRUE)
	
	## retrieve the possible section line index
	sectionInd <- grep('^\\[.*\\]', allData[,1], ignore.case=TRUE)
    
	if (length(sectionInd) > 0) {
		if (is.na(version)) verion <- 3
		otherData <- allData[sectionInd[1]:nrow(allData), ]
		## we assume the first section is the expression data section
		allData <- allData[1:(sectionInd[1]-1),, drop=FALSE]
		## remove the all NA columns, which can be produced when saved in Excel
		naCol <- apply(allData, 2, function(x) all(is.na(x) | x == ''))
		allData <- allData[,!naCol]
    
		## process otherData
		sectionInd <- sectionInd - sectionInd[1] + 1
		sectionName <- otherData[sectionInd, 1]
    
		## retrieve the control data
		controlInd <- grep('^\\[Control.*\\]', sectionName, ignore.case=TRUE)
		if (length(controlInd) > 1) {
			ind <- grep('^\\[Control probe.*\\]', sectionName[controlInd], ignore.case=TRUE)
			if (length(ind) > 0) {
				controlInd <- controlInd[ind[1]]
			} else {
				controlInd <- controlInd[1]
			}
		}
		if (length(controlInd) > 0) {
			startRow <- sectionInd[controlInd] + 1
			if (length(sectionInd) > controlInd)
				endRow <- sectionInd[controlInd + 1] - 1
			else 
				endRow <- nrow(otherData)
			controlData <- otherData[startRow:endRow,]
			## remove the all NA columns, which can be produced when save in Excel
			naCol <- apply(controlData, 2, function(x) all(is.na(x) | x == ''))
			controlData <- controlData[,!naCol, drop=FALSE]
			colnames(controlData) <- controlData[1,]
			controlData <- controlData[-1,]
		} else {
			controlData <- data.frame()
		}
    
		## retrieve the Sample Table
		summaryInd <- grep('^\\[Sample.*Table\\]', sectionName, ignore.case=TRUE)
		if (length(summaryInd) > 0) {
			startRow <- sectionInd[summaryInd] + 1
			if (length(sectionInd) > summaryInd)
				endRow <- sectionInd[summaryInd + 1] - 1
			else 
				endRow <- nrow(otherData)
			sampleSummary <- otherData[startRow:endRow,]
			## remove the all NA columns, which can be produced when save in Excel
			naCol <- apply(sampleSummary, 2, function(x) all(is.na(x) | x == ''))
			sampleSummary <- sampleSummary[,!naCol]
			colnames(sampleSummary) <- sampleSummary[1,]
			sampleSummary <- sampleSummary[-1,]
		} else {
			sampleSummary <- data.frame()
		}
	} else {
		controlData <- sampleSummary <- data.frame()
	}
	header <- names(allData)
    
	## Get Id. The ProbeID (by default it is the second column) is preferred if provided, 
	# otherwise the TargetID (by default it is the first column) is used.
	targetID <- as.character(as.vector(allData[,1]))
	if (length(grep('ProbeID', header[2], ignore.case=TRUE)) > 0) {
		id <- as.character(as.vector(allData[,2]))
		idName <- header[2]
	} else {
		id <- targetID
		idName <- header[1]
	}
    
	## identify where the signal column exists
	ind <- grep(columnGrepPattern['exprs'], header, ignore.case=TRUE)
	if (length(ind) == 0) stop('Input data format unrecognizable!\nThere is no column name contains "AVG_SIGNAL"!')
	exprs <- as.matrix(allData[,ind])
	if (!is.numeric(exprs[1])) {
		exprs <- matrix(as.numeric(exprs), nrow=nrow(allData))
	} 
	colnames(exprs) <- header[ind]
	## identify where the signal standard deviation column exists 
	ind <- grep(columnGrepPattern['se.exprs'], header, ignore.case=TRUE)
	if (length(ind) == 0) stop('Input data format unrecognizable!\nThere is no column name contains "BEAD_STDEV"!')
	se.exprs <- as.matrix(allData[,ind])
	if (!is.numeric(se.exprs[1])) {
		se.exprs <- matrix(as.numeric(se.exprs), nrow=nrow(allData))
	}
	colnames(se.exprs) <- header[ind]
	## identify the detection columns
	ind <- grep(columnGrepPattern['detection'], header, ignore.case=TRUE)
	if (length(ind) == 0) {
		detection <- NULL
	} else {
		detection <- as.matrix(allData[,ind])
		if (!is.numeric(detection[1])) {
			detection <- matrix(as.numeric(detection), nrow=nrow(allData))
		}
		colnames(detection) <- header[ind]
    
		if (length(grep("Detection Pval", header, ignore.case=TRUE)) == 0) {
			detection <- 1 - detection
		}
	}
	## identify the bead number columns
	ind <- grep(columnGrepPattern['beadNum'], header, ignore.case=TRUE)
	if (length(ind) == 0) {
		beadNum <- NULL
	} else {
	    beadNum <- as.matrix(allData[,ind])
		if (!is.numeric(beadNum[1])) {
			beadNum <- matrix(as.numeric(beadNum), nrow=nrow(allData))
		} 
		colnames(beadNum) <- header[ind]
	}
    
	## check for possible duplicated ids
	dupId <- unique(id[duplicated(id)])
	if (length(dupId) > 0) {
		warning('Duplicated IDs found and were merged!')
		rmInd <- NULL
		for (dupId.i in dupId) {
			selInd.i <- which(id == dupId.i)
			exprs[selInd.i[1],] <- colMeans(exprs[selInd.i,])
			if (is.null(beadNum)) {
				se.exprs[selInd.i[1],] <- colMeans(se.exprs[selInd.i,])
			} else {
				totalBead.i <- colSums(beadNum[selInd.i,])
				beadNum[selInd.i[1],] <- totalBead.i				
				temp <- colSums(se.exprs[selInd.i,]^2 * (beadNum[selInd.i,] - 1))
				temp <- temp / (totalBead.i - length(selInd.i))
				se.exprs[selInd.i[1],] <- sqrt(temp * (colSums(1/beadNum[selInd.i,])))
			}
			if (!is.null(detection)) {
				detection[selInd.i[1],] <- apply(detection[selInd.i,], 2, max)				
			}
			rmInd <- c(rmInd, selInd.i[-1])
		}
		## remove duplicated
		exprs <- exprs[-rmInd,,drop=FALSE]
		se.exprs <- se.exprs[-rmInd,,drop=FALSE]
		id <- id[-rmInd]
		if (!is.null(detection)) detection <- detection[-rmInd,,drop=FALSE]
		if (!is.null(beadNum)) beadNum <- beadNum[-rmInd,,drop=FALSE]
	}
    
	if (na.rm) {
		## remove the probes with all of them as NA
	    keepInd <- apply(is.na(exprs), 1, sum) == 0
	    exprs <- exprs[keepInd,,drop=FALSE]
	    se.exprs <- se.exprs[keepInd,,drop=FALSE]	## bead measurement standard variance
	    beadNum <- beadNum[keepInd,,drop=FALSE]
	    detection <- detection[keepInd,,drop=FALSE]
	    id <- id[keepInd]
		targetID <- targetID[keepInd]
	}
	rownames(exprs) <- rownames(se.exprs) <- rownames(beadNum) <- rownames(detection) <- id
    
	# get sample information
	pattern <- paste('[^[:alnum:]]*', columnGrepPattern['exprs'], '[^[:alnum:]]*', sep='')
	sampleName <-  sub(pattern, '', colnames(exprs), ignore.case=TRUE) 
	sampleNameInfo <- strsplit(sampleName, split="_")
	sampleID <- NULL
	label <- NULL
	temp <- lapply(sampleNameInfo, function(x) {sampleID <<- c(sampleID, x[1]); label <<- c(label, x[2])})
	if (any(is.na(label))) {
		sampleID <- sampleName
		label <- sampleName
	}
    
	## reportInfo save the id information
	if (!is.null(detection)) {
		presentCount <- apply(detection, 1, function(x) sum(x <= detectionTh))
		reporterInfo <- data.frame(id, presentCount)
		names(reporterInfo) <- c(idName, 'presentCount')
		varMetadata <- data.frame(labelDescription=c('The Illumina microarray identifier', 
			'The number of detectable measurements of the gene'))
		rownames(varMetadata) <- c(idName, 'presentCount')
	} else {
		reporterInfo <- data.frame(id)
		names(reporterInfo) <- idName
		varMetadata <- data.frame(labelDescription='The Illumina microarray identifier')
		rownames(varMetadata) <- idName
	}
	## if ProbeID is used as id, then also keep the TargetID information in the featureData
	if (idName != header[1]) {
		reporterInfo <- data.frame(reporterInfo, TargetID=targetID)
		varMetadata <- rbind(varMetadata, data.frame(labelDescription='The Illumina TargetID'))
		rownames(varMetadata)[nrow(varMetadata)] <- 'TargetID'
	}
	rownames(reporterInfo) <- id
	featureData <- new("AnnotatedDataFrame", data=reporterInfo, varMetadata=varMetadata)
    
	## set the colnames as the label or sampleName
	if (length(unique(label)) == length(label) & length(label) > 0) {
		colName <- label
	} else {
		colName <- sampleName
	}
	
	## check the dimensions of the input data
	if (ncol(exprs) == ncol(se.exprs)) {
		colnames(exprs) <- colnames(se.exprs) <- colName
	} else {
		stop('Different column numbers of exprs and se.exprs! Please check the input data format.')
	}
	if (ncol(beadNum) == length(colName)) {
		colnames(beadNum) <- colName
	} else {
		warning('The number of beadNum columns does not match! Please check the input data format.')
		if (ncol(beadNum) > length(colName)) beadNum <- beadNum[,1:length(colName)]
		if (ncol(beadNum) < length(colName)) {
			for (i in 1:(length(colName) - ncol(beadNum)))
				beadNum <- cbind(beadNum, rep(NA, nrow(beadNum)))
		}
	}
	if (ncol(detection) == length(colName)) {
		colnames(detection) <- colName
	} else {
		warning('The number of detection columns does not match! Please check the input data format.')
		if (ncol(detection) > length(colName)) beadNum <- beadNum[,1:length(colName)]
		if (ncol(detection) < length(colName)) {
			for (i in 1:(length(colName) - ncol(detection)))
				detection <- cbind(detection, rep(NA, nrow(detection)))
		}
	}

	## produce the phenoData object
	pData <- data.frame(sampleID=sampleID, label=label)
	rownames(pData) <- colName
	#pdata <- new("phenoData", pData=pData, varLabels=list('sampleID', 'label'))
	varMetadata <- data.frame(labelDescription=c('The unique Illumina microarray Id', 
		'The label of the sample'))
	rownames(varMetadata) <- c('sampleID', 'label')
	pdata <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)

    x.lumi <- new("LumiBatch", exprs=exprs, se.exprs=se.exprs, beadNum=beadNum, detection=detection, 
              featureData=featureData,  phenoData=pdata)
	x.lumi@controlData <- controlData
	x.lumi@QC <- list(BeadStudioSummary=sampleSummary)
	
	if (!is.null(info)) {
		info <- gsub('\t+', '\t', info)
	}
	experimentData(x.lumi)@other <- list(info)

    # history tracking
    history.finished <- as.character(Sys.time())
	#history.command <- match.call()
    history.command <- capture.output(print(match.call(lumiR)))  
	
	## replace with the real file name
	if (length(grep(',', history.command)) > 0) {
		history.command <- sub('\\(.+,', paste('("', fileName, '",', sep=''), history.command)
	} else {
		history.command <- sub('\\(.+\\)', paste('("', fileName, '")', sep=''), history.command)
	}

	lumiVersion <- packageDescription('lumi')$Version
	x.lumi@history<- data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion)

	## initialize the QC slot in the LumiBatch object
	x.lumi <- lumiQ(x.lumi)

	## Add nuID if the annotation library is provided
	if (!is.null(lib))  x.lumi <- addNuId2lumi(x.lumi, lib=lib)

	## resume the old settings
	options(stringsAsFactors = oldSetting)
    
    return(x.lumi)
}

