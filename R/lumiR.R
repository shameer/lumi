`lumiR` <-
function(fileName, sep=NULL, detectionTh=0.99, na.rm=TRUE) {
    history.submitted <- as.character(Sys.time())

	## ---------------------------------------
	## identify the Metadata lines 
	info <- readLines(file(fileName), n=20)    # take the first 20 lines to have a taste
    
	## Use "AVG_SIGNAL" as an indicator of Where the metaData stops
	##   intelligently find nMetaDataLines  
	nMetaDataLines <- grep("AVG_SIGNAL", info, ignore.case=TRUE) - 1
    
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
    # get metaData
    info <- readLines(file(fileName), n=nMetaDataLines)

	## check the format
	ind <- grep("Illumina Inc. BeadStudio version", info, ignore.case=TRUE) 
	## find where is the row index
	if (length(ind) == 0) {
	    warning("The data file is not in the Illumia BeadStudio output format.")
	}
    
	## must not normalized before
	ind <- grep("Normalization", info, ignore.case=TRUE)  # find where is the row index
	normalization <- strsplit(info, split='=')[[ind]][2]
	normalization <- gsub(pattern=" |,", replace="", normalization) # remove space or ","
	if (length(grep("none", normalization, ignore.case=TRUE)) == 0) {
	    warning("The raw data should not be normalized in BeadStudio.")
	}
	
    allData <- read.table(file=fileName, header = TRUE, sep = sep, quote = "\"",
               dec = ".", as.is = TRUE,
               na.strings = "NaN", colClasses = NA, 	#nrows = 60000,
               skip = nMetaDataLines, check.names = FALSE, 
               strip.white = FALSE, blank.lines.skip = TRUE,
               comment.char = "", allowEscapes = FALSE, flush = FALSE)
	header <- names(allData)

	# get target Id
	targetID <- as.character(as.vector(allData[,1]))
	rownames(allData) <- targetID
	id <- targetID

	## identify where the signal column exists
	ind <- grep("AVG_SIGNAL", header, ignore.case=TRUE)
	if (length(ind) == 0) stop('Input data format unrecognizable!\nThere is no column name contains "AVG_SIGNAL"!')
	exprs <- as.matrix(allData[,ind])
	## identify where the signal standard deviation column exists   
	ind <- grep("BEAD_STDEV", header, ignore.case=TRUE)
	if (length(ind) == 0) stop('Input data format unrecognizable!\nThere is no column name contains "BEAD_STDEV"!')
	se.exprs <- as.matrix(allData[,ind])
	ind <- grep("Detection", header, ignore.case=TRUE)
	if (length(ind) == 0) {
		detection <- NULL
	} else {
		detection <- as.matrix(allData[,ind])
	}
	ind <- grep("Avg_NBEADS", header, ignore.case=TRUE)
	if (length(ind) == 0) {
		beadNum <- NULL
	} else {
	    beadNum <- as.matrix(allData[,ind])
	}
	
	## check for possible duplicated ids
	dupId <- unique(id[duplicated(id)])
	if (length(dupId) > 0) {
		warning('Duplicated Ids found!')
		for (dupId.i in dupId) {
			selInd.i <- id == dupId.i
			exprs[dupId.i,] <- colMeans(exprs[selInd.i,])
			if (is.null(beadNum)) {
				se.exprs[dupId.i,] <- colMeans(se.exprs[selInd.i,])
			} else {
				totalBead.i <- colSums(beadNum[selInd.i,])
				beadNum[dupId.i,] <- totalBead.i				
				temp <- colSums(se.exprs[selInd.i,]^2 * (beadNum[selInd.i,] - 1))
				temp <- temp / (totalBead.i - length(which(selInd.i)))
				se.exprs[dupId.i,] <- sqrt(temp * (colSums(1/beadNum[selInd.i,])))
			}
			if (!is.null(detection)) {
				detection[dupId.i,] <- apply(detection[selInd.i,], 2, max)				
			}
		}
		## remove duplicated
		dupId <- duplicated(id)
		exprs <- exprs[!dupId,]
		se.exprs <- se.exprs[!dupId,]
		id <- id[!dupId]
		if (!is.null(detection)) detection <- detection[!dupId,]
		if (!is.null(beadNum)) beadNum <- beadNum[!dupId,]
	}
	
    if (na.rm) {
		## remove the probes with all of them as NA
        keepInd <- apply(is.na(exprs), 1, sum) == 0
        exprs <- exprs[keepInd,]
        se.exprs <- se.exprs[keepInd,]	## bead measurement standard variance
        beadNum <- beadNum[keepInd,]
        detection <- detection[keepInd,]
        id <- id[keepInd]
		targetID <- targetID[keepInd]
    }
    
    # get sample information
	sampleName <-  sub('AVG_SIGNAL.', '', colnames(exprs), ignore.case=TRUE) 
    sampleNameInfo <- strsplit(sampleName, split="_")
	sampleID <- NULL
	label <- NULL
	temp <- lapply(sampleNameInfo, function(x) {sampleID <<- c(sampleID, x[1]); label <<- c(label, x[2])})
	
	## reportInfo save the id information
	if (!is.null(detection)) {
		presentCount <- apply(detection, 1, function(x) sum(x >= detectionTh))
		reporterInfo <- data.frame(targetID=targetID, presentCount=presentCount)	
		varMetadata <- data.frame(labelDescription=c('The Illumina microarray gene ID', 
			'The number of detectable measurements of the gene'))
		rownames(varMetadata) <- c('targetID', 'presentCount')
	} else {
		reporterInfo <- data.frame(targetID=targetID)
		varMetadata <- data.frame(labelDescription='The Illumina microarray gene ID')
		rownames(varMetadata) <- 'targetID'
	}
	rownames(reporterInfo) <- id
	featureData <- new("AnnotatedDataFrame", data=reporterInfo, varMetadata=varMetadata)
	
    ## set the colnames as the label or sampleName
	if (length(unique(label)) == length(label) & length(label) > 0) {
		colName <- label
	} else {
		colName <- sampleName
	}
	colnames(exprs) <- colnames(se.exprs) <- colnames(beadNum) <- colnames(detection) <- colName

	## produce the phenoData object
	pData <- data.frame(sampleID=sampleID, label=label)
	rownames(pData) <- colName
	#pdata <- new("phenoData", pData=pData, varLabels=list('sampleID', 'label'))
	varMetadata <- data.frame(labelDescription=c('The unique Illumina microarray Id', 
		'The label of the sample'))
	rownames(varMetadata) <- c('sampleID', 'label')
	pdata <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)

    # history tracking
    history.finished <- as.character(Sys.time())
	#history.command <- match.call()
    history.command <- capture.output(print(match.call(lumiR)))  
	
	## replace with the real file name
	#history.command <- sub('[^(]fileName', paste('"', fileName, '"', sep=''), history.command)
	if (length(grep(',', history.command)) > 0) {
		history.command <- sub('\\(.+,', paste('("', fileName, '",', sep=''), history.command)
	} else {
		history.command <- sub('\\(.+\\)', paste('("', fileName, '")', sep=''), history.command)
	}

    x.lumi <- new("LumiBatch", exprs=exprs, se.exprs=se.exprs, beadNum=beadNum, detection=detection, 
              featureData=featureData,  phenoData=pdata)

	info <- gsub('\t+', '\t', info)
	experimentData(x.lumi)@other <- list(info)
    x.lumi@history<- rbind(x.lumi@history,
                       data.frame(submitted=history.submitted, finished=history.finished, command=history.command))
    
    return(x.lumi)
}

