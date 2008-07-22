lumiR.batch <- function(fileList, convertNuID = TRUE, lib.mapping = NULL, detectionTh = 0.01, QC = TRUE, transform = c('none', 'vst', 'log2', 'cubicRoot'), sampleInfoFile = NULL, verbose = TRUE, ...) {

	oldDir <- getwd()
	dirMode <- FALSE
	transform <- match.arg(transform)
	if (file.exists(fileList[1])) {
		if (file.info(fileList[1])[1,'isdir']) {
			dirMode <- TRUE
			setwd(fileList)
		}
	}
	if (dirMode && length(fileList) == 1) {
		fileList <- dir(fileList, pattern='.csv')
		if (length(fileList) == 0) stop('No data files were found!\n')
	}
	
	history.submitted <- as.character(Sys.time())

	if (verbose) {
		cat('Inputting the data ...\n')
		if (transform != 'none') cat(paste('Transformation', transform, 'will be performed for each data file ...\n'))
	}
	for (i in 1:length(fileList)) {
		file.i <- fileList[i]
		if (transform != 'none') {
			x.lumi.i <- lumiR(file.i, parseColumnName=FALSE, convertNuID = FALSE, verbose=FALSE, ...)
			x.lumi.i <- lumiT(x.lumi.i, method=transform, simpleOutput=TRUE, verbose=FALSE)
		} else {
			x.lumi.i <- lumiR(file.i, parseColumnName=FALSE, convertNuID = FALSE, QC = FALSE, verbose=FALSE, ...)
		}
		if (i == 1) {
			x.lumi <- x.lumi.i
		} else {
			x.lumi <- combine(x.lumi, x.lumi.i)
		}
	}
	if (!convertNuID) lib.mapping <- NULL	
	if (!is.null(lib.mapping) || convertNuID) {
		if (verbose) cat('\nAdding nuID to the data ...\n')
		x.lumi <- addNuID2lumi(x.lumi, lib.mapping=lib.mapping)
	}

	if (!is.null(sampleInfoFile)) {

		if (is.character(sampleInfoFile) || class(sampleInfoFile)[1] == 'file') {
			if (file.exists(sampleInfoFile)) {
				sampleInfo <- read.table(sampleInfoFile, head=TRUE, sep='\t', colClasses='character', comment='')
			} else {
				warning('The provided sampleInfoFile does not exist\n!')
				setwd(oldDir)
				return(x.lumi)
			}
		} else if (is.data.frame(sampleInfoFile)) {
			sampleInfo <- sampleInfoFile
		}
		## force the names to be capitalized
		colName <- toupper(names(sampleInfo))
		ind <- grep('ID$', colName, ignore.case=TRUE)
		if (length(ind) == 0) {
			ID <- sampleInfo[,1]
			if (any(duplicated(ID))) {
				warning('In sampleInfoFile, the ID column is required or the first column should be unique!\n')
				setwd(oldDir)
				return(x.lumi)
			}
			ind <- 1
		} else {
			ID <- sampleInfo[, ind[1]]			
		}
		rownames(sampleInfo) <- ID
		colnames(sampleInfo)[ind[1]] <- 'sampleID'

		sampleName <- sampleNames(x.lumi)
		ID <- ID[ID %in% sampleName]
		if (nrow(sampleInfo) != length(ID)) {
			warning('Some IDs provided in the sampleInfoFile do not exist the data file!\n')
			if (length(ID) == 0) {
				stop('The IDs provided in the sampleInfoFile do not match the data file!\n')
			} 
			
		} 
		x.lumi <- x.lumi[, ID]
		
		if (is.null(pData(phenoData(x.lumi)))) {
			pData <- sampleInfo[ID,]			
		} else {
			pData <- data.frame(pData(phenoData(x.lumi))[!(toupper(names(pData(phenoData(x.lumi)))) %in% c(toupper(names(sampleInfo)), 'ID', 'SAMPLEID'))], sampleInfo[ID,])
		}
		label <- sampleInfo[ID, colName == 'LABEL']
		if (length(label) == length(ID)) {
			rownames(pData) <- label
			sampleNames(x.lumi) <- label
		}

		pdata <- new("AnnotatedDataFrame", data=pData)
		phenoData(x.lumi) <- pdata
	} 
	
	## ----------------------------------------
	# Add history track
	if (is(x.lumi, 'LumiBatch')) {
	    history.finished <- as.character(Sys.time())
		#history.command <- match.call()
	    history.command <- capture.output(print(match.call(lumiR.batch)))  

		## replace with the real file name
		if (length(fileList) > 1) {
			fileList <- paste('c("', paste(fileList, collapse='","', sep=''), '")', sep='')
		} else {
			fileList <- paste('"', fileList, '"', sep='')
		}
		if (length(grep(',', history.command)) > 0) {
			history.command <- sub('\\(.+,', paste('(', fileList, ',', sep=''), history.command)
		} else {
			history.command <- sub('\\(.+\\)', paste('(', fileList, ')', sep=''), history.command)
		}

		lumiVersion <- packageDescription('lumi')$Version
		x.lumi@history<- data.frame(submitted=history.submitted, 
				finished=history.finished, command=history.command, lumiVersion=lumiVersion)		
	}
    
	## initialize the QC slot in the LumiBatch object
	if (QC) x.lumi <- lumiQ(x.lumi, detectionTh=detectionTh, verbose=verbose)

	setwd(oldDir)
	return(x.lumi)
}