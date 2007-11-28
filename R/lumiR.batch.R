lumiR.batch <- function(fileList, lib = NULL, transform = c('none', 'vst', 'log2', 'cubicRoot'), sampleInfoFile = NULL, ...) {

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
		if (length(fileList) == 0) stop('No data files were found!')
	}

	print('Inputting the data ...')
	for (i in 1:length(fileList)) {
		file.i <- fileList[i]
		x.lumi.i <- lumiR(file.i, parseColumnName=FALSE, ...)
		# x.lumi.i <- lumiR(file.i, parseColumnName=FALSE)
		if (transform != 'none') {
			x.lumi.i <- lumiT(x.lumi.i, method=transform, simpleOutput=TRUE)
		}
		if (i == 1) {
			x.lumi <- x.lumi.i
		} else {
			x.lumi <- combine(x.lumi, x.lumi.i)
		}
	}
	if (!is.null(lib)) {
		print('Adding nuID to the data ...')
		x.lumi <- addNuId2lumi(x.lumi, lib=lib)
	}

	if (!is.null(sampleInfoFile)) {

		if (is.character(sampleInfoFile) || class(sampleInfoFile)[1] == 'file') {
			if (file.exists(sampleInfoFile)) {
				sampleInfo <- read.table(sampleInfoFile, head=TRUE, sep='\t', colClasses='character', comment='')
			} else {
				warning('The provided sampleInfoFile does not exist!')
				setwd(oldDir)
				return(x.lumi)
			}
		} else if (is.data.frame(sampleInfoFile)) {
			sampleInfo <- sampleInfoFile
		}

		ID <- sampleInfo$ID
		if (is.null(ID)) {
			warning('The ID column in sampleInfoFile is required!')
			setwd(oldDir)
			return(x.lumi)
		} 
		rownames(sampleInfo) <- ID

		sampleName <- sampleNames(x.lumi)
		if (!all(sampleName %in% ID)) {
			warning('Some sample informatin is not provided!')
		} 
		ID <- intersect(ID, sampleName)
		sampleNames(x.lumi) <- sampleInfo[sampleName, 'Label']
		ind <- 1:ncol(x.lumi)
		names(ind) <- sampleName
		x.lumi <- x.lumi[, ind[ID]]

		pData <- sampleInfo[ID,]
		label <- sampleInfo[ID, 'Label']
		if (!is.null(label)) {
			label <- ID
			rownames(pData) <- label
			sampleNames(x.lumi) <- label
		}

		pdata <- new("AnnotatedDataFrame", data=pData)
		phenoData(x.lumi) <- pdata
	} 
	setwd(oldDir)
	return(x.lumi)
}