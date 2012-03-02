`produceMethylationGEOSubmissionFile` <-
function(methyLumiM, methyLumiM.raw=NULL, lib.mapping=NULL, idType='Probe', sampleInfo=NULL, fileName='GEOSubmissionFile.txt', supplementaryRdata=FALSE, ...) {
	if (missing(methyLumiM)) stop('Please provide all required input parameters!\n')
	options(stringsAsFactors=FALSE)
	if (is(methyLumiM, "MethyLumiM")) expr.norm <- signif(estimateBeta(methyLumiM, returnType='matrix'),5)
	if (is.null(methyLumiM.raw)) {
		detect <- signif(detection(methyLumiM),5)
		methyData <- signif(methylated(methyLumiM),5)
		unmethyData <- signif(unmethylated(methyLumiM),5)
		expr <- NULL
	} else {
		detect <- signif(detection(methyLumiM.raw),5)
		methyData <- signif(methylated(methyLumiM.raw),5)
		unmethyData <- signif(unmethylated(methyLumiM.raw),5)
		expr <- signif(estimateBeta(methyLumiM.raw, returnType='matrix'),5)
	}
		
	if (is.null(sampleInfo)) {
		sampleInfo <- produceGEOSampleInfoTemplate(methyLumiM, lib.mapping=lib.mapping, fileName=NULL)
	} else if (length(sampleInfo) == 1 && is.character(sampleInfo)) {
		sampleInfo <- read.table(sampleInfo, sep='\t', colClasses='character', skip=1, header=TRUE, strip.white=TRUE, quote='')
	} else if (is.null(nrow(sampleInfo))) {
		stop('Please provide correct sample information (a data.frame, matrix, or sampleInfo file)!\n')
	}
	sampleInfoTitle <- colnames(sampleInfo)
	if (any(sapply(sampleInfo[,-1, drop=F], nchar) == 0)) stop('No blank fields are allowed in the sampleInfo table!\nYou can check some example submissions, like GSM296418, at the GEO website.\n')
	if (supplementaryRdata) sampleInfo[, "Sample_supplementary_file"] <- 'supplementaryData.Rdata'
	nuID <- featureNames(methyLumiM)
	if (!supplementaryRdata) {
	  rm(methyLumiM.raw, methyLumiM)
	  gc()
	}
	
	probeId <- nuID
	if (length(which(is.nuID(sample(nuID, 100)))) < 20) {
		nuID <- NULL
	} else {
		if (!is.null(lib.mapping)) {
			probeId <- nuID2IlluminaID(nuID, lib.mapping=lib.mapping, idType=idType, ...)
		} else {
			nuID <- NULL
		}
	}
	
	sampleID <- sampleInfo[, "sampleID"]
	sampleTitle <- sampleInfo[,'Sample_title']
	# outputFile <- file(fileName, "w")
	outputFile <- fileName
	for (i in seq(sampleID)) {
	    cat('Processing sample', i, '\n')
		if (i == 1) {
			cat('^SAMPLE =', sampleTitle[i], '\n', sep='', file=outputFile, append=FALSE)
		} else {
			cat('^SAMPLE =', sampleTitle[i], '\n', sep='', file=outputFile, append=TRUE)			
		}
		sampleInfo.i <- paste('!', sampleInfoTitle[-1], ' = ', sampleInfo[i,-1], '\n', sep='', collapse='')
		sampleInfo.i <- gsub("'", "\\'", sampleInfo.i)
		cat(sampleInfo.i, file=outputFile, sep='', append=TRUE)
		tableHead <- "ID_REF"
		cat("#ID_REF = Illumina ID\n", file=outputFile, append=TRUE)
		if (!is.null(nuID)) {
			cat("#nuID = nucleotide universal IDentifier (nuID), convertible to and from probe sequence. See Bioconductor lumi package for more details.\n", file=outputFile, append=TRUE)
			tableHead <- c(tableHead, "nuID")
		}
		cat("#VALUE = Beta-value\n", file=outputFile, append=TRUE)
		if (!is.null(expr)) cat("#RAW_VALUE = raw Beta-value\n", file=outputFile, append=TRUE)
		tableHead <- c(tableHead, "VALUE")
		if (!is.null(expr)) tableHead <- c(tableHead, "RAW_VALUE")
		if (!is.null(methyData)) {
			cat("#METHYLATED = the intensities measured by methylated probes\n", file=outputFile, append=TRUE)
			tableHead <- c(tableHead, "METHYLATED")
		}
		if (!is.null(unmethyData)) {
			cat("#UNMETHYLATED = the intensities measured by unmethylated probes\n", file=outputFile, append=TRUE)
			tableHead <- c(tableHead, "UNMETHYLATED")
		}
		if (!is.null(detect)) {
			cat("#Detection_Pval = the detection p-value of the probe\n", file=outputFile, append=TRUE)
			tableHead <- c(tableHead, "Detection_Pval")
		}
		sampleTable.i <- probeId
		if (!is.null(nuID)) sampleTable.i <- cbind(sampleTable.i, nuID)
		sampleTable.i <- cbind(sampleTable.i, expr.norm[,sampleID[i], drop=FALSE])
	
		if (!is.null(expr)) sampleTable.i <- cbind(sampleTable.i, expr[,sampleID[i], drop=FALSE])
		if (!is.null(methyData)) sampleTable.i <- cbind(sampleTable.i, methyData[,sampleID[i], drop=FALSE])
		if (!is.null(unmethyData)) sampleTable.i <- cbind(sampleTable.i, unmethyData[,sampleID[i], drop=FALSE])
		if (!is.null(detect)) sampleTable.i <- cbind(sampleTable.i, detect[,sampleID[i], drop=FALSE])
		sampleTable.i <- rbind(tableHead, sampleTable.i)
		cat("!sample_table_begin\n", file=outputFile, append=TRUE)
		write.table(sampleTable.i, sep='\t', quote=FALSE, file=outputFile, append=TRUE, col.names=FALSE, row.names=FALSE)
		cat("!sample_table_end\n", file=outputFile, append=TRUE)
	}
	#close(outputFile)
	
	if (supplementaryRdata) {
		methyLumiM <- methyLumiM[,sampleID]
		if (!is.null(methyLumiM.raw)) {
			methyLumiM.raw <- methyLumiM.raw[,sampleID]
			save(methyLumiM, methyLumiM.raw, sampleInfo, file='supplementaryData.Rdata')
		} else {
			save(methyLumiM, sampleInfo, file='supplementaryData.Rdata')
		}
	}
}

