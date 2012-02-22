`produceMethylationGEOSubmissionFile` <-
function(methyLumiM, methyLumiM.raw=NULL, lib.mapping=NULL, idType='Probe', sampleInfo=NULL, fileName='GEOSubmissionFile.txt', supplementaryRdata=TRUE, ...) {
	if (missing(methyLumiM)) stop('Please provide all required input parameters!\n')
	if (is(methyLumiM, "MethyLumiM")) expr.norm <- estimateBeta(methyLumiM)
	if (is.null(methyLumiM.raw)) {
		detect <- detection(methyLumiM)
		methyData <- methylated(methyLumiM)
		unmethyData <- unmethylated(methyLumiM)
		expr <- NULL
	} else {
		detect <- detection(methyLumiM.raw)
		methyData <- methylated(methyLumiM.raw)
		unmethyData <- unmethylated(methyLumiM.raw)
		expr <- estimateBeta(methyLumiM.raw)
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
	for (i in seq(sampleID)) {
		if (i == 1) {
			cat('^SAMPLE =', sampleTitle[i], '\n', sep='', file=fileName, append=FALSE)
		} else {
			cat('^SAMPLE =', sampleTitle[i], '\n', sep='', file=fileName, append=TRUE)			
		}
		sampleInfo.i <- paste('!', sampleInfoTitle[-1], ' = ', sampleInfo[i,-1], '\n', sep='', collapse='')
		sampleInfo.i <- gsub("'", "\\'", sampleInfo.i)
		cat(sampleInfo.i, file=fileName, append=TRUE, sep='')
		tableHead <- "ID_REF"
		cat("#ID_REF = Illumina ID\n", file=fileName, append=TRUE)
		if (!is.null(nuID)) {
			cat("#nuID = nucleotide universal IDentifier (nuID), convertible to and from probe sequence. See Bioconductor lumi package for more details.\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "nuID")
		}
		cat("#VALUE = Beta-value\n", file=fileName, append=TRUE)
		if (!is.null(expr)) cat("#RAW_VALUE = raw Beta-value\n", file=fileName, append=TRUE)
		tableHead <- c(tableHead, "VALUE")
		if (!is.null(expr)) tableHead <- c(tableHead, "RAW_VALUE")
		if (!is.null(methyData)) {
			cat("#METHYLATED = the intensities measured by methylated probes\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "METHYLATED")
		}
		if (!is.null(unmethyData)) {
			cat("#UNMETHYLATED = the intensities measured by unmethylated probes\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "UNMETHYLATED")
		}
		if (!is.null(detect)) {
			cat("#Detection_Pval = the detection p-value of the probe\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "Detection_Pval")
		}
		sampleTable.i <- probeId
		if (!is.null(nuID)) sampleTable.i <- cbind(sampleTable.i, nuID)
		sampleTable.i <- cbind(sampleTable.i, expr.norm[,sampleID[i]])
		if (!is.null(expr)) sampleTable.i <- cbind(sampleTable.i, expr[,sampleID[i]])
		if (!is.null(methyData)) sampleTable.i <- cbind(sampleTable.i, methyData[,sampleID[i]])
		if (!is.null(unmethyData)) sampleTable.i <- cbind(sampleTable.i, unmethyData[,sampleID[i]])
		if (!is.null(detect)) sampleTable.i <- cbind(sampleTable.i, detect[,sampleID[i]])
		sampleTable.i <- rbind(tableHead, sampleTable.i)
		cat("!sample_table_begin\n", file=fileName, append=TRUE)
		write.table(sampleTable.i, sep='\t', quote=FALSE, file=fileName, append=TRUE, col.names=FALSE, row.names=FALSE)
		cat("!sample_table_end\n", file=fileName, append=TRUE)
	}
	
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

