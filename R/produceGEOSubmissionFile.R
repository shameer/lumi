`produceGEOSubmissionFile` <-
function(lumiNormalized, lumiRaw, lib.mapping, sampleInfo=NULL, fileName='GEOSubmissionFile.txt', supplementaryRdata=TRUE, ...) {
	if (missing(lumiNormalized) || missing(lumiRaw) || missing(lib.mapping)) stop('Please provide all required input parameters!\n')
	expr.norm <- exprs(lumiNormalized)
	detect <- detection(lumiRaw)
	se.expr <- se.exprs(lumiRaw)
	expr <- exprs(lumiRaw)
	beadNum <- beadNum(lumiRaw)
	if (is.null(sampleInfo)) {
		sampleInfo <- produceGEOSampleInfoTemplate(lumiNormalized, lib.mapping=lib.mapping, fileName=NULL)
	} else if (length(sampleInfo) == 1 && is.character(sampleInfo)) {
		sampleInfo <- read.table(sampleInfo, sep='\t', colClasses='character', skip=1, head=TRUE, strip.white=TRUE, quote='')
	} else if (is.null(nrow(sampleInfo))) {
		stop('Please provide correct sample information (a data.frame, matrix, or sampleInfo file)!\n')
	}
	sampleInfoTitle <- colnames(sampleInfo)
	if (any(sapply(sampleInfo[-1,-1], nchar) == 0)) stop('No blank fields are allowed in the sampleInfo table!\nYou can check some example submissions, like GSM296418, at the GEO website.\n')
	if (supplementaryRdata) sampleInfo[, "Sample_supplementary_file"] <- 'supplementaryData.Rdata'
	nuID <- featureNames(lumiNormalized)
	if (!all(is.nuID(sample(nuID, 100)))) {
		probeId <- nuID
		nuID <- NULL
	} else {
		probeId <- nuID2probeID(nuID, lib=lib.mapping, ...)		
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
		cat("#ID_REF = \n", file=fileName, append=TRUE)
		if (!is.null(nuID)) {
			cat("#nuID = nucleotide universal IDentifier (nuID), convertible to and from probe sequence. See Bioconductor lumi package for more details.\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "nuID")
		}
		cat("#VALUE = normalized signal intensity\n", file=fileName, append=TRUE)
		cat("#RAW_VALUE = raw signal intensity\n", file=fileName, append=TRUE)
		tableHead <- c(tableHead, "VALUE")
		tableHead <- c(tableHead, "RAW_VALUE")
		if (!is.null(se.expr)) {
			cat("#BEAD_STDERR = the standard error of the probe measurements\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "BEAD_STDERR")
		}
		if (!is.null(detect)) {
			cat("#Detection_Pval = the detection p-value of the probe\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "Detection_Pval")
		}
		if (!is.null(beadNum)) {
			cat("#Avg_NBEADS = Number of beads for the probe\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "Avg_NBEADS")
		}
		sampleTable.i <- probeId
		if (!is.null(nuID)) sampleTable.i <- cbind(sampleTable.i, nuID)
		sampleTable.i <- cbind(sampleTable.i, expr.norm[,sampleID[i]], expr[,sampleID[i]])
		if (!is.null(se.expr)) sampleTable.i <- cbind(sampleTable.i, se.expr[,sampleID[i]])
		if (!is.null(detect)) sampleTable.i <- cbind(sampleTable.i, detect[,sampleID[i]])
		if (!is.null(beadNum)) sampleTable.i <- cbind(sampleTable.i, beadNum[,sampleID[i]])
		sampleTable.i <- rbind(tableHead, sampleTable.i)
		cat("!sample_table_begin\n", file=fileName, append=TRUE)
		write.table(sampleTable.i, sep='\t', quote=FALSE, file=fileName, append=TRUE, col.names=FALSE, row.names=FALSE)
		cat("!sample_table_end\n", file=fileName, append=TRUE)
	}
	
	if (supplementaryRdata) {
		lumiNormalized <- lumiNormalized[,sampleID]
		lumiRaw <- lumiRaw[,sampleID]
		save(lumiNormalized, lumiRaw, sampleInfo, file='supplementaryData.Rdata')
	}
}

