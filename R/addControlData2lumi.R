`addControlData2lumi` <-
function(controlData, x.lumi) 
{
	if (missing(x.lumi) || missing(controlData)) stop('Both controlData and x.lumi are required!')
	if (is.character(controlData)) {
		controlData <- getControlData(controlData, type='data.frame')
	}
	if (is.matrix(controlData)) controlData <- as.data.frame(controlData)
	if (is(controlData, 'data.frame')) {
		## match the column names of controlData and LumiBatch object
		sampleID <- as.character(pData(phenoData(x.lumi))$sampleID)
		if (length(sampleID) == 0) sampleID <- sampleNames(x.lumi)
		controlSampleID <- names(controlData)
		if ('TargetID' %in% controlSampleID) {
			names(controlData)[controlSampleID == 'TargetID'] <- 'controlType'
		} else {
			names(controlData)[1] <- 'controlType'
		}
		probeId.pos <- grep('Probe.?ID', controlSampleID, ignore.case=TRUE)
		if (length(probeId.pos) > 0) {
			names(controlData)[probeId.pos] <- 'ProbeID'
		} else {
			controlData$ProbeID <- NA
		}
		retrieveColName <- c('controlType', 'ProbeID')
		
		if (all(sampleID %in% controlSampleID)) {
			x.lumi@controlData <- controlData[, c(retrieveColName, sampleID)]
		} else {
			sampleIDInfo <- strsplit(sampleID, split="_")
			newID <- NULL
			temp <- lapply(sampleIDInfo, function(x) {
				newID <<- c(newID, paste(x[1:2], collapse="_"))
			})
			if (all(newID %in% controlSampleID)) {
				controlData(x.lumi) <- controlData[, c(retrieveColName, newID)]
			} else {
				stop('SampleID does not match up between controlData and x.lumi!')				
			}
		}
		names(controlData(x.lumi)) <- c(retrieveColName, sampleNames(x.lumi))		
	} else {
		stop('Input data type is not supported!')
	}
	return(x.lumi)
}

