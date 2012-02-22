`getControlData` <-
function(x, type=c('data.frame', 'LumiBatch'), ...) 
{
	type <- match.arg(type)
	if (is.character(x)) {
		if (type == 'data.frame') {
			if ('columnNameGrepPattern' %in% names(list(...))) {
				allControlInfo <- lumiR.batch(x, lib.mapping=NULL, checkDupId=FALSE, convertNuID=FALSE, QC=FALSE, ...)
			} else {
				allControlInfo <- lumiR.batch(x, lib.mapping=NULL, checkDupId=FALSE, convertNuID=FALSE, columnNameGrepPattern = list(exprs='AVG_SIGNAL', se.exprs=NA, detection=NA, beadNum=NA), QC=FALSE, ...)				
			}
		} else {
			allControlInfo <- lumiR.batch(x, lib.mapping=NULL, checkDupId=FALSE, convertNuID=FALSE, QC=FALSE, ...)
			return(allControlInfo)
		}
		x <- allControlInfo
	} 

	if (is(x, 'ExpressionSet')) {
		if (type == 'LumiBatch') {
			return(x)
		} else {
			if ('controlData' %in% slotNames(x)) {
				if (nrow(x@controlData) > 0) {
					controlData <- x@controlData
					return(controlData)
				}
			} 
			controlData <- as.data.frame(exprs(x))
			controlType <- pData(featureData(allControlInfo))[,'TargetID']
			if (length(which(toupper(controlType) == 'NEGATIVE')) > 10) {
				ProbeID <- pData(featureData(allControlInfo))$ProbeID
				controlNames <- names(controlData)
				controlData <- data.frame(controlType=as.character(controlType), ProbeID=as.character(ProbeID), controlData)
				names(controlData) <- c('controlType', 'ProbeID', controlNames)
			} else {
				controlNames <- names(controlData)
				controlData <- data.frame(controlType=as.character(controlType), controlData)
				names(controlData) <- c('controlType', controlNames)
			}
		}		
	} else {
		stop('Input data should be a control data file or a LumiBatch object!')
	}
	return(controlData)
}

