`getControlData` <-
function(x) 
{
	if (is.character(x)) {
		allControlInfo <- lumiR.batch(x, lib=NULL, checkDupId=FALSE)
		controlData <- as.data.frame(exprs(allControlInfo))
		controlType <- pData(featureData(allControlInfo))$TargetID
		ProbeID <- pData(featureData(allControlInfo))$ProbeID
		controlData <- data.frame(controlType=controlType, ProbeID=ProbeID, controlData)
	} else if (is(x, 'LumiBatch')) {
		controlData <- x@controlData
	} else {
		stop('Input data should be a control data file or a LumiBatch object!')
	}
	return(controlData)
}

