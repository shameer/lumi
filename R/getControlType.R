`getControlType` <-
function(controlData) 
{
	if (is(controlData, 'LumiBatch')) {
		controlData <- controlData@controlData
	} 
	if (is(controlData, 'data.frame')) {
		return(unique(controlData$controlType))
	} else {
		return(NA)
	}
}

