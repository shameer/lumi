`detectionCall` <-
function(x.lumi, Th = 0.01) {
	if (is(x.lumi, 'matrix')) {
		detect <- x.lumi
	} else {
		if (!is(x.lumi, 'LumiBatch')) 
			stop('The object should be class "LumiBatch"!')
		detect <- detection(x.lumi)
	} 
	if (!is.null(detect)) {
		AP <- colSums(detect<= Th) / nrow(detect)
		attr(AP, 'threshold') <- Th
	} else {
		AP <- rep(NA, ncol(x.lumi))
	}
	return(AP)
}

