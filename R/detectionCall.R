`detectionCall` <-
function(x.lumi, Th = 0.01, type=c('probe', 'sample')) {
	type <- match.arg(type)
	if (is(x.lumi, 'matrix')) {
		detect <- x.lumi
	} else {
		if (!is(x.lumi, 'LumiBatch')) 
			stop('The object should be class "LumiBatch"!')
		detect <- detection(x.lumi)
	} 
	if (Th > 0.8) {
		detect <- 1 - detect
		Th <- 1 - Th
	}
	if (!is.null(detect)) {
		if (type == 'sample') AP <- colSums(detect<= Th)
		if (type == 'probe') AP <- rowSums(detect<= Th)
		attr(AP, 'threshold') <- Th
	} else {
		AP <- rep(NA, ncol(x.lumi))
	}
	return(AP)
}

