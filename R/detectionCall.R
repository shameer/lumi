`detectionCall` <-
function(x.lumi, Th=0.99) {

	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	detect <- detection(x.lumi)
	if (!is.null(detect)) {
		AP <- colSums(detect >= Th) / nrow(detect)
		attr(AP, 'threshold') <- Th
	} else {
		AP <- NULL
	}
	return(AP)
}

