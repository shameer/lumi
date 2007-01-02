`estimateLumiCV` <-
function(x.lumi, ifPlot=FALSE) {

	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	cv <- se.exprs(x.lumi) / exprs(x.lumi)
	if (ifPlot) {
		plotDensity(cv, xlab='cv')
	}
	return(cv)
}

