`lumiExpresso` <- 
function (lumiBatch, bg.correct = TRUE, bgcorrect.param = list(method='bgAdjust'), variance.stabilize = TRUE, 
	varianceStabilize.param = list(), normalize=TRUE, normalize.param = list(), 
	QC.evaluation = TRUE, QC.param = list(), verbose = TRUE) 
{
	if (verbose) {
		if (bg.correct) {
			bgMethod <- ifelse(is.null(bgcorrect.param$method), 'bgAdjust', bgcorrect.param$method)
			cat("Background Correction:", bgMethod, "\n")
		}
		if (variance.stabilize) {
			vstMethod <- ifelse(is.null(varianceStabilize.param$method), 'vst', varianceStabilize.param$method)
			cat("Variance Stabilizing Transform method:", vstMethod, "\n")
		}
		if (normalize) {
			normMethod <- ifelse(is.null(normalize.param$method), 'quantile', normalize.param$method)
			cat("Normalization method:", normMethod, "\n")
		}
		cat('\n')
	}
	if (bg.correct) {
		if (verbose) cat("\nBackground correction ...\n")
		lumiBatch <- do.call("lumiB", c(alist(lumiBatch), bgcorrect.param))
		if (verbose) cat("done.\n")
	}
	if (variance.stabilize) {
		if (verbose) cat("\nVariance stabilizing ...\n")
		lumiBatch <- do.call("lumiT", c(alist(lumiBatch), varianceStabilize.param))
		if (verbose) cat("done.\n")
	}
	if (normalize) {
		if (verbose) cat("\nNormalizing ...\n")
		lumiBatch <- do.call("lumiN", c(alist(lumiBatch), normalize.param))
		if (verbose) cat("done.\n")
	}
	if (QC.evaluation) {
		if (verbose) cat("\nQuality control after preprocessing ...\n")
		lumiBatch <- do.call("lumiQ", c(alist(lumiBatch), QC.param))
		if (verbose) cat("done.\n")
	}
	return(lumiBatch)
}