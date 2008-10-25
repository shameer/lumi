bgAdjust <- function(lumiBatch, probs=0.5, ...) {
	if (!is(lumiBatch, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	## by subtract an offset, which is estimated based on the quantile of the control probes
	if (min(exprs(lumiBatch)) <= 1) {
		cat('The data has already been background adjusted!\n')
		return(lumiBatch)
	}
	control <- lumiBatch@controlData
	if (is.null(control) || nrow(control) == 0) {
		cat('There is no control probe information in the LumiBatch object!\n No background adjustment will be performed.\n')
		return(lumiBatch)
	}
	colName <- colnames(control)
	sampleName <- sampleNames(lumiBatch)
	if (!all(sampleName %in% colName)) {
		sampleID <- pData(phenoData(lumiBatch))$sampleID
		if (!all(sampleID %in% colName)) {
			cat('Column names of controlData does not match with the LumiBatch object!\n No background adjustment will be performed.\n')
			return(lumiBatch)			
		} else {
			control <- control[, sampleID]
		}
	} else {
		control <- control[, sampleName]		
	}
	probeType <- rownames(control)
	if ('negative' %in% probeType) control <- control[probeType == 'negative',]
	quantile.ctrl <- apply(control, 2, quantile, probs=probs, ...)
	exprs(lumiBatch) <- exprs(lumiBatch) - matrix(rep(1, nrow(lumiBatch)), ncol=1) %*% quantile.ctrl
	return(lumiBatch)
}