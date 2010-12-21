`detectionCall` <-
function(x.lumi, Th = 0.01, type=c('probe', 'sample', 'matrix')) {
	type <- match.arg(type)
	if (Th > 0.5)  Th <- 1 - Th
	if (is(x.lumi, 'matrix')) {
		detect <- x.lumi
	} else {
		if (!assayDataValidMembers(assayData(x.lumi), "detection")) 
			stop('The object should be class "eSet" inherited classes and include "detection" element in the assayData!')
		detect <- assayDataElement(x.lumi, 'detection')
		if (is.null(detect)) {
			warning('No detection slot found!')
			return(NULL)
		}
		
		## check the detection is p-values or not
		if (class(x.lumi) == "MethyLumiM") {
			expr <- estimateIntensity(x.lumi, returnType='matrix')
		} else {
			expr <- exprs(x.lumi)
		}
		# low <- mean(expr[detect[,1] > 0.9,1])
		# high <- mean(expr[detect[,1] < 0.1,1])
		low <- expr[which.max(detect[,1]), 1]
		high <- expr[which.min(detect[,1]), 1]
		if (low > high) detect <- 1 - detect
	}
	
	if (!is.null(detect)) {
		if (type == 'sample') AP <- colSums(detect<= Th)
		if (type == 'probe') AP <- rowSums(detect<= Th)
		if (type == 'matrix') {
			AP <- detect
			AP[detect <= Th] <- 'P'
			AP[detect > Th] <- 'A'
		}
		attr(AP, 'threshold') <- Th
	} else {
		AP <- rep(NA, ncol(x.lumi))
	}
	return(AP)
}

