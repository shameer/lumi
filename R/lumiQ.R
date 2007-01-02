`lumiQ` <-
function(x.lumi, logMode=TRUE, sampleRelation.param=list(), detectOutLier.param=list()) {
	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	history.submitted <- as.character(Sys.time())

	exprs <- exprs(x.lumi)
	if (any(is.na(exprs))) {
		naInd <- apply(exprs, 1, function(x) any(is.na(x)))
		exprs <- exprs[!naInd,]
	}

	if (logMode & (max(exprs, na.rm=TRUE) > 50)) {
		exprs <- log2(exprs)
	} 
	sampleName <- colnames(exprs)

	## mean, variance, 'cv', 'AP', 'density', 'correlation', 'sample relation'
	mm <- colMeans(exprs)
	std <- apply(exprs, 2, sd)
	## Correlation between sample
	sampleCor <- cor(exprs)
	## calculate correlation of variance
	cv <- estimateLumiCV(x.lumi, ifPlot=FALSE)
	## AP calls
	detectionRate <- detectionCall(x.lumi)

	## relation of samples
	#sampleRelation <- getSampleRelation(exprs, ifPlot=FALSE, ...)
	sampleRelation <- do.call("getSampleRelation", c(alist(exprs, ifPlot=FALSE), sampleRelation.param))

	## detect outlier
	#outlier <- detectOutlier(exprs, ifPlot=FALSE, ...)
	outlier <- do.call("detectOutlier", c(alist(exprs, ifPlot=FALSE), detectOutLier.param))

	## record history
    history.finished <- as.character(Sys.time())
    history.command <- capture.output(print(match.call(lumiQ)))  
	history<- rbind(x.lumi@history, c(history.submitted, history.finished, history.command))

	## create a Lumi.QC object
	names(mm) <- names(std) <- names(detectionRate) <- names(outlier) <- sampleName
	lumiQC <- new('LumiQC', exprs=exprs(x.lumi), cv=cv, mean=mm, std=std, sampleCor=sampleCor,  
	detectionRate=detectionRate, sampleRelation=list(sampleRelation), outlier=outlier, history=history)
	return(lumiQC)
}

