`lumiQ` <-
function(x.lumi, logMode=TRUE, detectionTh=0.01) {
	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	history.submitted <- as.character(Sys.time())

	exprs <- exprs(x.lumi)
	if (any(is.na(exprs))) {
		naInd <- apply(exprs, 1, function(x) any(is.na(x)))
		exprs <- exprs[!naInd,]
	}

	if (logMode & (max(exprs, na.rm=TRUE) > 50)) {
		if (min(exprs) < 0) {
			warning('Negative values found in the expression values!')
			exprs <- exprs + abs(min(exprs)) + 1
		}
		exprs <- log2(exprs)
	} 
	sampleName <- colnames(exprs)

	## mean, variance, 'cv', 'AP', 'density', 'correlation', 'sample relation'
	mm <- colMeans(exprs)
	std <- apply(exprs, 2, sd)
	## AP calls
	detectionRate <- detectionCall(x.lumi, Th=detectionTh)

	## detect outlier
	center <- rowMeans(exprs)
	profile <- cbind(center, exprs)
	colnames(profile) <- c('Center', colnames(exprs))
	distCenter <- as.matrix(dist(t(profile), method="euclidean"))

	sampleSummary <- rbind(mm, std, detectionRate, distCenter[2:nrow(distCenter),1])
	rownames(sampleSummary) <- c('mean', 'standard deviation', paste('detection rate(', detectionTh, ')', sep=''), 'distance to sample mean')
	sampleSummary <- signif(sampleSummary, 4)

	## record history
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiQ)))

	if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
	lumiVersion <- packageDescription('lumi')$Version
	x.lumi@history <- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))

	## create a QC slot
	QC <- x.lumi@QC
	QC$sampleSummary <- sampleSummary
	QC$history <- x.lumi@history
	
	x.lumi@QC <- QC
	return(x.lumi)
}

