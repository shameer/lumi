`plotSampleRelation` <-
function(x, selProbe=NULL, cv.Th=0.1, standardize=TRUE, method=c('cluster', 'mds'), dimension=c(1,2), color=NULL, main=NULL, ...) {
	if (is(x, 'ExpressionSet')) {
		dataMatrix <- exprs(x)
	} else if (is.matrix(x)) {
		dataMatrix <- x
	} else {
		stop('The class of "x" should be matrix or LumiBatch!')
	}
	
	## Standardize each sample
	if (standardize) dataMatrix <- scale(dataMatrix)

	if (is.null(selProbe)) {
		## Filter the genes with most of the experiments "Absent"
		probeList <- rownames(dataMatrix)
		if (cv.Th > 0) {
			cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
			selProbe <- probeList[abs(cv.gene) > cv.Th]
			if (is.null(main)) main <- paste('Sample relations based on', length(selProbe), 'genes with sd/mean >', cv.Th)
		} else {
			selProbe <- probeList
			if (is.null(main)) main <- paste('Sample relations based on', length(selProbe), 'genes')
		}
	} else {
		if (is.null(main)) main <- paste('Sample relations based on', length(selProbe), 'selected genes')
	}

	dd <- dist(t(dataMatrix[selProbe,]))
	method <- match.arg(method)
	if (method == 'cluster') {
		hc = hclust(dd, 'ave')
		plot(hc, xlab='Sample', main=main, ...)
		attr(hc, 'geneNum') <- length(selProbe)
		attr(hc, 'threshold') <- cv.Th
		return(invisible(hc))	
	} else {
		## Multi-Dimension Scaling
		mds.result <- cmdscale(dd, k=max(dimension), eig=TRUE)
		ppoints <- mds.result$points
		eig <- mds.result$eig
		percent <- round(eig/sum(eig) * 100, 1)

		if (is.null(color)) {
			color <- 1
		} else {
			if (!is.numeric(color)) {
				allColor <- colors()
				if (!all(is.element(color, allColor))) {
					color <- as.numeric(factor(color, levels=unique(color)))
				} 
			}
		}
		plot(ppoints[,dimension[1]], ppoints[,dimension[2]], type='n', xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), main=main, ...)
		text(ppoints[,dimension[1]], ppoints[,dimension[2]], col=color, labels=colnames(dataMatrix), cex=1)
		attr(ppoints, 'geneNum') <- length(selProbe)
		attr(ppoints, 'threshold') <- cv.Th
		return(invisible(ppoints))	
	}
}

