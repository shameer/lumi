`plotSampleRelation` <-
function(x, subset=NULL, cv.Th=0.1, standardize=TRUE, method=c('cluster', 'mds'), dimension=c(1,2), color=NULL, main=NULL, pch=NULL, addLegend=TRUE, ...) {
	if (is(x, 'ExpressionSet')) {
		dataMatrix <- exprs(x)
	} else if (is.matrix(x)) {
		dataMatrix <- x
	} else {
		stop('The class of "x" should be matrix or LumiBatch!')
	}
	
	## Standardize each sample
	if (standardize) dataMatrix <- scale(dataMatrix)

	if (is.null(subset)) {
		## Filter the genes with most of the experiments "Absent"
		probeList <- rownames(dataMatrix)
		if (is.null(probeList)) probeList <- 1:nrow(dataMatrix)
		if (cv.Th > 0) {
			cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
			subset <- probeList[abs(cv.gene) > cv.Th]
			if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'genes with sd/mean >', cv.Th)
		} else {
			subset <- probeList
			if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'genes')
		}
	} else {
		if (length(subset) == 1 && is.numeric(subset)) {
			subset <- sample(1:nrow(dataMatrix), min(subset, nrow(dataMatrix)))
		}
		if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'selected genes')
	}

	dd <- dist(t(dataMatrix[subset,]))
	method <- match.arg(method)
	if (method == 'cluster') {
		hc = hclust(dd, 'ave')
		plot(hc, xlab='Sample', main=main, ...)
		attr(hc, 'geneNum') <- length(subset)
		attr(hc, 'threshold') <- cv.Th
		return(invisible(hc))	
	} else {
		## Multi-Dimension Scaling
		mds.result <- cmdscale(dd, k=max(dimension), eig=TRUE)
		ppoints <- mds.result$points
		eig <- mds.result$eig
		percent <- round(eig/sum(eig) * 100, 1)
		
		colorLegend <- NULL
		if (is.null(color)) {
			color <- 1
		} else {
			if (!is.numeric(color)) {
				allColor <- colors()
				if (!all(is.element(color, allColor))) {
					colorLegend <- unique(color)
					color <- as.numeric(factor(color, levels=colorLegend))
				} 
			}
		}
		if (missing(pch)) {
			plot(ppoints[,dimension[1]], ppoints[,dimension[2]], type='n', 
				xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),
				ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), 
				main=main, ...)
			text(ppoints[,dimension[1]], ppoints[,dimension[2]], col=color, labels=colnames(dataMatrix), cex=1)
		} else {
			plot(ppoints[,dimension[1]], ppoints[,dimension[2]],  
				xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),
				ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), 
				main=main, col=color, pch=pch, ...)
		}
		attr(ppoints, 'geneNum') <- length(subset)
		attr(ppoints, 'threshold') <- cv.Th
		
		## add legend if color is a factor
		if (!is.null(colorLegend) && addLegend) {
			if (!missing(pch)) {
				legend('topleft', legend=colorLegend, col=unique(color), pch=unique(pch))
			} else {
				legend('topleft', legend=colorLegend, text.col=1:length(colorLegend))
			}
		}
		
		return(invisible(ppoints))	
	}
}

