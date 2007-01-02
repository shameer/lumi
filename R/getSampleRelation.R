`getSampleRelation` <-
function(x, selProbe=NULL, Th=0.1, standardize=TRUE, method=c('mds', 'cluster'), dimension=c(1,2), col=NULL, ifPlot=FALSE) {
	if (is(x, 'ExpressionSet')) {
		dataMatrix <- exprs(x)
	} else if (is.matrix(x)) {
		dataMatrix <- x
	} else {
		stop('The class of "x" should be matrix or LumiBatch!')
	}

	if (standardize) dataMatrix <- scale(dataMatrix)

	if (is.null(selProbe)) {
		## Filter the genes with most of the experiments "Absent"
		cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
		probeList <- rownames(dataMatrix)
		selProbe <- probeList[cv.gene > Th]
	}
	dd <- dist(t(dataMatrix[selProbe,]))
	if (ifPlot) {
		method <- match.arg(method)
		if (method == 'cluster') {
			hc = hclust(dd, 'ave')
			plot(hc, xlab='Sample', main=paste('Clusters of the samples based on', length(selProbe), 'genes with sd/mean >', Th))
		} else {
			## Multi-Dimension Scaling
			a1 <- cmdscale(dd, k=max(dimension))
			if (is.null(col)) {
				color <- 1
			} else {
				if (!is.numeric(col)) {
					allColor <- colors()
					if (!all(is.element(col, allColor))) {
						color <- as.numeric(factor(col))						
					} else {
						color <- col
					}
				} else {
					color <- col
				}
			}
			plot(a1[,dimension[1]],a1[,dimension[2]], type='n', xlab=paste('Principle component', dimension[1]),ylab=paste('Principle component', dimension[2]))
			text(a1[,dimension[1]],a1[,dimension[2]], col=color, labels=colnames(dataMatrix), cex=1)
			title(paste('Sample relations based on', length(selProbe), 'genes with sd/mean >', Th))
		}
	} else {
		attr(dd, 'geneNum') <- length(selProbe)
		attr(dd, 'threshold') <- Th
		return(dd)
	}
}

