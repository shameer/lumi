`plotDensity` <-
function(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l",  
		index.highlight=NULL, color.highlight=2, symmetry=NULL, addLegend=TRUE, ...) 
{
	if (is(x, 'ExpressionSet')) {
	    index <- round(seq(1, nrow(x), len=5000))
	    mat <- exprs(x)[index, , drop = FALSE]
	} else if (is.numeric(x)) {
		mat <- as.matrix(x)
	} else {
		stop('Un-supported class of x!')
	}
	
    if (logMode & (max(mat, na.rm=TRUE) > 50)) {
        mat <- log2(mat)
        if (is.null(xlab)) 
            xlab <- "log2 intensity"
    } else if (is.null(xlab)) 
        xlab <- "intensity"

	if (!is.null(symmetry)) {
		x.range <- range(mat)
		if (symmetry > x.range[1] & symmetry < x.range[2]) {
			warning('symmetry point should not be within the range of x!')
			symmetry <- NULL
		} else {
			mat <- rbind(mat, 2*symmetry - mat)
		}
	}
	x.density <- apply(mat, 2, density)
    all.x <- do.call("cbind", lapply(x.density, function(x) x$x))
    all.y <- do.call("cbind", lapply(x.density, function(x) x$y))

	if (!is.null(symmetry)) {
		nr <- nrow(all.x)
		if (all.x[1,1] >= x.range[1] & all.x[1,1] <= x.range[2]) {
			all.x <- all.x[1:round(nr/2),]
			all.y <- all.y[1:round(nr/2),]
		} else {
			all.x <- all.x[round(nr/2):nr,]
			all.y <- all.y[round(nr/2):nr,]
		}
	}
    matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=1:ncol(all.x), 
        lty=1:ncol(all.x), ...)
	if (!is.null(index.highlight)) {
		if (index.highlight > ncol(all.x) | index.highlight < 1) {
			warning('Highlight index out of range!')
			index.highlight <- 1
		}
		lines(all.x[,index.highlight], all.y[,index.highlight], col=color.highlight, lwd=2, lty=1)
	}
	## add legend
	if (addLegend) {
		labels <- colnames(mat)
		if (is.null(labels)) labels <- as.character(1:ncol(mat))

		col <- 1:ncol(all.x)
		lwd <- rep(1, ncol(all.x))
		lty <- col
		if (!is.null(index.highlight)) {
			col[index.highlight] <- color.highlight
			lwd[index.highlight] <- 2
			lty[index.highlight] <- 1
		}
		x.pos <- (max(all.x) - min(all.x)) * 2/3 + min(all.x)
		y.pos <- max(all.y)
		legend(x.pos, y.pos, legend=labels, col=col, lwd=lwd, lty=lty)
	}
	
    invisible(list(all.x = all.x, all.y = all.y))
}

