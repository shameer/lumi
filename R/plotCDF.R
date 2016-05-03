## cdf plot in reverse direction (from high value to low value)
plotCDF <- function(x, reverse=TRUE, logMode=TRUE, xlab = NULL, ylab = "Cumulative density", col=1:dim(x)[2], 
	lwd=1, xlim=NULL, index.highlight=NULL, color.highlight=2, addLegend=TRUE,  main='',...) 
{
	if (is(x, 'ExpressionSet')) {
	    expr <- exprs(x)
	} else if (is.numeric(x)) {
		expr <- as.matrix(x)
	} else {
		stop('Un-supported class of x.')
	}
		
    if (logMode && (max(expr, na.rm=TRUE) > 50)) {
		# remove the negative values
		if (min(expr) < 0) {
			rMin <- rowMin(expr)
			expr <- expr[rMin > 0, , drop=FALSE]
		}
		expr <- log2(expr)
		if (is.null(xlab)) 
			xlab <- "log2 intensity"
    } else if (is.null(xlab)) 
        xlab <- "intensity"

	index <- 1:nrow(expr)
	expr <- expr[index,,drop=FALSE]
	
	if (reverse) {
		x.cdf <- apply(expr, 2, function(x) ecdf(-x))
	} else {
		x.cdf <- apply(expr, 2, function(x) ecdf(x))
	}
	usedCol <- usedLwd <- NULL
	for (i in 1:ncol(expr)) {
		col.i <- ifelse(i > length(col), col[1], col[i])
		lwd.i <- ifelse(i > length(lwd), lwd[1], lwd[i])
		if (i %in% index.highlight) {
			col.i <- color.highlight
			lwd.i <- 2
		}
		if (i == 1) {
			if (reverse) {
				x.range <- range(-expr)
			} else {
				x.range <- range(expr)
			}		
			plot(x.cdf[[i]], xlab=xlab, ylab = ylab, verticals = F, col=col.i,  lwd=lwd.i, xlim=x.range, axes=F, main=main, ...)
			tmp <- Axis(x=seq(x.range[1], x.range[2], len=100), side=1, labels=NA)
			if (reverse) {
				axis(side = 1, at = tmp, labels = -tmp)
			} else {
				axis(side = 1, at = tmp, labels = tmp)
			}
			Axis(x=seq(0,1,len=100), side=2)
			box()
		} else {
			plot(x.cdf[[i]], xlab=NULL, ylab = NULL, verticals = F, col=col.i, lwd=lwd.i, add = T, ...)
		}		
		usedCol <- c(usedCol, col.i)
		usedLwd <- c(usedLwd, lwd.i)
	}

	## add legend
	if (addLegend) {
		labels <- colnames(expr)
		if (is.null(labels)) labels <- as.character(1:ncol(expr))
		if (reverse) {
			x.range <- range(-expr)
			x.pos <- (max(x.range) - min(x.range)) * 1/10 + min(x.range)
			y.pos <- 0.95
		} else {
			x.range <- range(expr)
			x.pos <- max(x.range) - (max(x.range) - min(x.range)) * 1/3
			y.pos <- 0.5
		}
		legend(x.pos, y.pos, legend=labels, col=usedCol, lwd=usedLwd)
	}
	
    invisible(TRUE)
}


