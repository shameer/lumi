# estimate the intensity measured by Illumina Infinium methylation probes
# which basically is the sum of methylated and unmethylated probe intensity
estimateIntensity <- function(methyLumiM) {
	
	intensity.a <- assayDataElement(methyLumiM, 'unmethylated') 
	intensity.b <- assayDataElement(methyLumiM, 'methylated') 
	if (!is.null(intensity.a) && !is.null(intensity.b)) {
		intensity <- intensity.a + intensity.b
	} else {
		cat("The input data does not include methylated and unmethylated data information!\n")
		intensity <- NULL
	}
	return(intensity)
}

# convert beta-value to m-value
beta2m <- function(beta) {
	m <- log2(beta/(1-beta))
	return(m)
}

# convert m-value to beta-value 
m2beta <- function(m) {
	beta <- 2^m/(2^m + 1)
	return(beta)
}

# estimate the M-value based on methylated and unmethylated probe intensities
estimateM <- function(methyLumiM, minValue=100) {
	
	if (is(methyLumiM, "MethyLumiM")) return(exprs(methyLumiM))
	
	intensity.a <- assayDataElement(methyLumiM, 'unmethylated') 
	intensity.b <- assayDataElement(methyLumiM, 'methylated') 
	if (!is.null(intensity.a) && !is.null(intensity.b)) {
		intensity.a[intensity.a < 1] <- 1
		intensity.b[intensity.b < 1] <- 1
		mm <- min(c(intensity.a, intensity.b))
		if (mm < 100) {
			intensity.a <- intensity.a - mm + minValue
			intensity.b <- intensity.b - mm + minValue
		}
		M <- log2((intensity.b) / (intensity.a))
	} else {
		cat("The input data does not include methylated and unmethylated data information!\n")
		M <- NULL
	}
	return(M)
}

# estimate the Beta-value based on methylated and unmethylated probe intensities
estimateBeta <- function(methyLumiM, offset=100) {
	
	intensity.a <- assayDataElement(methyLumiM, 'unmethylated') 
	intensity.b <- assayDataElement(methyLumiM, 'methylated') 
	if (!is.null(intensity.a) && !is.null(intensity.b)) {
		intensity.a[intensity.a < 1] <- 1
		intensity.b[intensity.b < 1] <- 1
		intensity <- intensity.a + intensity.b
		beta <- intensity.b / (intensity + offset)
	} else {
		cat("The input data does not include methylated and unmethylated data information!\n")
		beta <- NULL
	}
	return(beta)
}


## boxplotColorBias
# boxplotColorBias(methyLumiM)
boxplotColorBias <- function(methyLumiM, logMode=TRUE, channel=c('both', 'unmethy', 'methy', 'sum'), grid=TRUE, main=NULL, mar=NULL, verbose=F, ...) {
	
	channel <- match.arg(channel)
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color bias adjustment because lack of COLOR_CHANNEL information!\n")
		return(invisible(FALSE))
	}
	redInd <- annotation$COLOR_CHANNEL == 'Red'
	grnInd <- annotation$COLOR_CHANNEL == 'Grn'
	intensity.a <- assayDataElement(methyLumiM, 'unmethylated') 
	intensity.b <- assayDataElement(methyLumiM, 'methylated') 
	
	if (logMode) {
		intensity.a[intensity.a < 1] <- 1
		intensity.b[intensity.b < 1] <- 1
		intensity.a <- log2(intensity.a)
		intensity.b <- log2(intensity.b)
	}
	
	nSample <- ncol(intensity.a)
	allRed <- allGrn <- NULL
	tmp <- lapply(1:nSample, function(i) {
		red.a.i <- intensity.a[redInd, i]
		red.b.i <- intensity.b[redInd, i]
		grn.a.i <- intensity.a[grnInd, i]
		grn.b.i <- intensity.b[grnInd, i]
		if (channel == 'both') {
			red <- c(red.a.i, red.b.i)
			grn <- c(grn.a.i, grn.b.i)
		} else if (channel == 'unmethy') {
			red <- red.a.i
			grn <- grn.a.i
		} else if (channel == 'methy') {
			red <- red.b.i
			grn <- grn.b.i
		} else {
			red <- red.a.i + red.b.i
			grn <- grn.a.i + grn.b.i
		}
		allRed <<- cbind(allRed, red)
		allGrn <<- cbind(allGrn, grn)
		return(NULL)
	})
	labels <- colnames(intensity.a)
	if (is.null(main)) {
		info <- switch(channel, 
			"both"="(both methylated and unmethylated probes)",
			"methy"="(methylated probes only)",
			"unmethy"="(unmethylated probes only)",
			"sum"="(Sum of methylated and unmethylated probes)")
		main <- paste("Boxplots of Red and Green color channels", info)
	}
	if (is.null(mar)) {
		mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
		mar[mar > 15] <- 15
		if (verbose) cat("mar:", mar, "\n")
	}
	old.par <- par(mar=mar, xaxt='n')
	ylab <- ifelse(logMode, "Log2( Intensity )", "Intensity")
	boxplot(allRed ~ col(allRed), col = "red", boxwex = 0.25, at = 1:nSample - 0.175,  ylab=ylab, xlab="", main=main, ...)
	boxplot(allGrn ~ col(allGrn), col = "green", boxwex = 0.25, at = 1:nSample + 0.175, axis=F, add=TRUE, ylab="", xlab="")
	par(xaxt='s')
	axis(1, at=1:ncol(allRed), labels=labels, tick=TRUE, las=2)
	par(old.par)
	if (grid) abline(v=1:(nSample-1) + 0.5, col = "lightgray", lty = "dotted")
	
	return(invisible(TRUE))
}


# plotColorBias1D(methyLumiM)
## plot either density or scatter plot of two color channels
plotColorBias1D <- function(methyLumiM, removeGenderProbes=FALSE, logMode=TRUE, channel=c('both', 'unmethy', 'methy', 'sum'), xlim=NULL, ...) {

	channel <- match.arg(channel)
	intensity.a <- assayDataElement(methyLumiM, 'unmethylated') 
	intensity.b <- assayDataElement(methyLumiM, 'methylated') 
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color bias adjustment because lack of COLOR_CHANNEL information!\n")
		return(methyLumiM)
	}
	if (removeGenderProbes && !is.null(annotation$CHR)) {
		redInd <- annotation$COLOR_CHANNEL == 'Red' & !(annotation$CHR %in% c('X', 'Y'))
		grnInd <- annotation$COLOR_CHANNEL == 'Grn' & !(annotation$CHR %in% c('X', 'Y'))
	} else {
		redInd <- annotation$COLOR_CHANNEL == 'Red'
		grnInd <- annotation$COLOR_CHANNEL == 'Grn'
	}
	
	nSample <- ncol(intensity.a)
	density.list <- lapply(1:nSample, function(i) {
		red.a.i <- intensity.a[redInd, i]
		red.b.i <- intensity.b[redInd, i]
		grn.a.i <- intensity.a[grnInd, i]
		grn.b.i <- intensity.b[grnInd, i]
		if (channel == 'both') {
			red <- c(red.a.i, red.b.i)
			grn <- c(grn.a.i, grn.b.i)
		} else if (channel == 'unmethy') {
			red <- red.a.i
			grn <- grn.a.i
		} else if (channel == 'methy') {
			red <- red.b.i
			grn <- grn.b.i
		} else {
			red <- red.a.i + red.b.i
			grn <- grn.a.i + grn.b.i
		}
		if (logMode) {
			dd.red <- density(log2(red))
			dd.grn <- density(log2(grn), ...)
		} else {
			dd.red <- density(red, ...)
			dd.grn <- density(grn, ...)
		}
		return(list(red=dd.red, green=dd.grn))
	})
	
	mm.density <- max(sapply(density.list, function(x) max(c(x$red$y, x$green$y))))
	xrange <- range(sapply(density.list, function(x) range(c(x$red$x, x$green$x))))
	
	for (i in 1:nSample) {
		dd.red.i <- density.list[[i]]$red
		dd.grn.i <- density.list[[i]]$green
		xlab <- ifelse(logMode, "Intensity (log2)", "Intensity")
		if (is.null(xlim)) xlim = xrange
		if (i == 1) {
			plot(dd.red.i, type='l', col=2, ylim=c(0,mm.density), xlim=xlim, xlab=xlab, ylab='Density', main="Compare density distribution of two color channels")
			lines(dd.grn.i, col=3)
		} else {
			lines(dd.red.i, col=2, lty=i)
			lines(dd.grn.i, col=3, lty=i)
		}
	}
	return(invisible(TRUE))
}

# plotColorBias(methyLumiM, selSample=1)
plotColorBias2D <- function(methyLumiM, selSample=1, combineMode=F, layoutRatioWidth=c(0.75,0.25), layoutRatioHeight=c(0.25, 0.75), margins = c(5, 5, 2, 2), cex=1.25, logMode=TRUE, main='') {

	ff <- pData(featureData(methyLumiM))
	color.channel <- ff[,"COLOR_CHANNEL"]
	if (is.null(color.channel) && !combineMode) stop("No color channel information included in the data!\n")
	
	#intensity.a <- unmethylated(methyLumiM)[, selSample[1]]
	#intensity.b <- methylated(methyLumiM)[, selSample[1]]
	intensity.a <- assayData(methyLumiM)[['unmethylated']][, selSample[1]]
	intensity.b <- assayData(methyLumiM)[['methylated']][, selSample[1]]

	if (logMode) {
		intensity.a <- log2(intensity.a)
		intensity.b <- log2(intensity.b)
	}
	
	if (!combineMode) {
		grn.a <- intensity.a[color.channel == "Grn"]
		red.a <- intensity.a[color.channel == "Red"]
		grn.b <- intensity.b[color.channel == "Grn"]
		red.b <- intensity.b[color.channel == "Red"]
	}

    layout(matrix(c(2,1,0,3), nrow=2), widths = layoutRatioWidth, heights = layoutRatioHeight, respect = FALSE)
	# layout.show(3)
	## plot the scatter plot 
	oldpar <- par(mar = c(margins[1], margins[2], 0, 0))
	on.exit(oldpar)
	plottype <- ifelse(combineMode, "p", "n")
	if (logMode) {
		plot(intensity.a, intensity.b, xlim=range(intensity.a), ylim=range(intensity.b), type=plottype, pch='.', xlab='Unmethylated Probe Intensity (log2)', ylab='Methylated Probe Intensity (log2)')
	} else {
		plot(intensity.a, intensity.b, xlim=range(intensity.a), ylim=range(intensity.b), type=plottype, pch='.', xlab='Unmethylated Probe Intensity', ylab='Methylated Probe Intensity')
	}
	if (combineMode) {
		dd.a <- density(intensity.a)
		dd.b <- density(intensity.b)
		par(mar = c(1, margins[2], margins[3], 0))
		plot(dd.a, xlab='', ylab='Density', xlim=range(intensity.a), xaxt='n', col='black', type='l', main='')
		## plot the density plot of methylated probes
		par(mar = c(margins[1], 1, 0, margins[4]))
		plot(dd.b$y, dd.b$x, xlab='Density', ylab='', ylim=range(intensity.a), yaxt='n', col='black', type='l', main='')
	} else {
		points(red.a, red.b, pch='.', cex=cex, col='red')
		points(grn.a, grn.b, pch='.', cex=cex, col='green')		

		dd.grn.a <- density(grn.a)
		dd.grn.b <- density(grn.b)
		dd.red.a <- density(red.a)
		dd.red.b <- density(red.b)
		## plot the density plot of unmethylated probes
		par(mar = c(1, margins[2], margins[3], 0))
		plot(dd.red.a, xlab='', ylab='Density', xlim=range(intensity.a), ylim=range(c(dd.grn.a$y, dd.red.a$y)), xaxt='n', col='red', type='l', main='')
		lines(dd.grn.a, col='green')
		## plot the density plot of methylated probes
		par(mar = c(margins[1], 1, 0, margins[4]))
		plot(dd.red.b$y, dd.red.b$x, xlab='Density', ylab='', xlim=range(c(dd.grn.b$y, dd.red.b$y)), ylim=range(intensity.a), yaxt='n', col='red', type='l', main='')
		lines(dd.grn.b$y, dd.grn.b$x, col='green')
	}

	par(mfrow=c(1,1))

	return(invisible(TRUE))
}


colorBiasSummary <- function(methyLumiM, logMode=TRUE, channel=c('both', 'unmethy', 'methy', 'sum')) {

	channel <- match.arg(channel)
	intensity.a <- assayDataElement(methyLumiM, 'unmethylated') 
	intensity.b <- assayDataElement(methyLumiM, 'methylated') 
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color bias adjustment because lack of COLOR_CHANNEL information!\n")
		return(methyLumiM)
	}
	redInd <- annotation$COLOR_CHANNEL == 'Red'
	grnInd <- annotation$COLOR_CHANNEL == 'Grn'
	
	nSample <- ncol(intensity.a)
	summary.list <- lapply(1:nSample, function(i) {
		red.a.i <- intensity.a[redInd, i]
		red.b.i <- intensity.b[redInd, i]
		grn.a.i <- intensity.a[grnInd, i]
		grn.b.i <- intensity.b[grnInd, i]
		if (channel == 'both') {
			red <- c(red.a.i, red.b.i)
			grn <- c(grn.a.i, grn.b.i)
		} else if (channel == 'unmethy') {
			red <- red.a.i
			grn <- grn.a.i
		} else if (channel == 'methy') {
			red <- red.b.i
			grn <- grn.b.i
		} else {
			red <- red.a.i + red.b.i
			grn <- grn.a.i + grn.b.i
		}
		if (logMode) {
			ss.red <- summary(log2(red))
			ss.grn <- summary(log2(grn))
		} else {
			ss.red <- summary(red)
			ss.grn <- summary(grn)
		}
		return(list(red=ss.red, green=ss.grn))
	})
	red.summary.matrix <- sapply(summary.list, function(x) x$red)
	grn.summary.matrix <- sapply(summary.list, function(x) x$green)
	
	return(list(red=red.summary.matrix, green=grn.summary.matrix))
}
