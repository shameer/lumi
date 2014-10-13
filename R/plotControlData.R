`plotControlData` <-
function(controlData, type=NULL, slideIndex=NULL, logMode=FALSE, new=TRUE, ...) 
{
	methylationData <- FALSE
	if (is(controlData, "MethyLumiM")) {
		controlData <- controlData(controlData)
	}
	if (is(controlData, "MethyLumiQC")) {
		methylationData <- TRUE
		# For control data, methylated data corresponds to green channel
		# For control data, unmethylated data corresponds to red channel
		grnData <- assayDataElement(controlData, "methylated") 
		redData <- assayDataElement(controlData, "unmethylated")
		allControlType <- sapply(strsplit(featureNames(controlData), "\\."), function(x) x[1])
		uniControlType <- getControlType(controlData)
	} else {
		if (is(controlData, 'LumiBatch')) {
			sampleID <- pData(phenoData(controlData))$sampleID
			controlData <- controlData(controlData)
			if (nrow(controlData) == 0) stop('Slot controlData is empty!')
		} else {
			sampleID <- NULL
		}
		allControlType <- controlData$controlType
		uniControlType <- getControlType(controlData)
		controlData <- controlData[, -c(1,2)]

		if (is.null(slideIndex)) {
			if (is.null(sampleID)) sampleID <- colnames(controlData)
			sampleIDInfo <- strsplit(sampleID, split="_")
			chipID <- sapply(sampleIDInfo, function(x) x[1])
			ord <- order(chipID)
			chipID <- chipID[ord]
			controlData <- controlData[, ord]
			col <- as.numeric(as.factor(chipID))
			chipNum <- length(unique(chipID))
		} else {
			ord <- order(slideIndex)
			controlData <- controlData[, ord]
			slideIndex <- slideIndex[ord]
			col <- as.numeric(as.factor(slideIndex))
			chipNum <- length(unique(slideIndex))
		}

	}
	
	if (!is.null(type)) {
		type <- toupper(type)
		allControlType <- toupper(allControlType)
		if (!all(type %in% uniControlType)) {
			warning('"type" does not match "controlType". Please use getControlType function to view the available type!')
			type <- ''
		} 
	} else {
		type <- ''
	}
	
	if (methylationData) {
		## dealing with the methylation control data
		if (logMode) {
			if (max(c(redData, grnData)) > 50) {
				redData[redData < 1] <- 1
				grnData[grnData < 1] <- 1

				redData <- log2(redData)
				grnData <- log2(grnData)
			}
			ylab <- 'Intensity (log2)'
		} else {
			if (max(c(redData, grnData)) < 50) {
				redData <- 2^(redData)
				grnData <- 2^(grnData)
			} 
			ylab <- 'Intensity'
		}

		if (length(type) > 1) oldPar <- par(mfrow=c(ceiling(length(type)/2), 2))
		for (type.i in type) {
			if (type.i != '') {
				selRedData <- redData[allControlType == type.i, , drop=FALSE]
				selGrnData <- grnData[allControlType == type.i, , drop=FALSE]
			} else {
				selRedData <- redData
				selGrnData <- grnData
			}

			if (nrow(selRedData) > 1) {
				mm.red <- apply(selRedData, 2, mean)
				std.red <- apply(selRedData, 2, sd)
				mm.grn <- apply(selGrnData, 2, mean)
				std.grn <- apply(selGrnData, 2, sd)
			} else {
				mm.red <- selRedData
				mm.grn <- selGrnData
				std.red <- std.grn <- rep(0, length(mm))
			}

			labels <- colnames(selRedData)
			if (is.null(labels)) labels <- as.character(1:ncol(selRedData))
			## set the margin of the plot
			mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
			old.mar <- par('mar')
			old.xaxt <- par('xaxt')
			par(xaxt='n')
			par(mar=mar)

			ind <- ncol(selRedData)
			if (new) {
				range <- c(min(c(mm.red - std.red * 1.5, mm.grn - std.red * 1.5)), max(c(mm.red + std.red * 1.1, mm.grn + std.grn * 1.1)))
				plot(1:ind, mm.red, pch = 19, cex=1.8, type='o', col='red', xlab='', ylab=ylab, ylim=range, ...)
				points(1:ind, mm.grn, pch = 19, cex=1.8, type='o', col='green', ...)
				title(type.i)
			} else {
				points(1:ind, mm.red, pch = 19, cex=1.8, type='o', col='red', ...)
				points(1:ind, mm.grn, pch = 19, cex=1.8, type='o', col='green', ...)
			}
			if (nrow(selRedData) > 1) {
				arrows(1:ind, mm.red - std.red, 1:ind, mm.red + std.red, code = 3, col='red', angle = 90, length = .1)
				arrows(1:ind, mm.grn - std.grn, 1:ind, mm.grn + std.grn, code = 3, col='green', angle = 90, length = .1)
			}

			par(xaxt='s')
			axis(1, at=1:ncol(selRedData), labels=labels, tick=TRUE, las=2)
			par(mar=old.mar)
			par(xaxt=old.xaxt)
		}
		if (length(type) > 1) par(oldPar)
		
	} else {
		if (logMode) {
			if (max(controlData) > 50) {
				controlData[controlData < 1] <- 1
				controlData <- log2(controlData)
			}
			ylab <- 'Intensity (log2)'
		} else {
			if (max(controlData) < 50) controlData <- 2^(controlData)
			ylab <- 'Intensity'
		}

		if (length(type) > 1) oldPar <- par(mfrow=c(ceiling(length(type)/2), 2))
		for (type.i in type) {
			if (type.i != '') {
				selControlData <- controlData[allControlType == type.i, , drop=FALSE]
			} else {
				selControlData <- controlData
			}

			if (nrow(selControlData) > 1) {
				mm <- apply(selControlData, 2, mean)
				std <- apply(selControlData, 2, sd)
			} else {
				mm <- selControlData
				std <- rep(0, length(mm))
			}

			labels <- colnames(selControlData)
			if (is.null(labels)) labels <- as.character(1:ncol(selControlData))
			## set the margin of the plot
			mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
			old.mar <- par('mar')
			old.xaxt <- par('xaxt')
			par(xaxt='n')
			par(mar=mar)

			ind <- ncol(selControlData)
			if (new) {
				range <- c(min(mm - std * 1.5), max(mm + std * 1.1))
				plot(1:ind, mm, pch = 19, cex=1.8, type='o', col=col, xlab='', ylab=ylab, ylim=range, ...)
				title(type.i)
			} else {
				points(1:ind, mm, pch = 19, cex=1.8, type='o', col=col, ...)
			}
			if (nrow(selControlData) > 1) {
				arrows(1:ind, mm - std, 1:ind, mm + std, code = 3, col=col, angle = 90, length = .1)
			}
			# plot the vertical lines
			if (chipNum > 1) {
				abline(v=0.5 + which(diff(col) != 0), lty=2)
			}

			par(xaxt='s')
			axis(1, at=1:ncol(selControlData), labels=labels, tick=TRUE, las=2)
			par(mar=old.mar)
			par(xaxt=old.xaxt)
		}
		if (length(type) > 1) par(oldPar)
	}
	
	return(invisible(TRUE))
}

