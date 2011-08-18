lumiMethyR <- function(..., lib=NULL, controlData=NULL) {
	methyLumiSet <- methylumiR(...)
	methyLumiM <- as(methyLumiSet, "MethyLumiM")
	if (!is.null(lib)) {
		# methyLumiM <- addColorChannelInfo(methyLumiM, lib=lib)
		methyLumiM <- addAnnotationInfo(methyLumiM, lib=lib)
	}
	if (!is.null(controlData)) {
		if (is.character(controlData)) {
			controlData <- methylumiR(controlData)
		}
		if (is(controlData, "MethyLumiQC")) {
			controlData(methyLumiM) <- controlData
		} else {
			cat("Provided controlData is not supported!\n")
		}
	}
	return(methyLumiM)
}


addAnnotationInfo <- function(methyLumiM, lib=NULL, hgVersion=c('hg19', 'hg18'), annotationColumn=c('COLOR_CHANNEL', 'CHROMOSOME', 'POSITION')) {
	
	hgVersion <- match.arg(hgVersion)
	# retrieve feature data
	ff <- fData(methyLumiM)
	if (is.null(ff$COLORCHANNEL)) {
		if (is.null(lib)) stop("Please provide the annotation library!\n")
		if (!require(lib, character=TRUE)) stop(paste(lib, "is not available!\n"))
		
		obj <- get(paste(sub("\\.db$", "", lib), "COLORCHANNEL", sep=""))
		colorInfo <- as.list(obj[mappedkeys(obj)])
		colorInfo <- colorInfo[rownames(ff)]
		ff$COLOR_CHANNEL <- sapply(colorInfo, function(x) x[1])
		
	} 
	
	if (length(grep('27k', lib, ignore.case=T)) > 0) {
		fData(methyLumiM) <- ff	
		return(methyLumiM)
	}
	## Only for 450K data, chromosome information will be added.
	if (hgVersion == 'hg18') {
		if (all(c('CHROMOSOME_36', 'COORDINATE_36') %in% names(ff))) {
			ff$CHROMOSOME <- ff$CHROMOSOME_36
			ff$POSITION <- as.numeric(ff$COORDINATE_36)
		} else {
			if (is.null(lib)) {
				warning("Please provide the annotation library!\n")
				return(methyLumiM)
			}
			if (!require(lib, character=TRUE)) stop(paste(lib, "is not available!\n"))
		
			obj <- get(paste(sub("\\.db$", "", lib), "CHR36", sep=""))
			chr <- as.list(obj[mappedkeys(obj)])
			chr <- chr[rownames(ff)]
			## only retrieve the first element
			ff$CHROMOSOME <- sapply(chr, function(x) x[1])

			obj <- get(paste(sub("\\.db$", "", lib), "CPG36", sep=""))
			loc <- as.list(obj[mappedkeys(obj)])
			loc <- loc[rownames(ff)]
			ff$POSITION <- sapply(loc, function(x) x[1])
		}
	} else {
		if (all(c('CHR', 'MAPINFO') %in% names(ff))) {
			ff$CHROMOSOME <- ff$CHR
			ff$POSITION <- as.numeric(ff$MAPINFO)
		} else {
			if (is.null(lib)) {
				warning("Please provide the annotation library!\n")
				return(methyLumiM)
			}
			if (!require(lib, character=TRUE)) stop(paste(lib, "is not available!\n"))
		
			obj <- get(paste(sub("\\.db$", "", lib), "CHR37", sep=""))
			chr <- as.list(obj[mappedkeys(obj)])
			chr <- chr[rownames(ff)]
			ff$CHROMOSOME <- sapply(chr, function(x) x[1])

			obj <- get(paste(sub("\\.db$", "", lib), "CPG37", sep=""))
			loc <- as.list(obj[mappedkeys(obj)])
			loc <- loc[rownames(ff)]
			ff$POSITION <- sapply(loc, function(x) x[1])
		}
	}
	
	fData(methyLumiM) <- ff	
	return(methyLumiM)
}


addControlData2methyLumiM <- function(controlData, methyLumiM, checkConsistency = TRUE, ...) 
{
	if (missing(methyLumiM) || missing(controlData)) stop('Both controlData and methyLumiM are required!')
	if (is.character(controlData)) {
		controlData <- methylumiR(controlData, ...)
	}
	if (is(controlData, "MethyLumiQC")) {
		## Match the column names of controlData and LumiBatch object 
		## Only keep the samples matching methyLumiM
		if (checkConsistency) {
			sampleID <- sampleNames(methyLumiM)
			controlSampleID <- sampleNames(controlData)
			if (all(sampleID %in% controlSampleID)) {
				controlData(methyLumiM) <- controlData[, sampleID]
			} else {
				stop('SampleNames do not match up between controlData and methyLumiM!')
			}
		} else {
			controlData(methyLumiM) <- controlData[, sampleID]
		}
	} else {
		cat("Provided controlData is not supported!\n")
	}

	return(methyLumiM)
}



# normalization
lumiMethyN <- function(methyLumiM, method = c('ssn', 'quantile', 'none'), separateColor=FALSE, verbose=TRUE, ...) 
{
	if (!is.function(method)) method <- match.arg(method)

	if (!is(methyLumiM, 'MethyLumiM')) {
		stop('The object should be class "MethyLumiM" inherited!')
	}
	
	history.submitted <- as.character(Sys.time())
	if (!(is.function(method))) {
		if (!(method %in% c('ssn', 'rssn', 'quantile', 'none'))) {
			cat('This method is not supported!\n')
			return(methyLumiM)
		} else if (method == 'none') {
			return(methyLumiM)
		} 		
		if (verbose) cat(paste('Perform', method, 'normalization ...\n'))
	} else {
		if (verbose) cat('Perform user provided normalization ...\n')		
	}
			
	if (is.function(method)) {
		unmethy <- assayDataElement(methyLumiM, 'unmethylated')
		methy <- assayDataElement(methyLumiM, 'methylated')
		if (separateColor) {
			annotation <- pData(featureData(methyLumiM))
			if (is.null(annotation$COLOR_CHANNEL)) {
				cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
				combData <- rbind(methy, unmethy)
				processed.comb <- method(combData, ...)
				if (!is(processed.comb, 'matrix')) stop("The return of user defined method should be a matrix!\n")
				methy <- processed.comb[1:(nrow(processed.comb)/2), ]
				unmethy <- processed.comb[(nrow(processed.comb)/2+1):nrow(processed.comb), ]
			} else {
				allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
				allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')
				allRed <- rbind(methy[allRedInd,], unmethy[allRedInd,])
				allGrn <- rbind(methy[allGrnInd,], unmethy[allGrnInd,])
				processed.red <- method(allRed, ...)
				if (!is(processed.red, 'matrix')) stop("The return of user defined method should be a matrix!\n")
				processed.grn <- method(allGrn, ...)
				methy[allRedInd,] <- processed.red[1:(nrow(processed.red)/2), ]
				unmethy[allRedInd,] <- processed.red[(nrow(processed.red)/2+1):nrow(processed.red), ]
				methy[allGrnInd,] <- processed.grn[1:(nrow(processed.grn)/2), ]
				unmethy[allGrnInd,] <- processed.grn[(nrow(processed.grn)/2+1):nrow(processed.grn), ]
			}
		} else {
			combData <- rbind(methy, unmethy)
			processed.comb <- method(combData, ...)
			if (!is(processed.comb, 'matrix')) stop("The return of user defined method should be a matrix!\n")
			methy <- processed.comb[1:(nrow(processed.comb)/2), ]
			unmethy <- processed.comb[(nrow(processed.comb)/2+1):nrow(processed.comb), ]
		}
		methylated(methyLumiM) <- methy
		unmethylated(methyLumiM) <- unmethy
		methyLumiM <- estimateM(methyLumiM)
	} else {
		if (method == 'quantile') {
			methyLumiM <- normalizeMethylation.quantile(methyLumiM, separateColor=separateColor, ...)
		} else if (method == 'ssn') {
			methyLumiM <- normalizeMethylation.ssn(methyLumiM, separateColor=separateColor, ...)
		#} else if (method == 'rssn') {
		#	methyLumiM <- normalizeMethylation.robust.ssn(methyLumiM, separateColor=separateColor, ...)
		}
	}
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiMethyB)))
	
	lumiVersion <- packageDescription('lumi')$Version
	methyLumiM@history<- rbind(methyLumiM@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))
			
	return(methyLumiM)
}


# color balance adjustment
lumiMethyC <- function(methyLumiM, method = c('quantile', 'ssn', 'none'), verbose=TRUE, ...) 
{
	if (!is.function(method)) method <- match.arg(method)

	if (!is(methyLumiM, 'MethyLumiM')) {
		stop('The object should be class "MethyLumiM" inherited!')
	}
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL))
		cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
	
	history.submitted <- as.character(Sys.time())
	if (!(is.function(method))) {
		if (!(method %in% c('quantile', 'ssn', 'none'))) {
			cat('This method is not supported!\n')
			return(methyLumiM)
		} else if (method == 'none') {
			return(methyLumiM)
		} 		
		if (verbose) cat(paste('Perform', method, 'color balance adjustment ...\n'))
	} else {
		if (verbose) cat('Perform user provided color balance adjustment ...\n')		
	}
			
	if (is.function(method)) {
		unmethy <- assayDataElement(methyLumiM, 'unmethylated')
		methy <- assayDataElement(methyLumiM, 'methylated')
		allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
		allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')
		allRed <- rbind(methy[allRedInd,], unmethy[allRedInd,])
		allGrn <- rbind(methy[allGrnInd,], unmethy[allGrnInd,])

		processedData <- method(allRed, allGrn, ...)
		processed.red <- processedData$red
		if (is.null(processed.red)) stop("The return of user defined method should be a list including 'red' and 'green' matrix!\n")
		processed.grn <- processedData$green
		methy[allRedInd,] <- processed.red[1:(nrow(processed.red)/2), ]
		unmethy[allRedInd,] <- processed.red[(nrow(processed.red)/2+1):nrow(processed.red), ]
		methy[allGrnInd,] <- processed.grn[1:(nrow(processed.grn)/2), ]
		unmethy[allGrnInd,] <- processed.grn[(nrow(processed.grn)/2+1):nrow(processed.grn), ]
		methylated(methyLumiM) <- methy
		unmethylated(methyLumiM) <- unmethy
		methyLumiM <- estimateM(methyLumiM)
	} else {
		if (method == 'quantile') {
			methyLumiM <- adjColorBias.quantile(methyLumiM, ...)
		} else if (method == 'ssn') {
			methyLumiM <- adjColorBias.ssn(methyLumiM, ...)
		} 
	}
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiMethyB)))
	
	lumiVersion <- packageDescription('lumi')$Version
	methyLumiM@history<- rbind(methyLumiM@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))
			
	return(methyLumiM)
}


lumiMethyB <- function(methyLumiM, method = c('bgAdjust2C', 'forcePositive', 'none'), separateColor=FALSE, verbose=TRUE, ...) 
{
	if (!is.function(method)) method <- match.arg(method)

	if (is(methyLumiM, 'MethyLumiM')) {
		unmethy <- assayDataElement(methyLumiM, 'unmethylated')
		methy <- assayDataElement(methyLumiM, 'methylated')
	} else {
		stop('The object should be class "MethyLumiM" inherited!')
	}
	
	history.submitted <- as.character(Sys.time())
	if (!(is.function(method))) {
		if (!(method %in% c('bgAdjust2C', 'none', 'forcePositive'))) {
			cat('This method is not supported!\n')
			return(methyLumiM)
		} else if (method == 'none') {
			return(methyLumiM)
		} else if (method == 'forcePositive') {
			mm <- min(c(min(methy), min(unmethy)))
			if (mm > 0)	return(methyLumiM)
		}		
		if (verbose) cat(paste('Perform', method, 'background correction ...\n'))
	} else {
		if (verbose) cat('Perform user provided background correction ...\n')		
	}
			
	if (is.function(method)) {
		if (separateColor) {
			annotation <- pData(featureData(methyLumiM))
			if (is.null(annotation$COLOR_CHANNEL)) {
				cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
				combData <- rbind(methy, unmethy)
				bgAdj.comb <- method(combData, ...)
				if (!is(bgAdj.comb, 'matrix')) stop("The return of user defined method should be a matrix!\n")
				methy <- bgAdj.comb[1:(nrow(bgAdj.comb)/2), ]
				unmethy <- bgAdj.comb[(nrow(bgAdj.comb)/2+1):nrow(bgAdj.comb), ]
			} else {
				allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
				allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')
				allRed <- rbind(methy[allRedInd,], unmethy[allRedInd,])
				allGrn <- rbind(methy[allGrnInd,], unmethy[allGrnInd,])
				bgAdj.red <- method(allRed, ...)
				if (!is(bgAdj.red, 'matrix')) stop("The return of user defined method should be a matrix!\n")
				bgAdj.grn <- method(allGrn, ...)
				methy[allRedInd,] <- bgAdj.red[1:(nrow(bgAdj.red)/2), ]
				unmethy[allRedInd,] <- bgAdj.red[(nrow(bgAdj.red)/2+1):nrow(bgAdj.red), ]
				methy[allGrnInd,] <- bgAdj.grn[1:(nrow(bgAdj.grn)/2), ]
				unmethy[allGrnInd,] <- bgAdj.grn[(nrow(bgAdj.grn)/2+1):nrow(bgAdj.grn), ]
			}
		} else {
			combData <- rbind(methy, unmethy)
			bgAdj.comb <- method(combData, ...)
			if (!is(bgAdj.comb, 'matrix')) stop("The return of user defined method should be a matrix!\n")
			methy <- bgAdj.comb[1:(nrow(bgAdj.comb)/2), ]
			unmethy <- bgAdj.comb[(nrow(bgAdj.comb)/2+1):nrow(bgAdj.comb), ]
		}
		methylated(methyLumiM) <- methy
		unmethylated(methyLumiM) <- unmethy
	} else {
		if (method == 'bgAdjust2C') {
			methyLumiM <- bgAdjustMethylation(methyLumiM, separateColor=separateColor, ...)
		} else if (method == 'forcePositive') {
			x.matrix <- rbind(methy, unmethy)
			offset <- apply(x.matrix, 2, min, na.rm=TRUE)
			offset[offset <= 0] <- offset[offset <= 0] - 0.01
			offset[offset > 0] <- 0
			offset <- rep(1, nrow(x.matrix)) %*% t(offset)
			x.matrix <- x.matrix - offset
			methylated(methyLumiM) <- x.matrix[1:nrow(methy),]
			unmethylated(methyLumiM) <- x.matrix[(nrow(unmethy)+1):nrow(x.matrix),]
		} 
	}
	methyLumiM <- estimateM(methyLumiM)
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiMethyB)))
	
	lumiVersion <- packageDescription('lumi')$Version
	methyLumiM@history<- rbind(methyLumiM@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))
			
	return(methyLumiM)
}


bgAdjustMethylation <- function(methyLumiM, separateColor=FALSE, targetBGLevel=300, negPercTh=0.25) {
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	methy <- methylated(methyLumiM)
	unmethy <- unmethylated(methyLumiM)
	
	## check whether control data is available in the MethyLumiM object
	bglevel <- estimateMethylationBG(methyLumiM, separateColor=separateColor)
	if (is.null(colnames(bglevel))) separateColor <- FALSE
	if (separateColor) {
		bg.red <- bglevel[,'red']
		bg.grn <- bglevel[,'green']
		
		annotation <- pData(featureData(methyLumiM))
		allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
		allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')

		## For 450K chip, the Infinium II type only has one type of probe for each CpG site. 
		## The color channel will depend on the methylation status
		allBothInd <- which(annotation$COLOR_CHANNEL == 'Both' | annotation$COLOR_CHANNEL == '')
		
		methy[allRedInd, ] <- methy[allRedInd, ] - rep(1, length(allRedInd)) %*% t(bg.red)
		methy[allGrnInd, ] <- methy[allGrnInd, ] - rep(1, length(allGrnInd)) %*% t(bg.grn)
		unmethy[allRedInd, ] <- unmethy[allRedInd, ] - rep(1, length(allRedInd)) %*% t(bg.red)
		unmethy[allGrnInd, ] <- unmethy[allGrnInd, ] - rep(1, length(allGrnInd)) %*% t(bg.grn)
		
		if (length(allBothInd) > 0) {
			## the methylated probe has 'Grn' color, while the unmethylated probe has 'Red' color
			unmethy[allBothInd, ] <- unmethy[allBothInd,] - rep(1, length(allBothInd)) %*% t(bg.red) 
			methy[allBothInd, ] <- methy[allBothInd,] - rep(1, length(allBothInd)) %*% t(bg.grn) 
		}
	} else {
		if (is.matrix(bglevel)) bglevel <- rowMeans(bglevel)
		methy <- methy - rep(1, nrow(methy)) %*% t(bglevel)
		unmethy <- unmethy - rep(1, nrow(unmethy)) %*% t(bglevel)
	}

	## check possible error of BG estimation model
	negPerc.methy <- apply(methy, 2, function(x) length(which(x < 0))/nrow(methy))
	negPerc.unmethy <- apply(unmethy, 2, function(x) length(which(x < 0))/nrow(unmethy))
	negPerc <- pmax(negPerc.methy, negPerc.unmethy)

	if (any(negPerc > negPercTh)) {
		warning("Possible bad quality samples or possible error of Background estimation model: \n")
		cat("\t", paste(colnames(methy)[negPerc > negPercTh], collapse=", "), "\n")
	}
	methylated(methyLumiM) <- methy + targetBGLevel
	unmethylated(methyLumiM) <- unmethy + targetBGLevel

	if (is(methyLumiM, "MethyLumiM")) methyLumiM <- estimateM(methyLumiM)
	
	attr(methyLumiM, "EstimatedBG") <- list(rawBG=bglevel, adjBG=rep(targetBGLevel, ncol(methy)))	
	return(methyLumiM)
}



## ---------------------------------------------------------------
# methy.raw1.adj <- adjColorBias(methy.raw1)
# methy.raw2.adj <- adjColorBias(methy.raw2)
# methy.adj <- cbind(methylated(methy.raw1.adj), methylated(methy.raw2.adj))
# unmethy.adj <- cbind(unmethylated(methy.raw1.adj), unmethylated(methy.raw2.adj))
# M.adj <- log2(methy.adj/unmethy.adj)
adjColorBias.ssn <- function(methyLumiM, refChannel=c("green", "red", "mean")) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	if (storageMode(assayData(methyLumiM)) == "environment") storageMode(assayData(methyLumiM)) <- "lockedEnvironment"
	refChannel <- match.arg(refChannel)
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	if (is.null(unmethy) || is.null(methy)) stop("methylated or unmethylated data is not included in the dataset!\n")
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
		return(methyLumiM)
	}
	
	allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
	allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')

	## For 450K chip, the Infinium II type only has one type of probe for each CpG site. 
	## The color channel will depend on the methylation status
	allBothInd <- which(annotation$COLOR_CHANNEL == 'Both' | annotation$COLOR_CHANNEL == '')
	
	bg <- estimateMethylationBG(methyLumiM, separateColor=TRUE)
	bg.red <- bg[,"red"]
	bg.grn <- bg[,"green"]
	
	intensity.grn <- unmethy[allGrnInd,] + methy[allGrnInd,]
	intensity.red <- unmethy[allRedInd,] + methy[allRedInd,]
	m.int.grn <- colMeans(intensity.grn) - bg.grn
	m.int.red <- colMeans(intensity.red) - bg.red
	if (refChannel == 'green') {
		m.ref <- m.int.grn
		bg.ref <- bg.grn
	} else if (refChannel == 'red') {
		m.ref <- m.int.red
		bg.ref <- bg.red
	} else {
		m.ref <- (m.int.grn + m.int.red)/2
		bg.ref <- (bg.grn + bg.red)/2
	}

	unmethy[allRedInd, ] <- (unmethy[allRedInd,] - rep(1, length(allRedInd)) %*% t(bg.red)) * 
		(rep(1, length(allRedInd)) %*% t(m.ref/m.int.red)) + rep(1, length(allRedInd)) %*% t(bg.ref)
	unmethy[allGrnInd, ] <- (unmethy[allGrnInd,] - rep(1, length(allGrnInd)) %*% t(bg.grn)) * 
		(rep(1, length(allGrnInd)) %*% t(m.ref/m.int.grn)) + rep(1, length(allGrnInd)) %*% t(bg.ref)
	methy[allRedInd, ] <- (methy[allRedInd,] - rep(1, length(allRedInd)) %*% t(bg.red)) * 
		(rep(1, length(allRedInd)) %*% t(m.ref/m.int.red)) + rep(1, length(allRedInd)) %*% t(bg.ref)
	methy[allGrnInd, ] <- (methy[allGrnInd,] - rep(1, length(allGrnInd)) %*% t(bg.grn)) * 
		(rep(1, length(allGrnInd)) %*% t(m.ref/m.int.grn)) + rep(1, length(allGrnInd)) %*% t(bg.ref)
	
	if (length(allBothInd) > 0) {
		## the methylated probe has 'Grn' color, while the unmethylated probe has 'Red' color
		unmethy[allBothInd, ] <- (unmethy[allBothInd,] - rep(1, length(allBothInd)) %*% t(bg.red)) * 
			(rep(1, length(allBothInd)) %*% t(m.ref/m.int.red)) + rep(1, length(allBothInd)) %*% t(bg.ref)
		methy[allBothInd, ] <- (methy[allBothInd,] - rep(1, length(allBothInd)) %*% t(bg.grn)) * 
			(rep(1, length(allBothInd)) %*% t(m.ref/m.int.grn)) + rep(1, length(allBothInd)) %*% t(bg.ref)
	}
		
	methy.adj <- methyLumiM
	assayDataElement(methy.adj, 'unmethylated') <- unmethy
	assayDataElement(methy.adj, 'methylated') <- methy
	if (is(methy.adj, "MethyLumiM")) methy.adj <- estimateM(methy.adj)
	return(methy.adj)
}


adjColorBias.quantile <- function(methyLumiM, refChannel=c("green", "red"), logMode=TRUE, ...) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	if (storageMode(assayData(methyLumiM)) == "environment") storageMode(assayData(methyLumiM)) <- "lockedEnvironment"
	refChannel <- match.arg(refChannel)
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
		return(methyLumiM)
	}
	redInd <- annotation$COLOR_CHANNEL == 'Red'
	grnInd <- annotation$COLOR_CHANNEL == 'Grn'

	## For 450K chip, the Infinium II type only has one type of probe for each CpG site. 
	## The color channel will depend on the methylation status
	bothInd <- which(annotation$COLOR_CHANNEL == 'Both' | annotation$COLOR_CHANNEL == '')
	bandwidth <- ifelse(logMode, 0.3, 200)
	for (i in 1:ncol(unmethy)) {
		red.a.i <- unmethy[redInd, i]
		red.b.i <- methy[redInd, i]
		red.i <- c(red.a.i, red.b.i)
		grn.a.i <- unmethy[grnInd, i]
		grn.b.i <- methy[grnInd, i]
		grn.i <- c(grn.a.i, grn.b.i)

		if (refChannel == 'green') {
			## For 450K data
			if (length(bothInd) > 0) {
				## the methylated probe has 'Grn' color, while the unmethylated probe has 'Red' color
				y.out.i <- smoothQuantileNormalization(red.i, grn.i, adjData=c(red.i, unmethy[bothInd,i]), logMode=logMode, bandwidth=bandwidth, ...)
				unmethy[bothInd,i] <- y.out.i[(length(red.i)+1):length(y.out.i)]
				# y.out.i <- smoothQuantileNormalization(red.i, grn.i, adjData=c(red.i, methy[bothInd,i]), logMode=logMode, bandwidth=bandwidth, ...)
				# methy[bothInd,i] <- y.out.i[(length(red.i)+1):length(y.out.i)]
			} else {
				y.out.i <- smoothQuantileNormalization(red.i, grn.i, logMode=logMode, bandwidth=bandwidth, ...)
			}
			unmethy[redInd, i] <- y.out.i[1:length(red.a.i)]
			methy[redInd, i] <- y.out.i[(length(red.a.i)+1):length(red.i)]
		} else {
			if (length(bothInd) > 0) {
				## the methylated probe has 'Grn' color, while the unmethylated probe has 'Red' color
				y.out.i <- smoothQuantileNormalization(grn.i, red.i, adjData=c(grn.i, methy[bothInd,i]), logMode=logMode, bandwidth=bandwidth, ...)
				methy[bothInd,i] <- y.out.i[(length(grn.i)+1):length(y.out.i)]
			} else {
				y.out.i <- smoothQuantileNormalization(grn.i, red.i, logMode=logMode, bandwidth=bandwidth, ...)
			}
			unmethy[grnInd, i] <- y.out.i[1:length(grn.a.i)]
			methy[grnInd, i] <- y.out.i[(length(grn.a.i)+1):length(grn.i)]
		}
	}

	methy.adj <- methyLumiM
	assayDataElement(methy.adj, 'unmethylated') <- unmethy
	assayDataElement(methy.adj, 'methylated') <- methy
	if (is(methy.adj, "MethyLumiM")) methy.adj <- estimateM(methy.adj)
	return(methy.adj)
}


smoothQuantileNormalization <- function(dataMatrix, ref=NULL, adjData=NULL, logMode=TRUE, bandwidth=NULL, degree=1, ...) {
	
	# require(KernSmooth)
	bandwidth <- ifelse(logMode, 0.3, 200)
	if (!is.matrix(dataMatrix)) dataMatrix <- matrix(dataMatrix, ncol=1)
	if (is.null(ref)) {
		# normData <- normalize.quantiles.robust(dataMatrix + 0.0)
		normData <- normalize.quantiles(dataMatrix + 0.0)
		ref <- normData[,1]
	} 
	if (!is.null(adjData)) {
		if (!is.matrix(adjData)) adjData <- matrix(adjData, ncol=1)
		if (ncol(dataMatrix) != ncol(adjData))
			stop("The number of columns of adjData should be consistent with dataMatrix!") 
	}
	interpolationMode <- ifelse(length(ref) != nrow(dataMatrix), TRUE, FALSE) 
	
	if (logMode) {
		## remove those equal or less than 0, which are unreliable values
		if (interpolationMode) ref <- ref[ref > 0]
		if (min(ref) < 1) ref <- ref - min(ref) + 1
		ref <- log2(ref)
	}
	len <- nrow(dataMatrix)
	refLen <- length(ref)
	gridsize <- min(min(len, refLen)/2, 1000)

	# In the case of different lengthes between reference and data, interpolation will be performed.
	if (interpolationMode) {
		x <- (1:refLen)/refLen
		y <- sort(ref, decreasing=F)
	} else {
		y <- sort(ref, decreasing=F)
	}
	
	# smoothing the quantile normalization results
	if (is.null(adjData)) {
		normData <- dataMatrix
	} else {
		normData <- adjData
	}
	for (i in 1:ncol(dataMatrix)) {
		profile.i <- dataMatrix[,i]
		if (!is.null(adjData)) {
			adjData.i <- adjData[,i]
		} else {
			adjData.i <- profile.i
		}
		if (logMode) {
			mm.i <- min(c(profile.i, adjData.i))
			if (mm.i < 1) {
				adjData.i <- log2(adjData.i - mm.i + 1)
				profile.i <- log2(profile.i - mm.i + 1)
			} else {
				adjData.i <- log2(adjData.i)
				profile.i <- log2(profile.i)
			}
		}

		if (interpolationMode) {
			## remove those equal or less than 0, which are unreliable values
			selProfile.i <- profile.i[profile.i > 0]
			len <- length(selProfile.i)
			
			# perform linear interpolation when the length of two profiles different
			x.out.i <- rank(selProfile.i)/len 
			y.out.i <- approx(x=x, y=y, xout=x.out.i, method="linear", rule=2)$y		
		} else {
			selProfile.i <- profile.i
			y.out.i <- y[rank(selProfile.i)]
		}
		
		## remove outliers points at two ends
		boundary.i <- quantile(selProfile.i, c(0.3, 0.7))
		lowInd.i <- selProfile.i < boundary.i[1]
		highInd.i <- selProfile.i > boundary.i[2]
		## fit the low segment
		rlm.i <- rlm(selProfile.i[lowInd.i], y.out.i[lowInd.i])
		dd.i <- y.out.i[lowInd.i] - rlm.i$fitted.values
		outlier.ind.low.i <- which(lowInd.i)[which(abs(dd.i) > 5 * sd(dd.i))]
		## fit the high segment
		rlm.i <- rlm(selProfile.i[highInd.i], y.out.i[highInd.i])
		dd.i <- y.out.i[highInd.i] - rlm.i$fitted.values
		outlier.ind.high.i <- which(highInd.i)[which(abs(dd.i) > 5 * sd(dd.i))]
		outlier.ind.i <- c(outlier.ind.low.i, outlier.ind.high.i)
		if (length(outlier.ind.i) > 0) {
			tt <- locpoly(selProfile.i[-outlier.ind.i], y.out.i[-outlier.ind.i], degree=degree,  gridsize=gridsize, bandwidth=bandwidth, ...)
		} else {
			tt <- locpoly(selProfile.i, y.out.i, degree=degree,  gridsize=gridsize, bandwidth=bandwidth, ...)		
		}

		norm.i <- approx(x=tt$x, y=tt$y, xout=adjData.i, rule=2)$y

		# plot(selProfile.i, y.out.i, pch='.')
		# lines(tt, col=2)
		# if (max(norm.i) > 30) browser()

		if (logMode) {
			normData[,i] <- 2^norm.i
		} else {
			normData[,i] <- norm.i
		}
	}

	return(normData)
}


estimateMethylationBG <- function(methyLumiM, separateColor=FALSE, nbin=1000) {
		
	estimateBG <- function(dataMatrix, nbin=1000) {
		bg <- apply(dataMatrix, 2, function(x) {
			hh.x <- hist(x, nbin, plot=FALSE)
			Th <- hh.x$breaks[which.max(hh.x$counts) + 1] * 2
			dd.x <- density(x[x < Th], na.rm=TRUE)
			bg.x <- dd.x$x[which.max(dd.x$y)]		
		})
		return(bg)
	}
	
	## check whether the control data is available
	controlData <- NULL
	if (is(methyLumiM, "MethyLumiM")) {
		controlData <- controlData(methyLumiM)
	}
	if (is(methyLumiM, "MethyLumiSet")) {
		controlData <- QCdata(methyLumiM)
	}
	if (is(controlData, "MethyLumiQC")) {
		# For control data, methylated data corresponds to green channel, "Signal_Grn"
		# For control data, unmethylated data corresponds to red channel, "Signal_Red"
		grnData <- assayDataElement(controlData, "methylated") 
		redData <- assayDataElement(controlData, "unmethylated")
		allControlType <- sapply(strsplit(featureNames(controlData), "\\."), function(x) x[1])
		allControlType <- toupper(allControlType)
		neg.ind <- which(allControlType == "NEGATIVE") 
		if (length(neg.ind) > 0) {
			bg.grn <- apply(grnData[neg.ind, ], 2, median)
			bg.red <- apply(redData[neg.ind, ], 2, median)
			if (separateColor) {
				bg <- cbind(red=bg.red, green=bg.grn)
			} else {
				bg <- pmin(bg.red, bg.grn)
			}
			return(bg)
		}
	}
	
	## In the case the negative control data is not available, the background will be estimated based on the mode positions of unmethylated or methylated distribution (the smaller one)
	if (is(methyLumiM, "eSet")) {
		if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
			stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
		}
		unmethy <- assayDataElement(methyLumiM, 'unmethylated')
		methy <- assayDataElement(methyLumiM, 'methylated')
		if (separateColor) {
			annotation <- pData(featureData(methyLumiM))
			if (is.null(annotation$COLOR_CHANNEL)) {
				cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
				bg.methy <- estimateBG(methy, nbin=nbin)
				bg.unmethy <- estimateBG(unmethy, nbin=nbin)
				bg <- pmin(bg.methy, bg.unmethy)
			} else {
				allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
				allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')
				bg.methy.red <- estimateBG(methy[allRedInd,], nbin=nbin)
				bg.unmethy.red <- estimateBG(unmethy[allRedInd,], nbin=nbin)
				bg.red <- pmin(bg.methy.red, bg.unmethy.red)

				bg.methy.grn <- estimateBG(methy[allGrnInd,], nbin=nbin)
				bg.unmethy.grn <- estimateBG(unmethy[allGrnInd,], nbin=nbin)
				bg.grn <- pmin(bg.methy.grn, bg.unmethy.grn)
				bg <- cbind(red=bg.red, green=bg.grn)
			}
		} else {
			bg.methy <- estimateBG(methy, nbin=nbin)
			bg.unmethy <- estimateBG(unmethy, nbin=nbin)
			bg <- pmin(bg.methy, bg.unmethy)
		}
	} else if (is(methyLumiM, 'matrix')) {
		bg <- estimateBG(methyLumiM, nbin=nbin)
	} else {
		stop("The input data class is not supported!\n")
	}
	
	return(bg)
}



normalizeMethylation.ssn <- function(methyLumiM, separateColor=FALSE) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	if (storageMode(assayData(methyLumiM)) == "environment") storageMode(assayData(methyLumiM)) <- "lockedEnvironment"	
	unmethy <- assayDataElement(methyLumiM, 'unmethylated')
	methy <- assayDataElement(methyLumiM, 'methylated')
	
	# estimate the background based on methy
	bg <- estimateMethylationBG(methyLumiM, separateColor=separateColor)
	
	if (is.null(colnames(bg))) separateColor <- FALSE
	if (separateColor) {
		bg.red <- bg[,'red']
		bg.grn <- bg[,'green']
		meanBg.red <- mean(bg.red)
		meanBg.grn <- mean(bg.grn)

		annotation <- pData(featureData(methyLumiM))
		allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
		allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')

		intensity.red <- unmethy[allRedInd, ] + methy[allRedInd, ]
		totalIntensity.red <- colSums(intensity.red)
		meanTotalIntensity.red <- mean(totalIntensity.red)

		intensity.grn <- unmethy[allGrnInd, ] + methy[allGrnInd, ]
		totalIntensity.grn <- colSums(intensity.grn)
		meanTotalIntensity.grn <- mean(totalIntensity.grn)
		
		methy[allRedInd, ] <- (methy[allRedInd, ] - rep(1, length(allRedInd)) %*% t(bg.red)) * (rep(1, length(allRedInd)) %*% t(meanTotalIntensity.red/totalIntensity.red)) + meanBg.red
		methy[allGrnInd, ] <- (methy[allGrnInd, ] - rep(1, length(allGrnInd)) %*% t(bg.grn)) * (rep(1, length(allGrnInd)) %*% t(meanTotalIntensity.grn/totalIntensity.grn)) + meanBg.grn
		unmethy[allRedInd, ] <- (unmethy[allRedInd, ] - rep(1, length(allRedInd)) %*% t(bg.red)) * (rep(1, length(allRedInd)) %*% t(meanTotalIntensity.red/totalIntensity.red)) + meanBg.red
		unmethy[allGrnInd, ] <- (unmethy[allGrnInd, ] - rep(1, length(allGrnInd)) %*% t(bg.grn)) * (rep(1, length(allGrnInd)) %*% t(meanTotalIntensity.grn/totalIntensity.grn)) + meanBg.grn
	} else {
		meanBg <- mean(bg)

		intensity <- unmethy + methy
		totalIntensity <- colSums(intensity)
		meanTotalIntensity <- mean(totalIntensity)
		unmethy <- (unmethy - rep(1, nrow(intensity)) %*% t(bg)) * (rep(1, nrow(intensity)) %*% t(meanTotalIntensity/totalIntensity)) + meanBg
		methy <- (methy - rep(1, nrow(intensity)) %*% t(bg)) * (rep(1, nrow(intensity)) %*% t(meanTotalIntensity/totalIntensity)) + meanBg
	} 
	
	assayDataElement(methyLumiM, 'unmethylated') <- unmethy
	assayDataElement(methyLumiM, 'methylated') <- methy
	if (is(methyLumiM, "MethyLumiM")) methyLumiM <- estimateM(methyLumiM)
	return(methyLumiM)
}


normalizeMethylation.quantile <- function(methyLumiM, separateColor=FALSE, ...) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	if (storageMode(assayData(methyLumiM)) == "environment") storageMode(assayData(methyLumiM)) <- "lockedEnvironment"
	unmethy <- assayDataElement(methyLumiM, 'unmethylated')
	methy <- assayDataElement(methyLumiM, 'methylated')

	if (separateColor) {
		annotation <- pData(featureData(methyLumiM))
		if (is.null(annotation$COLOR_CHANNEL)) {
			cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
			separateColor <- FALSE
		} else {
			allRedInd <- which(annotation$COLOR_CHANNEL == 'Red')
			allGrnInd <- which(annotation$COLOR_CHANNEL == 'Grn')
		}
	}
	
	if (separateColor) {
		methy.n <- methy; unmethy.n <- unmethy
		methy.red <- methy[allRedInd, ]
		methy.grn <- methy[allGrnInd, ]
		unmethy.red <- unmethy[allRedInd, ]
		unmethy.grn <- unmethy[allGrnInd, ]
		x.matrix.red <- rbind(methy.red, unmethy.red)
		x.matrix.grn <- rbind(methy.grn, unmethy.grn)
		# Normalize the intensity using robust quantile normalization
		x.matrix.red <- normalize.quantiles.robust(x.matrix.red + 0.0, ...)
		x.matrix.grn <- normalize.quantiles.robust(x.matrix.grn + 0.0, ...)
		len.red <- length(allRedInd)
		len.grn <- length(allGrnInd)
		methy.n[allRedInd,] <- x.matrix.red[1:len.red,]
		unmethy.n[allRedInd,] <- x.matrix.red[(len.red+1):nrow(x.matrix.red),]
		methy.n[allGrnInd,] <- x.matrix.grn[1:len.grn,]
		unmethy.n[allGrnInd,] <- x.matrix.grn[(len.grn+1):nrow(x.matrix.grn),]
	} else {
		x.matrix <- rbind(methy, unmethy)
		# Normalize the intensity using robust quantile normalization
		x.matrix <- normalize.quantiles.robust(x.matrix + 0.0, ...)
		methy.n <- x.matrix[1:nrow(methy),]
		unmethy.n <- x.matrix[(nrow(unmethy)+1):nrow(x.matrix),]		
	}
	colnames(methy.n) <- colnames(unmethy.n) <- colnames(methy)
	rownames(methy.n) <- rownames(unmethy.n) <- rownames(methy)

	assayDataElement(methyLumiM, 'unmethylated') <- unmethy.n
	assayDataElement(methyLumiM, 'methylated') <- methy.n
	if (is(methyLumiM, "MethyLumiM")) methyLumiM <- estimateM(methyLumiM)
	return(methyLumiM)
}



# estimate the intensity measured by Illumina Infinium methylation probes
# which basically is the sum of methylated and unmethylated probe intensity
estimateIntensity <- function(methyLumiM, returnType=c("ExpressionSet", "matrix")) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	returnType <- match.arg(returnType)
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	if (!is.null(unmethy) && !is.null(methy)) {
		intensity <- unmethy + methy
	} else {
		cat("The input data does not include methylated and unmethylated data information!\n")
		intensity <- NULL
	}
	if (returnType == "matrix") {
		return(intensity)
	} else {
		methyLumiIntensity <- as(methyLumiM, "ExpressionSet")
		exprs(methyLumiIntensity) <- intensity
		return(methyLumiIntensity)
	}
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
estimateM <- function(methyLumiM, returnType=c("ExpressionSet", "matrix"), offset=100) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	returnType <- match.arg(returnType)
	
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	mm <- min(c(unmethy, methy))
	if (mm < 0.01) {
		unmethy[unmethy < 0.01] <- 0.01 
		methy[methy < 0.01] <- 0.01 
	}
	M <- log2((methy + offset) / (unmethy + offset))
	if (returnType == "matrix") {
		return(M)
	} else {
		exprs(methyLumiM) <- M
		return(methyLumiM)
	}
}

# estimate the Beta-value based on methylated and unmethylated probe intensities
estimateBeta <- function(methyLumiM, returnType=c("ExpressionSet", "matrix"), offset=100) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	returnType <- match.arg(returnType)
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	unmethy[unmethy < 1] <- 1
	methy[methy < 1] <- 1
	intensity <- unmethy + methy
	beta <- methy / (intensity + offset)
	if (returnType == "matrix") {
		return(beta)
	} else {
		methyLumiBeta <- as(methyLumiM, "ExpressionSet")
		exprs(methyLumiBeta) <- beta
		return(methyLumiBeta)
	}
}

# boxplot unique for MethyLumiM-class
setMethod("boxplot",signature(x="MethyLumiM"),
	function(x, main, prob=c(seq(10,90, by=10), 95), col=gray(rev(seq(prob)/length(prob))), ...) {
		
	# if (!require(hdrcde)) stop("Please install the required hdrcde package./n")  

  	tmp <- description(x)
  	if (missing(main) && (is(tmp, "MIAME")))
     	main <- tmp@title

	dataMatrix <- exprs(x)
	labels <- colnames(dataMatrix)
	if (is.null(labels)) labels <- as.character(1:ncol(dataMatrix))
	## set the margin of the plot
	mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
	old.mar <- par('mar')
	old.xaxt <- par('xaxt')
	par(xaxt='n')
	par(mar=mar)

	tmp <- lapply(1:ncol(dataMatrix), function(i) dataMatrix[,i])
	hdr.boxplot(tmp, main=main, xlab='', ylab='M-value', prob=prob, col=col, ...)
	par(xaxt='s')
	axis(1, at=1:ncol(dataMatrix), labels=labels, tick=TRUE, las=2)
	par(mar=old.mar)
	par(xaxt=old.xaxt)
})


## boxplotColorBias
# boxplotColorBias(methyLumiM)
boxplotColorBias <- function(methyLumiM, logMode=TRUE, channel=c('both', 'unmethy', 'methy', 'sum'), grid=TRUE, main=NULL, mar=NULL, verbose=F, subset=NULL, ...) {
	
	channel <- match.arg(channel)
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(methyLumiM), min(subset, nrow(methyLumiM)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(methyLumiM)]
		}
		methyLumiM <- methyLumiM[index,]
	} 

	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
		return(invisible(FALSE))
	}
	redInd <- annotation$COLOR_CHANNEL == 'Red'
	grnInd <- annotation$COLOR_CHANNEL == 'Grn'
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	
	if (logMode) {
		unmethy[unmethy < 1] <- 1
		methy[methy < 1] <- 1
		unmethy <- log2(unmethy)
		methy <- log2(methy)
	}
	
	nSample <- ncol(unmethy)
	allRed <- allGrn <- NULL
	tmp <- lapply(1:nSample, function(i) {
		red.a.i <- unmethy[redInd, i]
		red.b.i <- methy[redInd, i]
		grn.a.i <- unmethy[grnInd, i]
		grn.b.i <- methy[grnInd, i]
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
	labels <- colnames(unmethy)
	if (is.null(main)) {
		main <- "Boxplots of Red and Green color channels"
	}
	if (is.null(mar)) {
		mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
		mar[mar > 15] <- 15
		if (verbose) cat("mar:", mar, "\n")
	}
	old.par <- par(mar=mar, xaxt='n')
	info <- switch(channel, 
		"both"="intensity of both methylated and unmethylated probes",
		"methy"="intensity of methylated probes only",
		"unmethy"="intensity of unmethylated probes only",
		"sum"="CpG-site Intensity")
	if (logMode) {
		ylab <- paste("Log2", info)
	} else {
		substr(info, 1, 1) <- toupper(substr(info, 1, 1))
		ylab <- info
	}
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
plotColorBias1D <- function(methyLumiM, channel=c('both', 'unmethy', 'methy', 'sum'), colorMode=TRUE, removeGenderProbes=FALSE, logMode=TRUE, subset=NULL, ...) {

	channel <- match.arg(channel)
	otherPar <- list(...)
	densityPar <- otherPar[names(otherPar) %in% names(formals(density.default))]
	otherPar[names(otherPar) %in% names(formals(density.default))] <- NULL
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}

	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(methyLumiM), min(subset, nrow(methyLumiM)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(methyLumiM)]
		}
		methyLumiM <- methyLumiM[index,]
	} 
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
		return(methyLumiM)
	}
	if (removeGenderProbes && !is.null(annotation$CHR)) {
		redInd <- annotation$COLOR_CHANNEL == 'Red' & !(annotation$CHR %in% c('X', 'Y'))
		grnInd <- annotation$COLOR_CHANNEL == 'Grn' & !(annotation$CHR %in% c('X', 'Y'))
	} else {
		redInd <- annotation$COLOR_CHANNEL == 'Red'
		grnInd <- annotation$COLOR_CHANNEL == 'Grn'
	}
	
	nSample <- ncol(unmethy)
	density.list <- lapply(1:nSample, function(i) {
		red.a.i <- unmethy[redInd, i]
		red.b.i <- methy[redInd, i]
		grn.a.i <- unmethy[grnInd, i]
		grn.b.i <- methy[grnInd, i]
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
		if (colorMode) {
			if (logMode) {
				red[red < 1] <- 1
				grn[grn < 1] <- 1
				if (length(densityPar) > 0) {
					dd.red <- do.call('density', c(list(log2(red)), densityPar))
					dd.grn <- do.call('density', c(list(log2(grn)), densityPar))
				} else {
					dd.red <- density(log2(red))
					dd.grn <- density(log2(grn))
				}
			} else {
				if (length(densityPar) > 0) {
					dd.red <- do.call('density', c(list(red), densityPar))
					dd.grn <- do.call('density', c(list(grn), densityPar))
				} else {
					dd.red <- density(red)
					dd.grn <- density(grn)
				}
			}
		} else {
			pool <- c(red, grn)
			if (logMode) {
				pool[pool < 1] <- 1
				if (length(densityPar) > 0) {
					dd.pool <- do.call('density', c(list(log2(pool)), densityPar))
				} else {
					dd.pool <- density(log2(pool))
				}
			} else {
				if (length(densityPar) > 0) {
					dd.pool <- do.call('density', c(list(pool), densityPar))
				} else {
					dd.pool <- density(pool)
				}
			}
			dd.red <- dd.grn <- dd.pool
		}
		return(list(red=dd.red, green=dd.grn))
	})
	
	mm.density <- max(sapply(density.list, function(x) max(c(x$red$y, x$green$y))))
	xrange <- range(sapply(density.list, function(x) range(c(x$red$x, x$green$x))))

	info <- switch(channel, 
		"both"="intensity of both methylated and unmethylated probes",
		"methy"="intensity of methylated probes only",
		"unmethy"="intensity of unmethylated probes only",
		"sum"="CpG-site Intensity")
	if (logMode) {
		xlab <- paste("Log2", info)
	} else {
		substr(info, 1, 1) <- toupper(substr(info, 1, 1))
		xlab <- info
	}
	
	if (is.null(otherPar$xlab)) otherPar$xlab <- xlab
	if (is.null(otherPar$ylab)) otherPar$ylab <- "Density"
	if (is.null(otherPar$xlim)) otherPar$xlim <- xrange
	if (is.null(otherPar$ylim)) otherPar$ylim <- c(0,mm.density)
	if (is.null(otherPar$main)) otherPar$main <- ifelse(colorMode, "Compare density distribution of two color channels", "Compare density distribution")
	if (is.null(otherPar$type)) otherPar$type <- 'l'
	if (colorMode) {
		if (is.null(otherPar$col)) otherPar$col <- 2
	} else {
		if (is.null(otherPar$col)) {
			otherPar$col <- 1:nSample
		} else if (length(otherPar$col) < nSample)
			otherPar$col <- rep(otherPar$col[1], nSample)
	}
	
	
	for (i in 1:nSample) {
		dd.red.i <- density.list[[i]]$red
		dd.grn.i <- density.list[[i]]$green
		if (i == 1) {
			do.call('plot', c(list(dd.red.i), otherPar))
			if (colorMode) lines(dd.grn.i, col=3)
		} else {
			if (colorMode) {
				lines(dd.red.i, col=otherPar$col, lty=i)
				lines(dd.grn.i, col=3, lty=i)
			} else {
				lines(dd.red.i, col=otherPar$col[i], lty=i)
			}
		}
	}
	return(invisible(TRUE))
}


## plotColorBias2D(methyLumiM, selSample=1)
plotColorBias2D <- function(methyLumiM, selSample=1, combineMode=F, layoutRatioWidth=c(0.75,0.25), layoutRatioHeight=c(0.25, 0.75), 
			margins = c(5, 5, 2, 2), cex=1.25, logMode=TRUE, subset=NULL, ...) {

	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	otherPar <- list(...)
	
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(methyLumiM), min(subset, nrow(methyLumiM)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(methyLumiM)]
		}
		methyLumiM <- methyLumiM[index,]
	} 

	ff <- pData(featureData(methyLumiM))
	color.channel <- ff[,"COLOR_CHANNEL"]
	if (is.null(color.channel) && !combineMode) stop("No color channel information included in the data!\n Please add it using addAnnotationInfo function.\n")
	
	unmethy <- assayDataElement(methyLumiM, 'unmethylated')[, selSample[1]]
	methy <- assayDataElement(methyLumiM, 'methylated')[, selSample[1]]

	if (logMode) {
		unmethy[unmethy < 1] <- 1
		methy[methy < 1] <- 1
		unmethy <- log2(unmethy)
		methy <- log2(methy)
	}
	
	if (!combineMode) {
		grn.a <- unmethy[color.channel == "Grn"]
		red.a <- unmethy[color.channel == "Red"]
		grn.b <- methy[color.channel == "Grn"]
		red.b <- methy[color.channel == "Red"]
	}

	oldpar <- par(no.readonly = TRUE)
    layout(matrix(c(2,1,0,3), nrow=2), widths = layoutRatioWidth, heights = layoutRatioHeight, respect = FALSE)
	# layout.show(3)
	## plot the scatter plot 
	par(mar = c(margins[1], margins[2], 0, 0))
	plottype <- ifelse(combineMode, "p", "n")

	if (is.null(otherPar$pch)) otherPar$pch <- '.'
	if (is.null(otherPar$type)) otherPar$type <- plottype
	if (is.null(otherPar$xlim)) otherPar$xlim <- range(unmethy)
	if (is.null(otherPar$ylim)) otherPar$ylim <- range(methy)
	if (is.null(otherPar$xlab)) otherPar$xlab <- 'Unmethylated Probe Intensity'
	if (is.null(otherPar$ylab)) otherPar$ylab <- 'Methylated Probe Intensity'
	if (is.null(otherPar$main)) otherPar$main <- ''

	if (logMode) {
		if (is.null(otherPar$xlab)) otherPar$xlab <- 'Unmethylated Probe Intensity (log2)'
		if (is.null(otherPar$ylab)) otherPar$ylab <- 'Methylated Probe Intensity (log2)'
	} 
	do.call('plot', c(list(unmethy), list(methy), otherPar))

	if (combineMode) {
		dd.a <- density(unmethy)
		dd.b <- density(methy)
		par(mar = c(1, margins[2], margins[3], 0))
		plot(dd.a, xlab='', ylab='Density', xlim=otherPar$xlim, xaxt='n', col='black', type='l', main='')
		## plot the density plot of methylated probes
		par(mar = c(margins[1], 1, 0, margins[4]))
		plot(dd.b$y, dd.b$x, xlab='Density', ylab='', ylim=otherPar$ylim, yaxt='n', col='black', type='l', main='')
	} else {
		points(red.a, red.b, pch='.', cex=cex, col='red')
		points(grn.a, grn.b, pch='.', cex=cex, col='green')		

		dd.grn.a <- density(grn.a)
		dd.grn.b <- density(grn.b)
		dd.red.a <- density(red.a)
		dd.red.b <- density(red.b)
		## plot the density plot of unmethylated probes
		par(mar = c(1, margins[2], margins[3], 0))
		plot(dd.red.a, xlab='', ylab='Density', xlim=otherPar$xlim, ylim=range(c(dd.grn.a$y, dd.red.a$y)), xaxt='n', col='red', type='l', main='')
		lines(dd.grn.a, col='green')
		## plot the density plot of methylated probes
		par(mar = c(margins[1], 1, 0, margins[4]))
		plot(dd.red.b$y, dd.red.b$x, xlab='Density', ylab='', xlim=range(c(dd.grn.b$y, dd.red.b$y)), ylim=otherPar$ylim, yaxt='n', col='red', type='l', main='')
		lines(dd.grn.b$y, dd.grn.b$x, col='green')
	}

	par(oldpar)

	return(invisible(TRUE))
}


colorBiasSummary <- function(methyLumiM, logMode=TRUE, channel=c('both', 'unmethy', 'methy', 'sum')) {

	channel <- match.arg(channel)
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	annotation <- pData(featureData(methyLumiM))
	if (is.null(annotation$COLOR_CHANNEL)) {
		cat("No color balance adjustment because lack of COLOR_CHANNEL information!\n Please add it using addAnnotationInfo function.\n")
		return(methyLumiM)
	}
	redInd <- annotation$COLOR_CHANNEL == 'Red'
	grnInd <- annotation$COLOR_CHANNEL == 'Grn'
	
	nSample <- ncol(unmethy)
	summary.list <- lapply(1:nSample, function(i) {
		red.a.i <- unmethy[redInd, i]
		red.b.i <- methy[redInd, i]
		grn.a.i <- unmethy[grnInd, i]
		grn.b.i <- methy[grnInd, i]
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
	colnames(red.summary.matrix) <- colnames(grn.summary.matrix) <- sampleNames(methyLumiM)
	
	return(list(red=red.summary.matrix, green=grn.summary.matrix))
}


plotDensity <- function(dataMatrix, logMode=TRUE, addLegend=TRUE, legendPos="topright", subset=NULL, ...) {
	otherPar <- list(...)
	if (is(dataMatrix, 'MethyLumiM')) logMode <- FALSE

	if (is(dataMatrix, 'ExpressionSet')) {
	    dataMatrix <- exprs(dataMatrix)
	} else if (is.numeric(dataMatrix)) {
		dataMatrix <- as.matrix(dataMatrix)
	} else {
		stop('Un-supported class of dataMatrix.')
	}

	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(dataMatrix), min(subset, nrow(dataMatrix)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(dataMatrix)]
		}
		dataMatrix <- dataMatrix[index,,drop=FALSE]
	} 
	
	if (logMode) {
		if (any(dataMatrix <= 0)) dataMatrix[dataMatrix <= 0.1] <- 0.1 
		tmp <- apply(log2(dataMatrix), 2, density)
	} else {
		tmp <- apply(dataMatrix, 2, density)
	}
	xx <- sapply(tmp, function(x) x$x)
	yy <- sapply(tmp, function(x) x$y)

	if (is.null(otherPar$type)) otherPar$type <- 'l'
	if (is.null(otherPar$xlab)) otherPar$xlab <- "Intensity"
	if (is.null(otherPar$ylab)) otherPar$ylab <- 'Density'
	if (is.null(otherPar$lty)) otherPar$lty <-  1:ncol(dataMatrix)
	if (is.null(otherPar$col)) otherPar$col <-  1:ncol(dataMatrix)
	if (is.null(otherPar$lwd)) otherPar$lwd <-  1
	if (is.null(otherPar$main)) otherPar$main <- "Density plot"
	if (logMode)  otherPar$xlab <- paste(otherPar$xlab, " (log2)", sep="")

	do.call('matplot', c(list(xx), list(yy), otherPar))

	## add legend
	if (addLegend) {
		labels <- colnames(dataMatrix)
		if (is.null(labels)) labels <- as.character(1:ncol(dataMatrix))

		if (length(legendPos) > 1) {
			x.pos <- legendPos[1]
			y.pos <- legendPos[2]
		} else {
			x.pos <- legendPos[1]
			y.pos <- NULL
		}
		legend(x.pos, y.pos, legend=labels, box.lwd=0, col=otherPar$col, lty=otherPar$lty, lwd=otherPar$lwd)
	}
	return(invisible(TRUE))
}



produceMethylationGEOSubmissionFile <- function(methyLumiM, methyLumiM.raw=NULL, lib.mapping=NULL, idType='Probe', sampleInfo=NULL, fileName='GEOSubmissionFile.txt', supplementaryRdata=TRUE, ...) {
	
	if (missing(methyLumiM)) stop('Please provide all required input parameters!\n')
	if (is(methyLumiM, "MethyLumiM")) expr.norm <- estimateBeta(methyLumiM)
	if (is.null(methyLumiM.raw)) {
		detect <- detection(methyLumiM)
		methyData <- methylated(methyLumiM)
		unmethyData <- unmethylated(methyLumiM)
		expr <- NULL
	} else {
		detect <- detection(methyLumiM.raw)
		methyData <- methylated(methyLumiM.raw)
		unmethyData <- unmethylated(methyLumiM.raw)
		expr <- estimateBeta(methyLumiM.raw)
	}
	
	if (is.null(sampleInfo)) {
		sampleInfo <- produceGEOSampleInfoTemplate(methyLumiM, lib.mapping=lib.mapping, fileName=NULL)
	} else if (length(sampleInfo) == 1 && is.character(sampleInfo)) {
		sampleInfo <- read.table(sampleInfo, sep='\t', colClasses='character', skip=1, head=TRUE, strip.white=TRUE, quote='')
	} else if (is.null(nrow(sampleInfo))) {
		stop('Please provide correct sample information (a data.frame, matrix, or sampleInfo file)!\n')
	}
	sampleInfoTitle <- colnames(sampleInfo)
	if (any(sapply(sampleInfo[,-1, drop=F], nchar) == 0)) stop('No blank fields are allowed in the sampleInfo table!\nYou can check some example submissions, like GSM296418, at the GEO website.\n')
	if (supplementaryRdata) sampleInfo[, "Sample_supplementary_file"] <- 'supplementaryData.Rdata'
	nuID <- featureNames(methyLumiM)
	probeId <- nuID
	if (length(which(is.nuID(sample(nuID, 100)))) < 20) {
		nuID <- NULL
	} else {
		if (!is.null(lib.mapping)) {
			probeId <- nuID2IlluminaID(nuID, lib=lib.mapping, idType=idType, ...)
		} else {
			nuID <- NULL
		}
	}
	
	sampleID <- sampleInfo[, "sampleID"]
	sampleTitle <- sampleInfo[,'Sample_title']
	for (i in seq(sampleID)) {
		if (i == 1) {
			cat('^SAMPLE =', sampleTitle[i], '\n', sep='', file=fileName, append=FALSE)
		} else {
			cat('^SAMPLE =', sampleTitle[i], '\n', sep='', file=fileName, append=TRUE)			
		}
		sampleInfo.i <- paste('!', sampleInfoTitle[-1], ' = ', sampleInfo[i,-1], '\n', sep='', collapse='')
		sampleInfo.i <- gsub("'", "\\'", sampleInfo.i)
		cat(sampleInfo.i, file=fileName, append=TRUE, sep='')
		tableHead <- "ID_REF"
		cat("#ID_REF = Illumina ID\n", file=fileName, append=TRUE)
		if (!is.null(nuID)) {
			cat("#nuID = nucleotide universal IDentifier (nuID), convertible to and from probe sequence. See Bioconductor lumi package for more details.\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "nuID")
		}
		cat("#VALUE = Beta-value\n", file=fileName, append=TRUE)
		if (!is.null(expr)) cat("#RAW_VALUE = raw Beta-value\n", file=fileName, append=TRUE)
		tableHead <- c(tableHead, "VALUE")
		if (!is.null(expr)) tableHead <- c(tableHead, "RAW_VALUE")
		if (!is.null(methyData)) {
			cat("#METHYLATED = the intensities measured by methylated probes\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "METHYLATED")
		}
		if (!is.null(unmethyData)) {
			cat("#UNMETHYLATED = the intensities measured by unmethylated probes\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "UNMETHYLATED")
		}
		if (!is.null(detect)) {
			cat("#Detection_Pval = the detection p-value of the probe\n", file=fileName, append=TRUE)
			tableHead <- c(tableHead, "Detection_Pval")
		}
		sampleTable.i <- probeId
		if (!is.null(nuID)) sampleTable.i <- cbind(sampleTable.i, nuID)
		sampleTable.i <- cbind(sampleTable.i, expr.norm[,sampleID[i]])
		if (!is.null(expr)) sampleTable.i <- cbind(sampleTable.i, expr[,sampleID[i]])
		if (!is.null(methyData)) sampleTable.i <- cbind(sampleTable.i, methyData[,sampleID[i]])
		if (!is.null(unmethyData)) sampleTable.i <- cbind(sampleTable.i, unmethyData[,sampleID[i]])
		if (!is.null(detect)) sampleTable.i <- cbind(sampleTable.i, detect[,sampleID[i]])
		sampleTable.i <- rbind(tableHead, sampleTable.i)
		cat("!sample_table_begin\n", file=fileName, append=TRUE)
		write.table(sampleTable.i, sep='\t', quote=FALSE, file=fileName, append=TRUE, col.names=FALSE, row.names=FALSE)
		cat("!sample_table_end\n", file=fileName, append=TRUE)
	}
	
	if (supplementaryRdata) {
		methyLumiM <- methyLumiM[,sampleID]
		if (!is.null(methyLumiM.raw)) {
			methyLumiM.raw <- methyLumiM.raw[,sampleID]
			save(methyLumiM, methyLumiM.raw, sampleInfo, file='supplementaryData.Rdata')
		} else {
			save(methyLumiM, sampleInfo, file='supplementaryData.Rdata')
		}
	}
}


## ---------------------------------------------------------------------------------
# Functions related with estimating methylation status
## EM estimation of the parameters s1, s2 and theta
## E-step: determine the class of each probe
## M-step: 
## 1. estimation the proportion of two classes
## 2. Estimate of optimized s1 and s2 using optimization method
## 3. Estimate theta based on equation

# fittedGamma <- gammaFitEM(M[,1], initialFit=NULL, maxIteration=50, tol=0.0001, plotMode=T, verbose=T)
gammaFitEM <- function(M, initialFit=NULL, fix.k=NULL, weighted=TRUE, maxIteration=50, tol=0.0001, plotMode=FALSE, verbose=FALSE) {
	
	fix.theta=NULL
	eps <- 10^-5
	# if (!require(nleqslv)) 	stop("Please install nleqslv package!\n")
	if (!is.null(fix.k)) {
		if (length(fix.k) == 1) fix.k <- c(fix.k, fix.k)
	}
	if (!is.null(fix.theta)) {
		if (length(fix.theta) == 1) fix.theta <- c(fix.theta, fix.theta)
	}
		
	fs1 <- function(s1, mi, zi, theta, k, n) {
		nonNegativeInd <- which(mi > s1) # this setting will allow to shift the profile without considering the small part of the trailing points
		sum(zi[nonNegativeInd]/(mi[nonNegativeInd]-s1)) - n /((k-1)*theta)	
		# mi[mi < s1] <- s1 + eps	
		# sum(zi/(mi-s1)) - n /((k-1)*theta)	
		#sum(zi/abs(mi-s1)) - n /((k-1)*theta)		
	}

	fs2 <- function(s2, mi,  zi, theta, k, n) {
		nonNegativeInd <- which(mi < s2)  # this setting will allow to shift the profile without considering the small part of the trailing points
		sum(zi[nonNegativeInd]/(s2-mi[nonNegativeInd])) - n /((k-1)*theta)	
		#mi[mi > s2] <- s2 - eps
		#sum(zi/(s2-mi)) - n /((k-1)*theta)	
		# sum(zi/abs(s2-mi)) - n /((k-1)*theta)	
	}

	fk1 <- function(k1, mi, zi, s, n) {
		nonNegativeInd <- which(mi > s) # this setting will allow to shift the profile without considering the small part of the trailing points
		log(k1) - digamma(k1) - log(sum(zi[nonNegativeInd]*(mi[nonNegativeInd]-s)/n)) + sum(zi[nonNegativeInd]*log(mi[nonNegativeInd]-s)/n)
		# mi[mi < s] <- s + eps
		# log(k1) - digamma(k1) - log(sum(zi*(mi-s)/n)) + sum(zi*log(mi-s)/n)
	}

	fk2 <- function(k2, mi, zi, s, n) {
		nonNegativeInd <- which(mi < s) # this setting will allow to shift the profile without considering the small part of the trailing points
		log(k2) - digamma(k2) - log(sum(zi[nonNegativeInd]*(s-mi[nonNegativeInd])/n)) + sum(zi[nonNegativeInd]*log(s-mi[nonNegativeInd])/n)
		# mi[mi > s] <- s - eps
		# log(k2) - digamma(k2) - log(sum(zi*(s-mi)/n)) + sum(zi*log(s-mi)/n)
	}
			
	if (is.matrix(M)) {
		if (ncol(M) > 1) cat("Only the first column of the matrix was processed!\n")
		x <- M[,1]
	} else {
		x <- M
	}

	# Initial value estimation got from grid search
	initialFit <- .initialGammaEstimation(x, initialFit=initialFit)
	
	k <- initialFit$k
	theta <- initialFit$theta
	s <- initialFit$shift
	p <- initialFit$proportion
	Mode <- initialFit$mode
	k[k < 2] <- 2
	theta[theta < 0.2] <- 0.2

	if (!is.null(fix.k)) k[!is.na(fix.k)] <- fix.k[!is.na(fix.k)]
	if (!is.null(fix.theta)) theta[!is.na(fix.theta)] <- fix.theta[!is.na(fix.theta)]
	if (!is.null(fix.k) || !is.null(fix.theta)) s <- c(Mode[1] - (k[1]-1)*theta[1], Mode[2] + (k[2]-1)*theta[2])

	if (verbose) {
		if (plotMode) plotGammaFit(x, k=k, theta=theta, shift=s, proportion=p)
		cat("\nInitial estimation:\n")
		cat("k:", k, "\n")
		cat("s:", s, "\n")
		cat("theta:", theta, "\n")
		cat("p:", p, "\n")
		f1 <- dgamma(x-s[1], shape=k[1], scale=theta[1])
		f2 <- dgamma(s[2]-x, shape=k[2], scale=theta[2])
		logLikelihood <- sum(log(p[1] * f1 + p[2]*f2))
		cat("logLikelihood:", logLikelihood, "\n")
	}
	
	iter <- 1
	N <- length(x)
	while (iter < maxIteration) {
		if (any(theta < eps) || any(k < 2)) {
			cat('It does not converge based on the current initial values!\n')
			logLikelihood <- -Inf
			return(list(logLikelihood=logLikelihood, k=k, theta=theta, shift=s, proportion=p, mode=Mode, probability=NULL))
		}
		
		Mode <- c(s[1] + (k[1]-1)*theta[1], s[2] - (k[2]-1)*theta[2])

		## E-step
		f1 <- dgamma(x-s[1], shape=k[1], scale=theta[1])
		f2 <- dgamma(s[2]-x, shape=k[2], scale=theta[2])
		## TO DO: add weighted function instead of truncate 10/25/2010
		
		if (weighted) {
			# down weight the long tails of two component densities beyond their modes
			w1 <- rep(1, length(f1)); w2 <- rep(1, length(f2))
			w1[x > Mode[2]] <- f2[x > Mode[2]] / max(f2)
			w2[x < Mode[1]] <- f1[x < Mode[1]] / max(f1)
			f1 <- f1 * w1
			f2 <- f2 * w2
			
			# adjust the density values because of the trunction
			# f1 <- f1 / (pgamma(Mode[2]-s[1], shape=k[1], scale=theta[1], lower.tail = TRUE))
			# f2 <- f2 / (pgamma(s[2]-Mode[1], shape=k[2], scale=theta[2], lower.tail = TRUE))
			# z1 <- p[1] * f1 / (p[1] * f1 + p[2] * f2) # + eps * 10)
			# z1[x > Mode[2]] <- 0
			# z1[x < Mode[1]] <- 1
			# f1[x > Mode[2]] <- 0
			# f2[x < Mode[1]] <- 0
		} 
		z1 <- p[1] * f1 / (p[1] * f1 + p[2] * f2) # + eps * 10)  # posterior probability of unmethylated
		z2 <- 1 - z1

		if (verbose && iter > 1) {
			# logLikelihood <- sum((z1*log(p[1] * f1))[z1 > 0]) + sum((z2* log(p[2]*f2))[z2 > 0])
			ll <- p[1] * f1 + p[2]*f2
			ll[ll < eps] <- eps
			logLikelihood <- sum(log(ll))		
			cat("logLikelihood:", logLikelihood, "\n")			
		}
		
		## M-step (estimate parameters: s, theta and k)
		# update the proportion of each class
		n1 <- sum(z1); n2 <- N - n1
		p.new <- n1 / N
		p.new <- c(p.new, 1 - p.new)
    	if (abs(p.new[1]-p[1]) < tol) break
		p <- p.new

		## Estimation of s, theta and k
		iter.inter <- 1
		while (iter.inter <= 10) {
			# set the range of s parameter		fss(s1, mi, zi, theta, k, n)
			s1.new <- nleqslv(s[1], fs1, mi=x, zi=z1, theta=theta[1], k=k[1], n=n1)$x
			s2.new <- nleqslv(s[2], fs2, mi=x, zi=z2, theta=theta[2], k=k[2], n=n2)$x
			s.new <- c(s1.new, s2.new)

			## Estimation of k # fk1(k1, mi, zi, s, n)
			if (is.null(fix.k)) {
				k1.new <- nleqslv(k[1], fk1, mi=x, zi=z1, s=s.new[1], n=n1)$x
				k2.new <- nleqslv(k[2], fk2, mi=x, zi=z2, s=s.new[2], n=n2)$x
				k.new <- c(k1.new, k2.new)
			} else {
				k.new <- fix.k
				if (any(is.na(k.new))) {
					if (is.na(fix.k[1])) k.new[1] <- nleqslv(k[1], fk1, mi=x, zi=z1, s=s.new[1], n=n1)$x
					if (is.na(fix.k[2])) k.new[2] <- nleqslv(k[2], fk2, mi=x, zi=z2, s=s.new[2], n=n2)$x
				}
			}

			## Estimation of theta
			if (is.null(fix.theta)) {
				theta.new <- c(sum(z1*(x - s.new[1]))/(k.new[1]*n1), sum(z2*(s.new[2] - x))/(k.new[2]*n2))
			} else {
				theta.new <- fix.theta
				if (any(is.na(theta.new))) {
					if (is.na(fix.theta[1])) theta.new[1] <- sum(z1*(x - s.new[1]))/(k.new[1]*n1)
					if (is.na(fix.theta[2])) theta.new[2] <- sum(z2*(s.new[2] - x))/(k.new[2]*n2)
				}
			}

			if (min(abs(s.new - s)) < tol && min(abs(k.new - k)) < tol && min(abs(theta.new - theta)) < tol)  break
			# if (verbose) cat("s:", s, " k:", k, " theta:", theta, "\n")
			s <- s.new; k <- k.new; theta <- theta.new
			# if (any(theta < 0.2) || any(k < 2)) {
			if (any(theta < eps) || any(k < 2)) {
				cat('It does not converge based on the current initial values!\n')
				logLikelihood <- -Inf
				return(list(logLikelihood=logLikelihood, k=k, theta=theta, shift=s, proportion=p, mode=Mode, probability=NULL))
			}
			iter.inter <- iter.inter + 1
		}
 
		if (verbose) {
			cat("\nIteration ", iter, "\n")
			cat("k:", k, "\n")
			cat("s:", s, "\n")
			cat("theta:", theta, "\n")
			cat("p:", p, "\n")
			cat("Mode:", Mode, "\n")
			if (plotMode) plotGammaFit(x, k=k, theta=theta, shift=s, proportion=p)
		}
		iter <- iter + 1
	}
	if (iter == maxIteration) {
		Mode <- c(s[1] + (k[1]-1)*theta[1], s[2] - (k[2]-1)*theta[2])
		f1 <- dgamma(x-s[1], shape=k[1], scale=theta[1])
		f2 <- dgamma(s[2]-x, shape=k[2], scale=theta[2])
		if (truncate) {
			f1 <- f1 / (pgamma(Mode[2]-s[1], shape=k[1], scale=theta[1], lower.tail = TRUE))
			f2 <- f2 / (pgamma(s[2]-Mode[1], shape=k[2], scale=theta[2], lower.tail = TRUE))
			z1 <- p[1] * f1 / (p[1] * f1 + p[2] * f2) # + eps)
			z1[x > Mode[2]] <- 0
			z1[x < Mode[1]] <- 1
			f1[x > Mode[2]] <- 0
			f2[x < Mode[1]] <- 0
		} else {
			z1 <- p[1] * f1 / (p[1] * f1 + p[2] * f2) # + eps)  # posterior probability of unmethylated
		}
		z2 <- 1 - z1
	}

	# estimate log-likelihood 
	ll <- p[1] * f1 + p[2]*f2
	ll[ll < eps] <- eps
	logLikelihood <- sum(log(ll))		
	if (iter == maxIteration && verbose) cat("logLikelihood:", logLikelihood, "\n")

	# return the methylation/unmethylation estimation of each probe
	## check class index to make sure those extreme probe belonging to the corresponding group
	probability <- cbind(z1, z2)
	colnames(probability) <- c('unmethylated', 'methylated')
	Mode <- c(s[1] + (k[1]-1)*theta[1], s[2] - (k[2]-1)*theta[2])
	if (plotMode) plotGammaFit(x, k=k, theta=theta, shift=s, proportion=p)

	fitResult <- list(logLikelihood=logLikelihood, k=k, theta=theta, shift=s, proportion=p, mode=Mode, probability=probability)
	class(fitResult) <- 'gammaFit'

	return(fitResult)
}


# initial gamma parameters estimation
.initialGammaEstimation <- function(x, initialFit=NULL) {

	k <- theta <- s <- p <- Mode <- NULL
	if (!is.null(initialFit)) {
		k <- initialFit$k
		theta <- initialFit$theta
		s <- initialFit$shift
		p <- initialFit$proportion
		Mode <- initialFit$mode
	} 
	
	if (!is.null(k) && !is.null(theta) &&!is.null(s) && !is.null(p) && !is.null(Mode)) 	return(initialFit)
	
	# mode positions
	if (is.null(Mode)) {
		dd <- density(x)
		density.m <- dd$y
		density.x <- dd$x  # x-axis of density plot (methylation ratio)
		midpoint <- mean(quantile(x, c(0.01, 0.99)))
		unmethy.ind <- density.x <= midpoint
		P1 <- density.x[unmethy.ind][which.max(density.m[unmethy.ind])]
		if (abs(P1 - midpoint) < 0.1) P1 <- (midpoint + quantile(x, 0.01))/2
		methy.ind <- density.x >= midpoint
		P2 <- density.x[methy.ind][which.max(density.m[methy.ind])]		
		if (abs(P2 - midpoint) < 0.1) P2 <- (midpoint + quantile(x, 0.99))/2
		Mode <- c(P1, P2)
	}

	# update the class index
	unmethy.ind <- which(x < (Mode[1] + Mode[2])/2)
	methy.ind <- which(x > (Mode[1] + Mode[2])/2)

	# percentage of two classes
	if (is.null(p)) {
		p <- c(length(unmethy.ind)/length(x), length(methy.ind)/length(x))
	}

	mean.unmethy <- mean(x[unmethy.ind])
	mean.methy <- mean(x[methy.ind])
	var.unmethy <- var(x[unmethy.ind])
	var.methy <- var(x[methy.ind])
	if (is.null(theta)) {
		theta <- c(mean.unmethy - Mode[1], Mode[2] - mean.methy)
		theta[theta < 0.2] <- 0.2
	}
	if (is.null(k))  {
		k <- c(round(var.unmethy/theta[1]^2), round(var.methy/theta[2]^2))
		k[k < 2] <- 2
	}
	if (is.null(s))  s <- c(mean.unmethy - k[1]*theta[1], mean.methy + k[2]*theta[2])

	fitResult <- list(k=k, theta=theta, shift=s, proportion=p, mode=Mode)
	class(fitResult) <- 'gammaFit'
	return(fitResult)
}


# plot gammFit results
plotGammaFit <- function(x, gammaFit=NULL, k=NULL, theta=NULL, shift=NULL, proportion=NULL, plotType=c('histogram', 'density'), ...) {
	plotType <- match.arg(plotType)
	
	if (!is.null(gammaFit)) {
		if (class(gammaFit) != 'gammaFit') stop("gammaFit should be an object of 'gammaFit' class!")
		k <- gammaFit$k
		theta <- gammaFit$theta
		shift <- gammaFit$shift
		proportion <- gammaFit$proportion
	}
	x <- sort(x)

	if (is.null(k) || is.null(theta) || is.null(shift) || is.null(proportion)) stop("Information of parameters k, theta, shift and proportion is required!\n")
	y1 = dgamma(x-shift[1], shape=k[1], scale=theta[1]) 
	y2 = dgamma(shift[2]-x, shape=k[2], scale=theta[2])
	if (plotType == 'histogram') {
		hist(x, 50, probability=TRUE, main='Compare data histogram and fitted distribution', xlab='M value', ...)
	} else {
		plot(density(x), type='l', lwd=1.5, main='Compare data density and fitted distribution', xlab="M value", ...)
	}
	lines(x, y1*proportion[1] + y2 * proportion[2], col=2, lty=3, lwd=1.5)
	lines(x, y1*proportion[1], col=3, lty=2)
	lines(x, y2*proportion[2], col=4, lty=2)
	return(invisible(TRUE))
}


# estimate methylation call probability based on gamma fit parameters
methylationCall <- function(x, threshold=0.95, ...) {
		
	if (length(threshold) == 1) threshold <- rep(threshold, 2)
	probability <- NULL
	if (class(x) != 'gammaFit') {
		fit <- gammaFitEM(x, ...) 
	} else {
		fit <- x
	}
	probability <- fit$probability
	
	# methyCall <- probability > threshold
	methyCall <- apply(probability, 1, function(x) {
		cc <- "Margin"
		if (x[1] > threshold[1]) {
			cc <- "Unmethy"
		} else if (x[2] > threshold[2]) {
			cc <- "Methy"
		}
		return(cc)
	})

	attr(methyCall, "probability") <- probability[,"methylated"]
	return(methyCall)	
}


# methylation status
lumiMethyStatus <- function(methyLumiM, ...) 
{

	if (!is(methyLumiM, 'MethyLumiM')) {
		stop('The object should be class "MethyLumiM" inherited!')
	}
	
	history.submitted <- as.character(Sys.time())
	
	M <- exprs(methyLumiM)
	M.status <- M.prob <- NULL
	for (i in 1:ncol(M)) {
		status.i <- methylationCall(M[,1], ...)
		prob.i <- attr(status.i, "probability")
		M.status <- cbind(M.status, status.i)
		M.prob <- cbind(M.prob, prob.i)
	}
	rownames(M.status) <- rownames(M.prob) <- rownames(M)
	colnames(M.status) <- colnames(M.prob) <- colnames(M)

	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiMethyStatus)))
	
	lumiVersion <- packageDescription('lumi')$Version
	attr(M.status, "history") <- rbind(methyLumiM@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	
	attr(M.status, "probability") <- M.prob	
	return(M.status)
}

