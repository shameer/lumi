##=================================================
## Define LumiBatch object:

setClass('LumiBatch', 
	representation(history='data.frame', controlData='data.frame', QC='list'), 
	prototype=list(history=data.frame(
		submitted   = I(vector()),
		finished    = I(vector()),
		command     = I(vector()),
		lumiVersion = I(vector())
	), controlData = data.frame(), QC = list()),
	contains='ExpressionSet')

#
setMethod('initialize', 'LumiBatch', function(.Object, 
	exprs = new('matrix'),
	se.exprs = new('matrix'),		# standard deviation of the bead measurements of each probe
	detection = new('matrix'),
	beadNum = new('matrix'),
    ...,
    assayData)
{
	if (missing(assayData)) {
		cmd <- 'assayData <- assayDataNew(exprs=exprs, se.exprs=se.exprs'
		nSample <- ncol(exprs)
		if (ncol(detection) == nSample) cmd <- paste(cmd, ', detection=detection')
		if (ncol(beadNum) == nSample) cmd <- paste(cmd, ', beadNum=beadNum')
		cmd <- paste(cmd, ')')
		eval(parse(text=cmd))
	}
	else if (!missing(exprs) || !missing(se.exprs))
		stop("only one of 'assayData' or ('exprs' and 'se.exprs') allowed")
	callNextMethod(.Object, assayData=assayData, ...)
})


setValidity("LumiBatch", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "ExpressionSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("se.exprs")))
		## validate controlData slot, added by Mark Cowley
		msg <- Biobase:::validMsg(msg, valid_controlData(object))
    if (is.null(msg)) TRUE else msg
})


## validate controlData slot, added based on the code provided by Mark Cowley
valid_controlData <- function(object) {
	is(object, "LumiBatch") || return("object must be a LumiBatch")
	if ( nrow(controlData(object)) == 0 ) {
		return(TRUE)
	} else if ( ncol(controlData(object)) == (ncol(object) + 2) ) {
	 	if ( !identical(colnames(controlData(object))[1:2], c("controlType", "ProbeID")) ) {
		  cat("The first 2 extra columns in the controlData should be 'controlType' and 'ProbeID'.\n")
			return(FALSE)
		} 
		sampleName.ctrl <- colnames(controlData(object))[-c(1,2)]
	} else {
		sampleName.ctrl <- colnames(controlData(object))
	}

	if (!identical(sampleName.ctrl, sampleNames(object))) {
		cat("The sample names in the controlData don't match sampleNames(object).\n")
		return(FALSE)
	}
	
	return(TRUE)
}


##=================================================
## methods

setMethod("se.exprs", signature(object="ExpressionSet"), function(object) {
	if ('se.exprs' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"se.exprs"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("se.exprs", signature(object="ExpressionSet"), function(object, value) {
		if (is.null(value)) {
			assay <- assayData(object)
			if (exists('se.exprs', envir=assay)) {
				oldMode <- storageMode(assay)
				storageMode(assay) <- 'environment'
				rm(se.exprs, envir=assay)
				storageMode(assay) <- oldMode
				assayData(object) <- assay
			}
			return(object)
		} else {
			assayDataElementReplace(object, "se.exprs", value)
		}
	})


setMethod("beadNum", signature(object="ExpressionSet"), function(object) {
	if ('beadNum' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"beadNum"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("beadNum", signature(object="ExpressionSet"), function(object, value) {
		if (is.null(value)) {
			assay <- assayData(object)
			if (exists('beadNum', envir=assay)) {
				oldMode <- storageMode(assay)
				storageMode(assay) <- 'environment'
				rm(beadNum, envir=assay)
				storageMode(assay) <- oldMode
				assayData(object) <- assay
			}
			return(object)
		} else {
			assayDataElementReplace(object, "beadNum", value)
		}
	})

setMethod("detection", signature(object="ExpressionSet"), function(object) {
	if ('detection' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"detection"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("detection", signature(object="ExpressionSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('detection', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(detection, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "detection", value)
	}
})	

setMethod("controlData", signature(object="LumiBatch"), function(object) {
	object@controlData
})


## added based on the code provided by Mark Cowley
setReplaceMethod("controlData", signature(object="LumiBatch"), function(object, value) {
	
	if (is.null(value)) value <- data.frame()
	
	if (is.character(value) && length(value) == 1) {
		if ( file.exists(value) ) {
			object <- addControlData2lumi(value, object)
			validObject(object)
		}
	} else if (is(value, 'data.frame')) {
		object@controlData <- value
		validObject(object)
	} else {
		stop("value must be one of: NULL, data.frame, or character(1) as the path to a control data file")
	}
	
	return(object)
})	


setMethod("getHistory",signature(object="LumiBatch"), function(object) object@history)


setReplaceMethod("sampleNames", signature(object="LumiBatch",value="ANY"),
	function(object, value) {

	 	object <- callNextMethod()
		ddim <- dim(object)
        
	 	## subsetting the QC information
	 	if (!is.null(object@QC)) {
	 		QC <- object@QC
	 		if (!is.null(QC$sampleSummary))
	 			if (ncol(QC$sampleSummary) == ddim[2])
	 				colnames(QC$sampleSummary) <- value
	 		if (!is.null(QC$BeadStudioSummary))
	 			if (nrow(QC$BeadStudioSummary) == ddim[2])
	 				rownames(QC$BeadStudioSummary) <- value
	 		object@QC <- QC
	 	}
	 	if (!is.null(attr(object, 'vstParameter'))) {
	 		vstParameter <- attr(object, 'vstParameter')
	 		if (!is.null(nrow(vstParameter))) {
	 			if (nrow(vstParameter) == ddim[2]) {
					rownames(vstParameter) <- value
					transformFun <- attr(object, 'transformFun')
					names(transformFun) <- value
	 				attr(object, 'vstParameter') <- vstParameter
	 				attr(object, 'transformFun') <- transformFun
	 			}
	 		}
	 	}
        
	 	## controlData information
	 	if (nrow(object@controlData) > 0) {
			if (ncol(object@controlData) == ddim[2]) {
				colnames(object@controlData) <- value
			} else if (ncol(object@controlData) == ddim[2] + 2) {
				colnames(controlData(object)) <- c("controlType", "ProbeID", value)
			} 				
	 	}
	 	    
	 	validObject(object)
	  return(object)
	})


setMethod("summary",signature(object="LumiBatch"), function(object, type=c('data', 'QC')) 
{
	type <- match.arg(type)
	if (type == 'data') {
		show(object)		
	} else {
		if (is.null(object@QC$sampleSummary)) {
			cat('Run Quality Control estimation ...\n')
			object <- lumiQ(object)
		}
		if (!is.null(object@QC$sampleSummary)) {
			QC <- object@QC
			dimen <- dim(exprs(object))
			cat(paste("Data dimension: ", paste(dimen[1], 'genes', 'x', dimen[2], 'samples', collapse=" "), '\n'))
			sampleSummary <- QC$sampleSummary
			cat(paste("\nSummary of Samples:\n", sep=''))
			print(sampleSummary, quote=FALSE)

			cat('\nMajor Operation History:\n')
			print(QC$history, quote=FALSE) 
		}
	}
})


setMethod("show",signature(object="LumiBatch"), function(object) 
{
	cat('Summary of data information:\n')
	note <- notes(object)
	n.names <- names(note)
	for (i in seq(note)) {
		if (!is.null(n.names)) cat('\t', paste(n.names[i], ':\n\t\t', sep=''))
		cat(note[[i]], sep='\n\t\t')		
	}
	cat('\nMajor Operation History:\n')
	hh <- getHistory(object)
	if (nrow(hh) > 5) {
		print(head(hh, 2))
		cat("...\n")
		print(tail(hh, 2))
	} else {
		print(hh)
	}	
	cat('\nObject Information:\n')
	callNextMethod()
	if (nrow(controlData(object))>0) {
		cat('Control Data: Available\n')
	} else {
		cat('Control Data: N/A\n')
	}
	cat("QC information: Please run summary(x, 'QC') for details!\n")
})


setMethod("[", "LumiBatch", function(x, i, j, ..., drop = FALSE) 
{
	if (missing(drop)) drop <- FALSE
   	history.submitted <- as.character(Sys.time())
	ddim <- dim(x)
	
	sampleName <- sampleNames(x)
	## convert the names as index to avoid some potential problems of name inconsistency
	if (!missing(j)) {
		if (is.character(j)) {
			ind <- seq(sampleName)
			names(ind) <- sampleName
			j <- ind[j]
		}	
	}
	## do default processing of 'ExpressionSet'
	x <- callNextMethod()

	if (!missing(i) && !missing(j)) {
		history.command <- paste('Subsetting', ddim[1], 'features and', ddim[2], 'samples.')		
	} else if (!missing(i)) {
		history.command <- paste('Subsetting', ddim[1], 'features.')
	} else if (!missing(j)) {
		history.command <- paste('Subsetting', ddim[2], 'samples.')
	} else {
		return(x)
	}

	## subsetting the QC information
	if (!missing(j)) {
		if (!is.null(x@QC)) {
			QC <- x@QC
			if (!is.null(QC$sampleSummary)) {
				sampleSummary <- QC$sampleSummary
				if (ncol(sampleSummary) == ddim[2]) {
					QC$sampleSummary <- sampleSummary[,j,drop=FALSE]
				} else if (ncol(sampleSummary) > 0) {
					if (all(sampleName %in% colnames(sampleSummary))) {
						QC$sampleSummary <- sampleSummary[,sampleName[j],drop=FALSE]
					} else {
						warning('The sampleSummary in QC slot does not match the sampleNames.\nThe subsetting did not execute on sampleSummary.\n')
					}
				}
			}
			if (!is.null(QC$BeadStudioSummary)) {
				BeadStudioSummary <- QC$BeadStudioSummary
				if (nrow(BeadStudioSummary) == ddim[2]) {
					QC$BeadStudioSummary <- BeadStudioSummary[j,,drop=FALSE]
				} else if (nrow(BeadStudioSummary) > 0) {
					if (all(sampleName %in% rownames(BeadStudioSummary))) {
						QC$BeadStudioSummary <- BeadStudioSummary[sampleName[j],,drop=FALSE]
					} else {
						warning('The BeadStudioSummary in QC slot does not match the sampleNames.\nThe subsetting did not execute on BeadStudioSummary.\n')
					}
				}
			}
			x@QC <- QC
		}
		if (!is.null(attr(x, 'vstParameter'))) {
			vstParameter <- attr(x, 'vstParameter')
			if (!is.null(nrow(vstParameter))) {
				if (nrow(vstParameter) == ddim[2]) {
					attr(x, 'vstParameter') <- attr(x, 'vstParameter')[j,,drop=FALSE]
					attr(x, 'transformFun') <- attr(x, 'transformFun')[j]
				}
			}
		}

		## controlData information
		if (nrow(controlData(x)) > 0) {
			if (is.logical(j) || is.numeric(j)) j <- sampleName[j]
			if (all(j %in% colnames(controlData(x)))) {
				if (all(c("controlType", "ProbeID") %in% colnames(controlData(x)))) {
				  j <- c("controlType", "ProbeID", j)
				}
				controlData(x) <- controlData(x)[,j, drop=FALSE]
			} else {
				warning('The controlData slot does not match the sampleNames.\nThe subsetting did not execute on controlData.\n')
			}
		}
	}

  # history tracking
  history.finished <- as.character(Sys.time())
	if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
	lumiVersion <- packageDescription('lumi')$Version
	x@history<- rbind(x@history, data.frame(submitted=history.submitted, finished=history.finished, 
			command=history.command, lumiVersion=lumiVersion))

  validObject(x)
	return(x)
})


setMethod("combine", signature=c(x="LumiBatch", y="ExpressionSet"), function(x, y, ...) 
{
	if (missing(y)) return(x)
	warning('The Lumibatch object was forced as ExpressionSet.')
	x <- as(x, 'ExpressionSet')
	if (length(list(...)) > 0) 
		return(combine(x, combine(y, ...)))
	return(combine(x, y))
})

setMethod("combine", signature=c(x="ExpressionSet", y="LumiBatch"), function(x, y, ...) 
{
	if (missing(y)) return(x)
	warning('The LumiBatch object was forced as ExpressionSet.')
	y <- as(y, 'ExpressionSet')
	if (length(list(...)) > 0) 
		return(combine(x, combine(y, ...)))
	return(combine(x, y))
})

setMethod("combine", signature=c(x="LumiBatch", y="LumiBatch"), function(x, y, ...) 
{
	if (missing(y)) return(x)
	if (length(list(...)) > 0) 
		return(combine(x, combine(y, ...)))
	if (class(x) != class(y))
		stop(paste("objects must be the same class, but are ",
                 class(x), ", ", class(y), sep=""))
	
	if (any(sort(featureNames(x)) != sort(featureNames(y)))) stop('Two data sets have different row names.')
	## determine whether there are duplicated sample names
	sampleName.x <- sampleNames(x)
	sampleName.y <- sampleNames(y)
	if (any(sampleName.x %in% sampleName.y)) {
		warning('Two data sets have some duplicated sample names.\n "_2" were attached to the duplicated sample names.')
		sampleName.x <- sampleNames(x)   # paste(sampleNames(x), '_1', sep='')
		sampleName.y <- paste(sampleNames(y), '_2', sep='')
		sampleNames(x) <- sampleName.x
		sampleNames(y) <- sampleName.y
	}
	featureName.x <- featureNames(x)
	featureName.y <- featureNames(y)
	featureName.com <- featureName.x[featureName.x %in% featureName.y]
	if (length(featureName.com) < length(featureName.x) || length(featureName.com) != length(featureName.y)) {
		if (length(featureName.com) > 0) {
			warning('Two data sets have different featureNames, only the common ones were used.')			
		} else {
			stop('Two data sets have totally different featureNames.')
		}
	}
	## make sure two data sets have the sample order of features
	x <- x[featureName.com,]
	y <- y[featureName.com,]
	history.submitted <- as.character(Sys.time())
	dimm.x <- dim(x); dimm.y <- dim(y)
	assayData(x) <- combine(assayData(x), assayData(y))
	experimentData(x) <- combine(experimentData(x),experimentData(y))

	## combine pheno data
	phenoData(x) <- combine(phenoData(x),phenoData(y))
	## combine protocolData 
	protocolData(x) <- combine(protocolData(x),protocolData(y))
	
	## combining the QC information
	if (length(x@QC) > 0 && length(y@QC) > 0) {
		if (!is.null(x@QC$BeadStudioSummary) && !is.null(y@QC$BeadStudioSummary)) {
			if (ncol(x@QC$BeadStudioSummary) == ncol(y@QC$BeadStudioSummary) && ncol(x@QC$BeadStudioSummary) > 0) {
				BeadStudioSummary <- rbind(x@QC$BeadStudioSummary, y@QC$BeadStudioSummary)
				x@QC$BeadStudioSummary <- BeadStudioSummary
			} else {
				x@QC <- list()
			}
		} else {
			x@QC <- list()
		}
		if (!is.null(x@QC$sampleSummary) && !is.null(y@QC$sampleSummary)) {
			if (nrow(x@QC$sampleSummary) == nrow(y@QC$sampleSummary)) {
				sampleSummary <- cbind(x@QC$sampleSummary, y@QC$sampleSummary)
				x@QC$sampleSummary <- sampleSummary
			} else {
				x@QC <- list()
			}
		} else {
			x@QC <- list()
		}
		if (!is.null(x@QC)) {
			history.x <- x@QC$history
			if (is.null(history.x)) history.x <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
			if (is.null(history.x$lumiVersion)) history.x$lumiVersion <- rep(NA, nrow(history.x))
			history.y <- y@QC$history
			if (is.null(history.y)) history.y <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
			if (is.null(history.y$lumiVersion)) history.y$lumiVersion <- rep(NA, nrow(history.y))
			x@QC$history <- rbind(history.x, history.y)
		} 
	} else {
		x@QC <- list()
	}
	
	## VST transformation parameters
	if (!is.null(attr(x, 'vstParameter')) && !is.null(attr(y, 'vstParameter'))) {
		vstParameter.x <- attr(x, 'vstParameter')
		vstParameter.y <- attr(y, 'vstParameter')
		if (is.null(nrow(vstParameter.x))) {
			vstParameter.x <- matrix(vstParameter.x, nrow=1)
		}
		if (is.null(nrow(vstParameter.y))) {
			vstParameter.y <- matrix(vstParameter.y, nrow=1)
		}
		if (nrow(vstParameter.x) != dimm.x[2] || nrow(vstParameter.y) != dimm.y[2]) {
			attr(x, 'vstParameter') <- attr(x, 'transformFun') <- NULL
		} else {
			attr(x, 'vstParameter') <- rbind(attr(x, 'vstParameter'), attr(y, 'vstParameter'))
			attr(x, 'transformFun') <- c(attr(x, 'transformFun'), attr(y, 'transformFun'))
		}
	}
	
	## controlData information
	if (nrow(controlData(x)) > 0) {
		if (nrow(controlData(x)) == nrow(controlData(y))) {
			if (c('controlType', 'ProbeID') %in% colnames(controlData(x)) && 
				  c('controlType', 'ProbeID') %in% colnames(controlData(y))) {
					if (all(controlData(x)$controlType == controlData(y)$controlType) &&
				      all(controlData(x)$ProbeID == controlData(y)$ProbeID)) {
								controlData <- cbind(controlData(x), controlData(y)[,-c(1,2)])
								controlData(x) <- as.data.frame(controlData)
					} else {
						warning("controlData slot not combined: controlType and ProbeID do not match.")
					}
			} else {
				warning("controlData slot not combined: controlType and ProbeID are required.")
			}
		} else {
    	warning("controlData slot not combined: different numbers of rows found.") 
    }
	}

	# history tracking
	history.finished <- as.character(Sys.time())
	# history.command <- capture.output(print(match.call(combine)))  
	history.command <- paste(deparse(match.call(combine)), collapse='') 

	x@history<- rbind(x@history, y@history)
	if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
	lumiVersion <- packageDescription('lumi')$Version
    x@history<- rbind(x@history, 
	       data.frame(submitted=history.submitted, finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	return(x)
})


##some special handling of main is needed
setMethod("boxplot",signature(x="ExpressionSet"),
	function(x, range=0, main, logMode=TRUE, subset=NULL, xlab='', ylab='Amplitude', ...) 
{
  	tmp <- description(x)
  	if (missing(main) && (is(tmp, "MIAME")))
     	main <- tmp@title
	
	if (class(x) == 'MethyLumiM') logMode <- FALSE
	expr <- exprs(x)
  	if (logMode && max(exprs(x), na.rm=TRUE) > 50) {
		## force the expression value as positive in the logMode
		# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
		# remove the negative values
		if (min(expr, na.rm=TRUE) < 0) {
			rMin <- rowMin(expr)
			expr <- expr[rMin > 0, , drop=FALSE]
		}
		expr <- log2(expr)
	} 
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(expr)]
		}
	} else {
		index <- 1:nrow(expr)
	}

	dataMatrix <- expr[index,]
	labels <- colnames(dataMatrix)
	if (is.null(labels)) labels <- as.character(1:ncol(dataMatrix))
	## set the margin of the plot
	mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
	old.mar <- par('mar')
	old.xaxt <- par('xaxt')
	par(xaxt='n')
	par(mar=mar)
	boxplot(dataMatrix ~ col(dataMatrix), main=main, range=range, xlab=xlab, ylab=ylab, ...)
	par(xaxt='s')
	axis(1, at=1:ncol(dataMatrix), labels=labels, tick=TRUE, las=2)
	par(mar=old.mar)
	par(xaxt=old.xaxt)
})


setMethod('hist', signature(x='ExpressionSet'), 
	function(x, ...) 
{
	density(x, ...)
})


setMethod('density', signature(x='ExpressionSet'), 
	function(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l", col=1:dim(x)[2], lty=1:dim(x)[2], 
	lwd=1, xlim=NULL, index.highlight=NULL, color.highlight=2, symmetry=NULL, addLegend=TRUE, legendPos="topright", subset=NULL, main='',...) 
{
	if (is(x, 'ExpressionSet')) {
	    expr <- exprs(x)
	} else {
		stop('Un-supported class of x.')
	}

	if (class(x) == 'MethyLumiM') logMode <- FALSE
		
    if (logMode && (max(expr, na.rm=TRUE) > 50)) {
		## force the expression value as positive in the logMode
		# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
		# remove the negative values
		if (min(expr, na.rm=TRUE) < 0) {
			rMin <- rowMin(expr)
			expr <- expr[rMin > 0, , drop=FALSE]
		}
		expr <- log2(expr)
		if (is.null(xlab)) 
			xlab <- "log2 intensity"
    } else if (is.null(xlab)) 
        xlab <- "intensity"

	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(expr)]
		}
		expr <- expr[index,,drop=FALSE]
	} 

	if (!is.null(symmetry)) {
		x.range <- range(expr)
		if (symmetry > x.range[1] && symmetry < x.range[2]) {
			warning('symmetry point should not be within the range of x.')
			symmetry <- NULL
		} else {
			expr <- rbind(expr, 2*symmetry - expr)
		}
	}
	x.density <- apply(expr, 2, density, ...)
    all.x <- do.call(cbind, lapply(x.density, function(x) x$x))
    all.y <- do.call(cbind, lapply(x.density, function(x) x$y))

	if (!is.null(symmetry)) {
		nr <- nrow(all.x)
		if (all.x[1,1] >= x.range[1] && all.x[1,1] <= x.range[2]) {
			all.x <- all.x[1:round(nr/2),]
			all.y <- all.y[1:round(nr/2),]
		} else {
			all.x <- all.x[round(nr/2):nr,]
			all.y <- all.y[round(nr/2):nr,]
		}
	}
	# matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, lty=lty, lwd=lwd, ...)
	matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, col=col, lty=lty, lwd=lwd, xlim=xlim, main=main)
	if (!is.null(index.highlight)) {
		if (index.highlight > ncol(all.x) || index.highlight < 1) {
			warning('Highlight index out of range.')
			index.highlight <- 1
		}
		lines(all.x[,index.highlight], all.y[,index.highlight], col=color.highlight, lwd=2, lty=1)
	}
	## add legend
	if (addLegend) {
		labels <- colnames(expr)
		if (is.null(labels)) labels <- as.character(1:ncol(expr))

		col <- col
		lwd <- lwd
		lty <- lty
		if (!is.null(index.highlight)) {
			col[index.highlight] <- color.highlight
			lwd[index.highlight] <- 2
			lty[index.highlight] <- 1
		}
		if (length(legendPos) > 1) {
			x.pos <- legendPos[1]
			y.pos <- legendPos[2]
		} else {
			x.pos <- legendPos[1]
			y.pos <- NULL
		}
		legend(x.pos, y.pos, legend=labels, col=col, lwd=lwd, lty=lty, box.col = par("bg"))
	}
	
    invisible(list(all.x = all.x, all.y = all.y))
})


setMethod("pairs", signature(x="ExpressionSet"), 
	function(x, ..., smoothScatter=FALSE, logMode=TRUE, subset=5000, fold=2, dotColor=1, highlight=NULL, highlightColor=2, main=NULL, checkTransform=TRUE) 
{
	upperPanel <- function(x, y) {
		if (smoothScatter) {
			par(new=TRUE)
			smoothScatter(x[subset], y[subset], main='')
		} else {
			if (length(dotColor) == length(x)) {
				points(x[subset], y[subset], col=dotColor[subset], pch='.', cex=3)
			} else {
				points(x[subset], y[subset], col=dotColor, pch='.', cex=3)
			}
			if (!is.null(highlight)) {
				highlight <- intersect(subset, highlight)
				points(x[highlight], y[highlight], col=highlightColor, pch='.', cex=3)
			}
		}
		abline(0, 1, col="red", lty=1)
		if (!is.null(fold) && !is.na(fold)) {
  		if (logMode) {
  			if (checkTransform) {
  				abline(log2(fold), 1, col="blue", lty=2)
  				abline(log2(1/fold), 1, col="blue", lty=2)
  			} else {
  				abline(fold, 1, col="blue", lty=2)
  				abline(-fold, 1, col="blue", lty=2)
  			}
  		} else {
  			abline(fold, 1, col="blue", lty=2)
  			abline(-fold, 1, col="blue", lty=2)
  		}
		}
	}

	lowerPanel <- function(x, y, cex=1.44) {
		txt <- paste("Cor =", as.character(round(cor(x,y),3)),"\n")
		ex <- par("fin")[1]*0.9
	  if (!is.null(fold) && !is.na(fold)) {
  		if (logMode) {
  			if (checkTransform) {
  				up <- length(which((x-y) > log2(fold)))
  				down <- length(which((y-x) > log2(fold)))
  			} else {
  				up <- length(which((x-y) > fold))
  				down <- length(which((y-x) > fold))
  			}
  		} else {
  			up <- length(which((x/y) > fold))
  			down <- length(which((y/x) > fold))
  		}
  		txt <- paste(txt, up, " (> ", fold, ", up)\n", sep="")
  		txt <- paste(txt, down, " (> ", fold, ", down)\n", sep="")
	  }
		text(mean(range(x)), mean(range(y)), labels=txt, cex=ex)
	}

	## put histograms on the diagonal
	diagPanel <- function(x, ...) {
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col="cyan")
	}

	if (smoothScatter) subset <- NULL
	if (is.character(highlight)) {
		# convert it as numeric index
		index <- 1:nrow(x)
		names(index) <- featureNames(x)
		highlight <- index[highlight]
	} 

	expr <- exprs(x)
	if (class(x) == 'MethyLumiM') logMode <- FALSE
	if (checkTransform) {
		if(logMode) {
			if (max(expr, na.rm=TRUE) > 50) {
				## force the expression value as positive in the logMode
				# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
				# remove the negative values
				if (min(expr, na.rm=TRUE) < 0) {
					rMin <- rowMin(expr)
					expr <- expr[rMin > 0, , drop=FALSE]
				}
				expr <- log2(expr)
			}
		} else if (class(x) != 'MethyLumiM') {
			if (max(expr, na.rm=TRUE) < 50) {
				expr <- 2^expr
			}
		}
	}


	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			subset <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			subset <- subset[subset > 0 & subset <= nrow(expr)]
		}
		if (is.null(main)) 
			main <- paste('Pair plot based on', length(subset), 'random sampling.')
	} else {
		subset <- 1:nrow(expr)
	}
	
	
	pairs(expr,upper.panel=upperPanel, diag.panel=diagPanel, 
			lower.panel=lowerPanel, main=main, ...)
})
  	

setMethod("MAplot", signature(object="ExpressionSet"), 
	function(object, ..., smoothScatter=FALSE, logMode=TRUE, subset=5000, main=NULL) 
{
	if (smoothScatter) subset <- NULL
	expr <- exprs(object)
	if (class(object) == 'MethyLumiM') logMode <- FALSE
	if(logMode) {
		if (max(expr, na.rm=TRUE) > 50) {
			## force the expression value as positive in the logMode
			# if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
			# remove the negative values
			if (min(expr, na.rm=TRUE) < 0) {
				rMin <- rowMin(expr)
				expr <- expr[rMin > 0, ,drop=FALSE]
			}
			expr <- log2(expr)
		} 
	} else {
		if (max(expr, na.rm=TRUE) < 50) {
			expr <- 2^expr
		}
	}
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric.')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
			index <- index[index > 0 & index <= nrow(expr)]
		}
	} else {
		index <- 1:nrow(expr)
	}
	if (length(subset) < nrow(expr) && is.null(main)) 
		main <- paste('Pair plot based on', length(subset), 'random sampling.')

	if (smoothScatter) {
		mva.pairs(expr[index, ], log.it=FALSE, plot.method='smoothScatter', main=main, ...)
	} else {
		mva.pairs(expr[index, ], log.it=FALSE, plot.method='normal', main=main, ...)
	}
})


setMethod('plot',
	signature('ExpressionSet', 'missing'),
	function(x, what=c('density', 'boxplot', 'pair', 'MAplot', 'sampleRelation', 'outlier', 'cv'), main, ...)
{
	object <- x
	if (!is(object, 'ExpressionSet')) stop('The object should be class "ExpressionSet".')
	what <- match.arg(what)

	if (what == 'density') {
		if (missing(main)) main <- 'Density plot of intensity'
		density(object, xlab="intensity", ylab="density", main=main, ...)
	} else if (what == 'boxplot') {
		if (missing(main)) main <- 'Boxplot of microarray intensity'
		boxplot(object, xlab='microarrays', ylab='intensity', main=main, ...)
	} else if (what == 'cv') {
		if (missing(main)) main <- 'Density plot of coefficient of variance'
		estimateLumiCV(object, ifPlot=TRUE, main=main, ...)
	} else if (what == 'sampleRelation') {
		plotSampleRelation(object, ...)
	} else if (what == 'pair') {
		if (missing(main)) main <- 'Pairwise plot with sample correlation'
		pairs(object, main=main, ...)
	} else if (what == 'MAplot') {
		if (missing(main)) main <- 'Pairwise MA plots between samples'
		MAplot(object, main=main, ...)
	} else if (what == 'outlier') {
		detectOutlier(object, ifPlot=TRUE, ...)
	} else {
		print('Unsupported .')
	}
	return(invisible(TRUE))	
})


## Example
# library(lumi)
# data(example.lumi)
# test <- asBigMatrix(example.lumi, nCol=6)
setMethod('asBigMatrix',
	signature(object='ExpressionSet'),
	function(object, rowInd=NULL, colInd=NULL, nCol=NULL, dimNames=NULL, saveDir='.', savePrefix=NULL, ...)
{
	## check whether the user just wants a physical copy of the data to a new location
	if (grepl('windows', R.Version()$platform, ignore.case=TRUE)) stop("Package bigmemoryExtras does not support Windows systems.")
	
	if (saveDir == '.') saveDir <- getwd()
	if (is.null(savePrefix)) {
		## use the variable name as the bigmatrix directory prefix
		savePrefix <- match.call(asBigMatrix)[['object']]  
	} 
	if (!file.exists(saveDir)) stop('saveDir does not exist!')
	saveDir <- file.path(saveDir, paste(savePrefix, 'bigmat', sep='_'))
	
	if (is.null(rowInd) && is.null(colInd) && is.null(nCol) && is.null(dimNames) && is(exprs(object), 'BigMatrix')) {
	  fieldnames <- ls(assayData(object)[[assayDataElementNames(object)[1]]])
	  if ('datapath' %in% fieldnames) {
  		oldDir <- dirname(assayData(object)[[assayDataElementNames(object)[1]]]$datapath)
	  } else if ('backingfile' %in% fieldnames) {
  		oldDir <- dirname(assayData(object)[[assayDataElementNames(object)[1]]]$backingfile)
	  }

		if (oldDir != saveDir) {
			if (!file.exists(saveDir)) 	dir.create(saveDir,showWarnings=FALSE)
			sapply(dir(oldDir, full.names=T), file.copy, to=saveDir, overwrite=TRUE, recursive=TRUE)  
		}
		object <- bigmemoryExtras::updateAssayDataElementPaths(object, saveDir)
		return(object)		
	}
	
	if (!is.null(dimNames)) {
		if (!is.list(dimNames) || length(dimNames) != 2) stop("dimNames should be a list with length 2!")
	}
	nRow <- nrow(object)
	if (is.null(nCol)) {
		if (is.null(dimNames)) {
			nCol <- ncol(object)
		} else {
			nCol <- length(dimNames[[2]])
			if (length(nCol) < ncol(object)) stop('The length of dimNames[[2]] should not be less than the columns of the input object!')
		}
	} 
	if (is.null(dimNames)) {
		dimNames <- list(featureNames(object), sampleNames(object))
		## append colnames if nCol is longer than dimNames[[2]]
		if (length(dimNames[[2]]) < nCol) {
			appLen <- nCol - length(dimNames[[2]])
			dimNames[[2]] <- c(dimNames[[2]], paste(rep('unknown', appLen), seq(appLen), sep='.'))
		}
	} else if (length(dimNames[[1]]) < nRow) {
		stop('The length of dimNames[[1]] should be the same as the rows of the input object!')
	}
	subsetMode <- FALSE
	extensionMode <- FALSE 
	if (is.null(rowInd)) {
		rowInd <- 1:nRow
	} else {
		nRow <- length(rowInd)
		dimNames[[1]] <- dimNames[[1]][rowInd]
		subsetMode <- TRUE
	}
	if (is.null(colInd)) {
		colInd <- 1:nCol
	} else {
		nCol <- length(colInd)
		dimNames[[2]] <- dimNames[[2]][colInd]
		subsetMode <- TRUE
	}
	if (nCol != ncol(object)) extensionMode <- TRUE
	
	for (ad.name in assayDataElementNames(object)) {
  	matrix.i <- assayDataElement(object, ad.name)
		if (is.null(matrix.i)) next
		backingfile <- file.path(saveDir, ad.name)
		
		if (!is(assayDataElement(object, ad.name), "BigMatrix") && !extensionMode) {
			x.mat <- bigmemoryExtras::BigMatrix(matrix.i[dimNames[[1]], dimNames[[2]]], backingfile=backingfile, nrow=nRow, ncol=nCol, dimnames=dimNames, ...)
		} else {
			x.mat <- bigmemoryExtras::BigMatrix(backingfile=backingfile, nrow=nRow, ncol=nCol, dimnames=dimNames, ...)
			for (i in 1:ncol(matrix.i)) {
				col.i <- colnames(matrix.i)[i]
				x.mat[dimNames[[1]], col.i] <- matrix.i[dimNames[[1]], col.i]
			}
		}
		assayDataElement(object, ad.name) <- x.mat
  }
	object <- bigmemoryExtras::updateAssayDataElementPaths(object, saveDir)

	if (extensionMode || subsetMode) {
		appLen <- nCol - nrow(pData(object))
		if (length(appLen) > 0) {
			pdata <- rbind(as.matrix(pData(object)), matrix(NA, nrow=appLen, ncol=ncol(pData(object))))
			pdata <- as.data.frame(pdata)
			rownames(pdata) <- dimNames[[2]]
		} else {
			pdata <- pData(object)
			if (length(colInd) < nrow(pdata)) pdata <- pdata[colInd,]
		}
		className <- class(object)
		object.new <- new(className, assayData=assayData(object), phenoData=AnnotatedDataFrame(pdata), 
				featureData=AnnotatedDataFrame(fData(object)[rowInd,,drop=FALSE]))
		object <- object.new
	}
	return(object)
})

