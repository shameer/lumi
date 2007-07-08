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


setMethod('initialize', 'LumiBatch', function(.Object, 
	exprs = new('matrix'),
	se.exprs = new('matrix'),		# standard deviation of the bead measurements of each gene
	...) 
{
	callNextMethod(.Object, exprs=exprs, se.exprs=se.exprs, ...)
})


setValidity("LumiBatch", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "ExpressionSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("se.exprs")))
    if (is.null(msg)) TRUE else msg
})


##=================================================
## methods
if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
if (is.null(getGeneric("beadNum"))) setGeneric("beadNum", function(object) standardGeneric("beadNum"))
if (is.null(getGeneric("beadNum<-"))) setGeneric("beadNum<-", function(object, value) standardGeneric("beadNum<-"))

if (is.null(getGeneric("detection"))) setGeneric("detection", function(object) standardGeneric("detection"))
if (is.null(getGeneric("detection<-"))) setGeneric("detection<-", function(object, value) standardGeneric("detection<-"))

if (is.null(getGeneric("summary"))) setGeneric("summary", function(object, ...) standardGeneric("summary"))
if (is.null(getGeneric("show"))) setGeneric("show", function(object) standardGeneric("show"))
if (is.null(getGeneric("combine"))) setGeneric("combine", function(x, y, ...) standardGeneric("combine"))
if (is.null(getGeneric("density"))) setGeneric("density", function(x, ...) standardGeneric("density"))

setMethod("se.exprs", signature(object="LumiBatch"),
          function(object) assayDataElement(object,"se.exprs"))

setReplaceMethod("se.exprs", signature(object="LumiBatch",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "se.exprs", value))

setMethod("beadNum", signature(object="LumiBatch"),
          function(object) assayDataElement(object,"beadNum"))

setReplaceMethod("beadNum", signature(object="LumiBatch",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "beadNum", value))

setMethod("detection", signature(object="LumiBatch"),
          function(object) assayDataElement(object,"detection"))

setReplaceMethod("detection", signature(object="LumiBatch",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "detection", value))

setMethod("getHistory",signature(object="LumiBatch"), function(object) object@history)


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
	cat('Summary of BeadStudio output:\n\t')
	cat(notes(object)[[1]], sep='\n\t')
	cat('Major Operation History:\n')
	print(getHistory(object)) 
	cat('\nObject Information:\n')
	callNextMethod()
})


setMethod("combine", signature=c(x="LumiBatch", y="LumiBatch"), function(x, y) 
{
	if (class(x) != class(y))
		stop(paste("objects must be the same class, but are ",
                 class(x), ", ", class(y), sep=""))
	
	if (any(sort(featureNames(x)) != sort(featureNames(y)))) stop('Two data sets have different row names!')
	## determine whether there are duplicated sample names
	sampleName.x <- sampleNames(x)
	sampleName.y <- sampleNames(y)
	if (any(sampleName.x %in% sampleName.y)) {
		warning('Two data sets have some duplicated sample names!\n "_1" and "_2" were attached to the sample names!')
		sampleName.x <- paste(sampleNames(x), '_1', sep='')
		sampleName.y <- paste(sampleNames(y), '_2', sep='')
	}

	history.submitted <- as.character(Sys.time())

	assayData(x) <- combine(assayData(x), assayData(y))
	experimentData(x) <- combine(experimentData(x),experimentData(y))

	## combine pheno data
	phenoData(x) <- combine(phenoData(x),phenoData(y))
	# if (!is.null(phenoData(x)) | !is.null(phenoData(y))) {
	# 	phenoData.x <- phenoData(x)
	# 	phenoData.y <- phenoData(y)
    # 
	# 	pData(phenoData.x) <- rbind(pData(phenoData.x), pData(phenoData.y))
	# 	metaInfo <- rbind(varMetadata(phenoData.x), varMetadata(phenoData.y))
	# 	varMetadata(phenoData.x) <- metaInfo[!duplicated(c(rownames(varMetadata(phenoData.x)),
	# 	 		rownames(varMetadata(phenoData.y)))), ,drop=FALSE]
	# 	phenoData(x) <- phenoData.x
	# }	
	
	## featureData(x) <- combine(featureData(x),featureData(y)) # very slow
	## combine feature data
	if (!is.null(featureData(x)) | !is.null(featureData(y))) {
		feature.x <- featureData(x)
		feature.y <- featureData(y)
    
		pData.x <- pData(feature.x)
		pData.y <- pData(feature.y)
		if (names(pData.x)[1] != names(pData.y)[1])	stop('The featureData of two objects are not incompatible!')
		repInfo <- merge(pData.x, pData.y, by=names(pData.x)[1], all=TRUE, suffixes = c(".x",".y"), sort=FALSE)
		pData(feature.x) <- repInfo
		
		metaInfo <- rbind(varMetadata(feature.x), varMetadata(feature.y)[-1,])
		rownames(metaInfo) <- names(repInfo)
		varMetadata(feature.x) <- metaInfo
		featureData(x) <- feature.x
	}

	## combining the QC information
	if (length(x@QC) > 0 | length(y@QC) > 0) {
		if (!is.null(x@QC$BeadStudioSummary) & !is.null(y@QC$BeadStudioSummary)) {
			if (ncol(x@QC$BeadStudioSummary) == ncol(y@QC$BeadStudioSummary))
			BeadStudioSummary <- rbind(x@QC$BeadStudioSummary, y@QC$BeadStudioSummary)
		} else {
			BeadStudioSummary <- x@QC$BeadStudioSummary
		}
		if (!is.null(x@QC$sampleSummary) & !is.null(y@QC$sampleSummary)) {
			if (nrow(x@QC$sampleSummary) == nrow(y@QC$sampleSummary))
			sampleSummary <- cbind(x@QC$sampleSummary, y@QC$sampleSummary)
		} else {
			sampleSummary <- x@QC$sampleSummary
		}
		
		x@QC$BeadStudioSummary <- BeadStudioSummary
		x@QC$sampleSummary <- sampleSummary
		history.x <- x@QC$history
		if (is.null(history.x)) history.x <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
		if (is.null(history.x$lumiVersion)) history.x$lumiVersion <- rep(NA, nrow(history.x))
		history.y <- y@QC$history
		if (is.null(history.y)) history.y <- data.frame(submitted=NA, finished=NA, command=NA, lumiVersion=NA)
		if (is.null(history.y$lumiVersion)) history.y$lumiVersion <- rep(NA, nrow(history.y))
		x@QC$history <- rbind(history.x, history.y)
	}
	
	## VST transformation parameters
	if (!is.null(attr(x, 'vstParameter')) & !is.null(attr(x, 'vstParameter'))) {
		attr(x, 'vstParameter') <- rbind(attr(x, 'vstParameter'), attr(y, 'vstParameter'))
		attr(x, 'transformFun') <- c(attr(x, 'transformFun'), attr(y, 'transformFun'))
	}
	
	sampleNames(x) <- c(sampleName.x, sampleName.y)

	# history tracking
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(combine)))  
	x@history<- rbind(x@history, y@history)
	if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
	lumiVersion <- packageDescription('lumi')$Version
    x@history<- rbind(x@history, 
	       data.frame(submitted=history.submitted, finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	return(x)
})


##some special handling of main is needed
setMethod("boxplot",signature(x="ExpressionSet"),
	function(x, range=0, main, logMode=TRUE, subset=5000, ...) 
{
  	tmp <- description(x)
  	if (missing(main) && (is(tmp, "MIAME")))
     	main <- tmp@title

	expr <- exprs(x)
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			index <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			index <- subset
		}
	} else {
		index <- 1:nrow(expr)
	}
  	if (logMode & max(exprs(x), na.rm=TRUE) > 50) {
		## force the expression value as positive in the logMode
		if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
		expr <- log2(expr)
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
	boxplot(dataMatrix ~ col(dataMatrix), main=main, range=range, xlab='', ylab='amplitude', ...)
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
	function(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l", index.highlight=NULL, 
	color.highlight=2, symmetry=NULL, addLegend=TRUE, subset=5000, ...) 
{
	if (is(x, 'ExpressionSet')) {
	    expr <- exprs(x)
	} else if (is.numeric(x)) {
		expr <- as.matrix(x)
	} else {
		stop('Un-supported class of x!')
	}
		
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			ind <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			ind <- subset
		}
	} else {
		ind <- 1:nrow(expr)
	}
	expr <- expr[ind,,drop=FALSE]

    if (logMode & (max(expr, na.rm=TRUE) > 50)) {
		## force the expression value as positive in the logMode
		if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1

		expr <- log2(expr)
		if (is.null(xlab)) 
			xlab <- "log2 intensity"
    } else if (is.null(xlab)) 
        xlab <- "intensity"

	if (!is.null(symmetry)) {
		x.range <- range(expr)
		if (symmetry > x.range[1] & symmetry < x.range[2]) {
			warning('symmetry point should not be within the range of x!')
			symmetry <- NULL
		} else {
			expr <- rbind(expr, 2*symmetry - expr)
		}
	}
	x.density <- apply(expr, 2, density)
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
		labels <- colnames(expr)
		if (is.null(labels)) labels <- as.character(1:ncol(expr))

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
})


setMethod("pairs", signature(x="ExpressionSet"), 
	function(x,..., logMode=TRUE, subset=5000) 
{
	upperPanel <- function(x, y, fold=2) {
		points(x[subset], y[subset])
		abline(0, 1, col="red", lty=1)
		if (logMode) {
			abline(log2(fold), 1, col="green", lty=2)
			abline(log2(1/fold), 1, col="green", lty=2)
		} else {
			abline(fold, 1, col="green", lty=2)
			abline(-fold, 1, col="green", lty=2)
		}
	}

	lowerPanel <- function(x, y, cex=1.44, fold=2) {
		if (logMode) {
			up <- length(which((x-y) > log2(fold)))
			down <- length(which((y-x) > log2(fold)))
		} else {
			up <- length(which((x/y) > fold))
			down <- length(which((y/x) > fold))
		}
		ex <- par("fin")[1]*0.9
		txt <- paste("Cor =", as.character(round(cor(x,y),2)),"\n")
		txt <- paste(txt, up, " (> ", fold, ", up)\n", sep="")
		txt <- paste(txt, down, " (> ", fold, ", down)\n", sep="")
		text(mean(range(x)), mean(range(x)), labels=txt, cex=ex)
	}

	## put histograms on the diagonal
	diagPanel <- function(x, ...) {
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
	}

	expr <- exprs(x)
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			subset <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} 
	} else {
		subset <- 1:nrow(expr)
	}
	
	if(logMode) {
		if (max(expr, na.rm=TRUE) > 50) {
			## force the expression value as positive in the logMode
			if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
			expr <- log2(expr)
		}
	} else {
		if (max(expr, na.rm=TRUE) < 50) {
			expr <- 2^expr
		}
	}
	pairs(expr,upper.panel=upperPanel, diag.panel=diagPanel, 
			lower.panel=lowerPanel, ...)
})


if(is.null(getGeneric("MAplot")))
	setGeneric("MAplot", function(object, ...)
		standardGeneric("MAplot"))
  	

setMethod("MAplot", signature(object="ExpressionSet"), 
	function(object, ..., logMode=TRUE, subset=5000) 
{
	expr <- exprs(object)
	if (!is.null(subset)) {
		if (!is.numeric(subset)) stop('subset should be numeric!')
		if (length(subset) == 1) {
			ind <- sample(1:nrow(expr), min(subset, nrow(expr)))
		} else {
			ind <- subset
		}
	} else {
		ind <- 1:nrow(expr)
	}
	if(logMode) {
		if (max(expr, na.rm=TRUE) > 50) {
			## force the expression value as positive in the logMode
			if (min(expr, na.rm=TRUE) < 0) expr <- expr - min(expr, na.rm=TRUE) + 1
			expr <- log2(expr)
		} 
		mva.pairs(expr[ind, ], log.it=FALSE, ...)
	} else {
		if (max(expr, na.rm=TRUE) < 50) {
			expr <- 2^expr
		}
		mva.pairs(expr[ind, ], log.it=FALSE, ...)
	}
})


setMethod("[", "LumiBatch", function(x, i, j, ..., drop = FALSE) 
{
	if (missing(drop)) drop <- FALSE
   	history.submitted <- as.character(Sys.time())
		
	## do default processing of 'ExpressionSet'
	x <- callNextMethod()

	ddim <- dim(x)
	if (!missing(i) & !missing(j)) {
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
			if (!is.null(QC$sampleSummary)) QC$sampleSummary <- QC$sampleSummary[,j,drop=FALSE]
			if (!is.null(QC$BeadStudioSummary)) QC$BeadStudioSummary <- QC$BeadStudioSummary[j,,drop=FALSE]
			x@QC <- QC
		}
		if (!is.null(attr(x, 'vstParameter'))) {
			attr(x, 'vstParameter') <- attr(x, 'vstParameter')[j,,drop=FALSE]
			attr(x, 'transformFun') <- attr(x, 'transformFun')[j]
		}
	}

    # history tracking
    history.finished <- as.character(Sys.time())
	if (is.null(x@history$lumiVersion)) x@history$lumiVersion <- rep(NA, nrow(x@history))
	lumiVersion <- packageDescription('lumi')$Version
	x@history<- rbind(x@history, data.frame(submitted=history.submitted, finished=history.finished, 
			command=history.command, lumiVersion=lumiVersion))

	return(x)
})


setMethod('plot',
	signature('LumiBatch', 'missing'),
	function(x, what=c('density', 'boxplot', 'pair', 'MAplot', 'sampleRelation', 'outlier', 'cv'), main, ...)
{

	object <- x
	if (!is(object, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
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
		print('Unsupported !')
	}
	return(invisible(TRUE))	
})
