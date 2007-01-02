##=================================================
## Define LumiQC
# new('LumiQC', mean=mm, std=std, sampleCor=sampleCor, cv=cv, AP=AP, sampleRelation=sampleRelation, outlier=outlier)
setClass("LumiQC", 
	representation(
		mean = "numeric",     		# mean of the sample
		std = "numeric",   			# standard deviation of the sample
		sampleCor = 'matrix',			# correlation between samples
		detectionRate = "numeric",	# the percentage of detected genes
		sampleRelation = "list",	# the sample relations
		outlier = 'logical',		# identification of outlier samples
		history= "data.frame"		# operation history
	),
	prototype = list(
		mean = numeric(),
		std = numeric(),   		
		sampleCor = matrix(nr=0, nc=0),		
		detectionRate = numeric(),
		sampleRelation = list(),
		outlier = logical(),
		history = data.frame(
		   submitted   = I(vector()),
		   finished    = I(vector()),
		   command     = I(vector())
		)      
	),
	contains="ExpressionSet"
)


setMethod('initialize', 'LumiQC', function(.Object, 
	exprs = new('matrix'),
	cv = new('matrix'),				# coefficients of variance of each measurement
	mean, std, sampleCor, detectionRate, sampleRelation, outlier, history) 
{
	if(!missing(mean)) .Object@mean <- mean
	if(!missing(std)) .Object@std <- std
	if(!missing(sampleCor)) .Object@sampleCor <- sampleCor
	if(!missing(detectionRate)) .Object@detectionRate <- detectionRate
	if(!missing(sampleRelation)) .Object@sampleRelation <- sampleRelation
	if(!missing(outlier)) .Object@outlier <- outlier
	if(!missing(history)) .Object@history <- history
	callNextMethod(.Object, exprs=exprs, cv=cv)
})


##=================================================
## Methods
setMethod("summary",signature(object="LumiQC"), function(object) 
{
	show(object)
})

setMethod("show",signature(object="LumiQC"), function(object) 
{
  	cat(paste("Class:", class(object)[1], '\n'))
	dimen <- dim(exprs(object))
  	cat(paste("Data dimension: ", paste(dimen[1], 'genes', 'x', dimen[2], 'samples', collapse=" "), '\n'))
	detectTh <- attr(object@detectionRate, 'threshold')
	outlier <- object@outlier
	distCenter <- attr(outlier, 'sampleDistance')
	temp <- rbind(object@mean, object@std, object@detectionRate, distCenter[2:nrow(distCenter),1])
	rownames(temp) <- c('mean', 'standard deviation', 'detection rate', 'distance to center')
  	cat(paste("\nSummary of Samples:\n", sep=''))
	print(signif(temp, 4), quote=FALSE)
	Th <- attr(outlier, 'threshold')
 	cat(paste('\nNote: the detection threshold is ', detectTh, ', "center" is the mean of all samples.\n', sep=''))

	if (length(which(outlier)) > 0) {
		cat(paste("\nDetected outliers based on distance-to-center threshold ", signif(Th, 3), ":\n", sep=""))
		sampleName <- names(outlier)
		attributes(outlier) <- NULL
		names(outlier) <- sampleName
		cat(sampleName[outlier])
	} else {
		cat(paste("\nNo outliers detected based on distance-to-center threshold ", signif(Th, 3), ".\n", sep=""))
	}

	cat('\nMajor Operation History:\n')
	print(getHistory(object)) 
})


if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object, ...)
		standardGeneric("getHistory"))
setMethod("getHistory",signature(object="LumiQC"), function(object) object@history)


setMethod("[", "LumiQC", function(x, i, j, ..., drop = FALSE) 
{
	if (missing(drop)) drop <- FALSE
   	history.submitted <- as.character(Sys.time())
	
	if (!missing(j)) {
		x@mean <- x@mean[j]
		x@std <- x@std[j]
		x@detectionRate <- x@detectionRate[j]
		x@sampleCor <- x@sampleCor[j,j]
		## subset the sample relation slot
		sampleRelation <- as.matrix(x@sampleRelation[[1]])
		x@sampleRelation <- list(as.dist(sampleRelation[j,j]))
		## subset the outlier slot
		outlier <- x@outlier
		sampleDist <- attr(outlier, 'sampleDistance')
		sampleDist <- sampleDist[c(1, 1+j), c(1,1+j)]
		outlier <- outlier[j]
		attr(outlier, 'sampleDistance') <- sampleDist
		attr(outlier, 'threshold') <- attr(x@outlier, 'threshold')
		x@outlier <- outlier		
	}
	
	## do default processing of 'ExpressionSet'
	x <- callNextMethod()

	ddim <- dim(x)
	if (!missing(i) & !missing(j)) {
		history.command <- paste('Subsetting', ddim[1], 'features and', ddim[2], 'samples.')		
	} else if (!missing(i)) {
		history.command <- paste('Subsetting', ddim[1], 'features.')
	} else if (!missing(i)) {
		history.command <- paste('Subsetting', ddim[2], 'samples.')
	} else {
		return(x)
	}
	
    # history tracking
    history.finished <- as.character(Sys.time())
    x@history<- rbind(x@history, c(history.submitted, history.finished, history.command))

	return(x)
})


setMethod('plot',
	signature('LumiQC', 'missing'),
	function(x, what=c('density', 'boxplot', 'pair', 'MAplot', 'sampleRelation', 'outlier', 'cv'), parameterList=NULL, main, ...)
{

	object <- x
	if (!is(object, 'LumiQC')) stop('The object should be class "LumiQC"!')
	what <- match.arg(what)
	
	if (what == 'density') {
		if (missing(main)) main <- 'Density plot of intensity'
		hist(object, xlab="intensity", ylab="density", main=main, ...)
	} else if (what == 'boxplot') {
		if (missing(main)) main <- 'Boxplot of microarray intensity'
		boxplot(object, xlab='microarrays', ylab='intensity', main=main, ...)
	} else if (what == 'cv') {
		if (missing(main)) main <- 'Density plot of coefficient of variance'
		cv <- assayDataElement(object,"cv")
		plotDensity(cv, xlab='coefficient of variance', main=main, ...)
	} else if (what == 'sampleRelation') {
		dd <- object@sampleRelation[[1]]
		Th <- attr(dd, 'threshold')
		geneNum <- attr(dd, 'geneNum')
		if (is.null(parameterList$method)) {
			method <- 'cluster'
		} else {
			method <- parameterList$method
		}
		if (method == 'cluster') {
			hc <- hclust(dd, 'ave')
			if (missing(main)) main <- paste('Clusters of the samples based on', geneNum, 'genes with sd/mean >', Th)
 			plot(hc, xlab='Sample', main=main, ...)
 		} else {
 			if (is.null(parameterList$dimension)) {
 				dimension <- c(1,2)
 			} else {
 				dimension <- parameterList$dimension
 			}
 			col <- parameterList$col		
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
 			plot(a1[,dimension[1]],a1[,dimension[2]], type='n', xlab=paste('Principle component', dimension[1]), 
				ylab=paste('Principle component', dimension[2]))
			text(a1[,dimension[1]],a1[,dimension[2]], col=color, labels=sampleNames(object), cex=1)
			title(paste('Sample relations based on', geneNum, 'genes with sd/mean >', Th))			
		}
	} else if (what == 'pair') {
		if (missing(main)) main <- 'Pairwise plot with sample correlation'
		pairs(object, main=main, ...)
	} else if (what == 'MAplot') {
		if (missing(main)) main <- 'Pairwise MA plots between samples'
		MAplot(object, main=main, ...)
	} else if (what == 'outlier') {
		outlier <- object@outlier
		d <- attr(outlier, 'sampleDistance')
		Th <- attr(outlier, 'threshold')
		hc <- hclust(as.dist(d), 'ave')
		if (missing(main)) main <- paste('Outlier detection based on sample distance to "Center"')
		plot(hc, xlab='Sample', ylab='Distance', main=main, ...)
		abline(h=Th, col=2, lty=2)
	} else {
		print('Unsupported !')
	}
	return(invisible(TRUE))	
})


