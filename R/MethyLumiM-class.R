## ---------------------------------------------------------------
## define a new class lumiMethyM
setClass('MethyLumiM', 
	representation(history='data.frame'), 
	prototype=list(history=data.frame(
		submitted   = I(vector()),
		finished    = I(vector()),
		command     = I(vector()),
		lumiVersion = I(vector())
	)),
	contains='ExpressionSet')

#
setMethod('initialize', 'MethyLumiM', function(.Object, 
	exprs = new('matrix'),
	methylated = new('matrix'),		
	unmethylated = new('matrix'),
	detection = new('matrix'),  # detection pvalues
    ...,
    assayData)
{
	if (missing(assayData)) {
		cmd <- 'assayData <- assayDataNew(exprs=exprs, methylated=methylated, unmethylated=unmethylated'
		nSample <- ncol(exprs)
		if (ncol(detection) == nSample) cmd <- paste(cmd, ', detection=detection')
		cmd <- paste(cmd, ')')
		eval(parse(text=cmd))
	} else if (!missing(exprs)) 
		stop("only one of 'assayData' or ('exprs', 'methylated' and 'unmethylated') allowed")

	callNextMethod(.Object, assayData=assayData, ...)
})


setValidity("MethyLumiM", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "ExpressionSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("methylated")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("unmethylated")))
    if (is.null(msg)) TRUE else msg
})


setAs("eSet", "MethyLumiM", function(from) {

  	history.submitted <- as.character(Sys.time())

	from <- asS4(from)
	
	if (exists('methylated', assayData(from)) && exists('methylated', assayData(from))) {
		M <- estimateM(from)
	} else {
		stop("Cann't convert as MethyLumiM object because methylated and unmethylated slots do not exist!\n")
	}
	aData <- assayData(from)
	storageMode(aData) = "environment"
	aData[['exprs']] = M
	storageMode(aData) = "lockedEnvironment"
	
	if ("history" %in% slotNames(from)) {
		history <- from@history
		if (is.null(history$lumiVersion)) history$lumiVersion <- rep(NA, nrow(history))
	} else {
		history <- data.frame(submitted = I(vector()), finished = I(vector()), command = I(vector()), lumiVersion = I(vector()))
	}

	to <- new("MethyLumiM", assayData=aData, phenoData=phenoData(from), featureData=featureData(from), annotation=annotation(from), experimentData=experimentData(from), protocolData=protocolData(from))
    history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(setAs)))  
	lumiVersion <- packageDescription('lumi')$Version
	to@history<- rbind(history, data.frame(submitted=history.submitted, finished=history.finished, 
			command=history.command, lumiVersion=lumiVersion))
	
	return(to)
})


if (is.null(getGeneric("methylated"))) setGeneric("methylated", function(object) standardGeneric("methylated"))
if (is.null(getGeneric("methylated<-"))) setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))

if (is.null(getGeneric("unmethylated"))) setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
if (is.null(getGeneric("unmethylated<-"))) setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))


setMethod("methylated", signature(object="MethyLumiM"), function(object) {
	if ('methylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"methylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("methylated", signature(object="MethyLumiM"), function(object, value) {
		if (is.null(value)) {
			assay <- assayData(object)
			if (exists('methylated', envir=assay)) {
				oldMode <- storageMode(assay)
				storageMode(assay) <- 'environment'
				rm(detection, envir=assay)
				storageMode(assay) <- oldMode
				assayData(object) <- assay
			}
			return(object)
		} else {
			assayDataElementReplace(object, "methylated", value)
		}
	})	


setMethod("unmethylated", signature(object="MethyLumiM"), function(object) {
	if ('unmethylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"unmethylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("unmethylated", signature(object="MethyLumiM"), function(object, value) {
		if (is.null(value)) {
			assay <- assayData(object)
			if (exists('methylated', envir=assay)) {
				oldMode <- storageMode(assay)
				storageMode(assay) <- 'environment'
				rm(detection, envir=assay)
				storageMode(assay) <- oldMode
				assayData(object) <- assay
			}
			return(object)
		} else {
			assayDataElementReplace(object, "unmethylated", value)
		}
	})	

setMethod("getHistory",signature(object="MethyLumiM"), function(object) object@history)


setMethod("boxplot",signature(x="MethyLumiM"),
	function(x, main, subset=NULL, ...) {
		
	if (!require(hdrcde)) stop("Please install the required hdrcde package./n")  

  	tmp <- description(x)
  	if (missing(main) && (is(tmp, "MIAME")))
     	main <- tmp@title

	expr <- exprs(x)

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

	tmp <- lapply(1:ncol(dataMatrix), function(i) dataMatrix[,i])
	hdr.boxplot(tmp, main=main, xlab='', ylab='M-value', ...)
	par(xaxt='s')
	axis(1, at=1:ncol(dataMatrix), labels=labels, tick=TRUE, las=2)
	par(mar=old.mar)
	par(xaxt=old.xaxt)
})

