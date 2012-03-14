## ---------------------------------------------------------------
## define a new class MethyGenoSet
setClass('MethyGenoSet', 
	representation(history='data.frame'), 
	prototype=list(history=data.frame(
		submitted   = I(vector()),
		finished    = I(vector()),
		command     = I(vector()),
		lumiVersion = I(vector())
	)), 
	contains='GenoSet')

#
setMethod('initialize', 'MethyGenoSet', function(.Object, 
	locData = new('RangedData'),
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

		if (ncol(detection) == nSample) {
		  cmd <- paste(cmd, ', detection=detection')
		}

		cmd <- paste(cmd, ')')
		eval(parse(text=cmd))
	} else if (!missing(exprs)) 
		stop("only one of 'assayData' or ('exprs', 'methylated' and 'unmethylated') allowed")

	# callNextMethod(.Object, locData=locData, assayData=assayData, ...)
})


setValidity("MethyGenoSet", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "eSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("methylated")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("unmethylated")))
    if (is.null(msg)) TRUE else msg
})


# setAs("eSet", "MethyGenoSet", function(from) {
# 
#   	history.submitted <- as.character(Sys.time())
# 	
# 	if (exists('methylated', assayData(from)) && exists('methylated', assayData(from))) {
# 		M <- estimateM(from, returnType="matrix")
# 	} else {
# 		stop("Cannot convert as MethyGenoSet object because methylated and unmethylated slots do not exist!\n")
# 	}
# 	
# 	if (.hasSlot(from, "history")) {
# 		history <- from@history
# 		if (is.null(history$lumiVersion)) history$lumiVersion <- rep(NA, nrow(history))
# 	} else {
# 		history <- data.frame(submitted = I(vector()), finished = I(vector()), command = I(vector()), lumiVersion = I(vector()))
# 	}
# 	
# 	aData <- assayData(from)
# 	detection <- NULL
# 	if (exists('pvals', assayData(from))) {
# 		detection <- aData[['pvals']]
# 		storageMode(aData) <- "environment"
# 		aData[['detection']] <- detection
# 		storageMode(aData) <- "lockedEnvironment"
# 	} else if (exists('detection', assayData(from))) {
# 		detection <- aData[['detection']]
# 	}
# 	
# 	methy <- aData[['methylated']]
# 	unmethy <- aData[['unmethylated']]
# 	storageMode(aData) <- "environment"
# 	ts <- ls(envir=aData)
# 	rm(list=ts, envir=aData)
# 	aData[['exprs']] <- M
# 	aData[['methylated']] <- methy
# 	aData[['unmethylated']] <- unmethy
# 	if (!is.null(detection)) aData[['detection']] <- detection
# 
# 	storageMode(aData) <- "lockedEnvironment"	
# 	to <- new("MethyGenoSet", assayData=aData, phenoData=phenoData(from), featureData=featureData(from), annotation=annotation(from), experimentData=experimentData(from), protocolData=protocolData(from))		
# 
# 	history.finished <- as.character(Sys.time())
# 	history.command <- capture.output(print(match.call(setAs)))  
# 	lumiVersion <- packageDescription('lumi')$Version
# 	to@history <- rbind(history, 
#                       data.frame(submitted=history.submitted, 
#                                  finished=history.finished, 
#                                  command=history.command, 
#                                  lumiVersion=lumiVersion))
# 	return(to)
# })

setMethod("exprs", signature(object="MethyGenoSet"), function(object) {
	if ('exprs' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"exprs"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("exprs", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('exprs', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(methylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "exprs", value)
	}
})	


setMethod("methylated", signature(object="MethyGenoSet"), function(object) {
	if ('methylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"methylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("methylated", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(methylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "methylated", value)
	}
})	


setMethod("unmethylated", signature(object="MethyGenoSet"), function(object) {
	if ('unmethylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"unmethylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("unmethylated", signature(object="MethyGenoSet"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(unmethylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "unmethylated", value)
	}
})	
	


setMethod("detection", signature(object="MethyGenoSet"), function(object) {
	if ('detection' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"detection"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("detection", signature(object="MethyGenoSet"), function(object, value) {
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


setMethod("getHistory",signature(object="MethyGenoSet"), function(object) object@history)


setMethod("combine", signature=c(x="MethyGenoSet", y="MethyGenoSet"), function(x, y, ...) {

   	history.submitted <- as.character(Sys.time())

	if (missing(y)) return(x)
	if (length(list(...)) > 0) 
	        return(combine(x, combine(y, ...)))

  	## do default processing of 'GenoSet'
  	x.comb <- callNextMethod()

    # history tracking
    history.finished <- as.character(Sys.time())
	#history.command <- match.call()
    history.command <- capture.output(print(match.call(combine)))  
	x.comb@history<- rbind(x@history, y@history)
	if (is.null(x.comb@history$lumiVersion) && nrow(x@history) > 0) {
		x.comb@history <- data.frame(x.comb@history, lumiVersion=rep(NA, nrow(x.comb@history)))
	} 
	lumiVersion <- packageDescription('lumi')$Version
	x.comb@history<- rbind(x.comb@history, data.frame(submitted=history.submitted,finished=history.finished,command=history.command, lumiVersion=lumiVersion))
	return(x.comb)
})


setMethod("[", "MethyGenoSet", function(x, i, j, ..., drop = FALSE)  {

	if (missing(drop)) drop <- FALSE
	history.submitted <- as.character(Sys.time())
	
	## do default processing of 'GenoSet'
	x <- callNextMethod()
	
	ddim <- dim(x)
	if (!missing(i) & !missing(j)) {
	  	history.command <- paste('Subsetting', ddim[1], 'features and', ddim[2], 'samples.')		
	} else if (!missing(i)) {
	  	history.command <- paste('Subsetting', ddim[1], 'features.')
	} else if (!missing(j)) {
	  	history.command <- paste('Subsetting', ddim[2], 'samples.')
	  	if (!is.null(controlData(x)))
	  		controlData(x) <- controlData(x)[,j]
	} else {
	  	return(x)
	}
	
	# history tracking
	history.finished <- as.character(Sys.time())
	if (is.null(x@history$lumiVersion) && nrow(x@history) > 0) {
		x@history <- data.frame(x@history, lumiVersion=rep(NA, nrow(x@history)))
	}
	lumiVersion <- packageDescription('lumi')$Version
	x@history<- rbind(x@history, data.frame(submitted=history.submitted,finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	
	return(x)
})




