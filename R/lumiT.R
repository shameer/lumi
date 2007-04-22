`lumiT` <-
function(x.lumi, method=c('vst', 'log2', 'cubicRoot'), ifPlot=FALSE, simpleOutput = TRUE, ...) {
	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')

	method <- match.arg(method)
	## check the negative values
	if (method %in% c('vst', 'log2')) {
		if (min(exprs(x.lumi)) < 0) {
			x.lumi <- lumiB(x.lumi, method='forcePositive')
		}
	}
	history.submitted <- as.character(Sys.time())

    new.lumi <- x.lumi 
	exprs <- exprs(x.lumi) 
	if (method == 'log2') {
		offset <- ifelse(min(exprs) < 0, -min(exprs), 0)
		exprs(new.lumi) <- log2(exprs + offset + 1)
	} else if (method == 'cubicRoot') {
		exprs(new.lumi) <- sign(exprs) * (abs(exprs))^1/3
	} else {
	    se.exprs <- se.exprs(x.lumi)
	    nArray <- ncol(exprs)
		transExpr <- NULL
		transPara <- NULL
		transFun <- NULL
	    temp <- lapply(1:nArray, function(i) {
				cat(as.character(Sys.time()), ", processing array ", i, "\n")
				if (method == 'vst') {
			        x <- vst(u=exprs[,i], std=se.exprs[,i], method='iterate', ifPlot=ifPlot, ...)
				} else {
					x <- vst(u=exprs[,i], std=se.exprs[,i], method='quadratic', ifPlot=ifPlot, ...)
				}
				transExpr <<- cbind(transExpr, x)
				transPara <<- rbind(transPara, attr(x, 'parameter'))
				transFun <<- c(transFun, attr(x, 'transformFun'))
		        return(TRUE)
			})
		rownames(transPara) <- colnames(exprs(x.lumi))
		names(transFun) <- colnames(exprs(x.lumi))
	    exprs(new.lumi) <- transExpr
	}
	colnames(exprs(new.lumi)) <- colnames(exprs(x.lumi))
	rownames(exprs(new.lumi)) <- rownames(exprs(x.lumi))

	if (simpleOutput) {
		storage.mode <- storageMode(new.lumi)
		if ("lockedEnvironment" == storage.mode) {
			aData <- copyEnv(assayData(new.lumi))
			rm(list=c('se.exprs', 'detection', 'beadNum'), envir=aData)
			lockEnvironment(aData, bindings = TRUE)
			assayData(new.lumi) <- aData
		} else {
			aData <- assayData(new.lumi)
			rm(list=c('se.exprs', 'detection', 'beadNum'), envir=aData)
			assayData(new.lumi) <- aData
		}
	}
	
	# history tracking
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiT)))

	if (is.null(new.lumi@history$lumiVersion)) new.lumi@history$lumiVersion <- rep(NA, nrow(new.lumi@history))
	lumiVersion <- packageDescription('lumi')$Version
	new.lumi@history<- rbind(new.lumi@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))

	if (method == 'vst') {
		attr(new.lumi, 'vstParameter') <- transPara
		attr(new.lumi, 'transformFun') <- transFun
	}

    return(new.lumi)
}

