lumiB <- function(x.lumi, method = c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy'), verbose=TRUE, ...) 
{
	# method <- match.arg(method)
	if (is(x.lumi, 'ExpressionSet')) {
	    # x.lumi is a LumiBatch object
	    x.matrix <- exprs(x.lumi)		
	} else if (is.numeric(x.lumi)) {
		x.matrix <- as.matrix(x.lumi)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}
	
	if (!(is.function(method))) {
		method <- method[1]
		if (!(method %in% c('bgAdjust', 'none', 'forcePositive', 'bgAdjust.affy'))) {
			cat('This method is not supported!\n')
			return(x.lumi)
		} else if (method == 'none') {
			return(x.lumi)
		} else if (method == 'forcePositive') {
			if (min(x.matrix, na.rm=TRUE) > 0)  return(x.lumi)
		}		
		if (verbose) cat(paste('Perform', method, 'background correction ...\n'))
	} else {
		if (verbose) cat('Perform user provided background correction ...\n')		
	}
	
	history.submitted <- as.character(Sys.time())
	if (is.function(method)) {
		x.matrix <- method(x.lumi, ...)
		if (is(x.matrix, 'ExpressionSet')) x.matrix <- exprs(x.matrix)
	} else {
		if (method == 'bgAdjust') {
			x.matrix <- exprs(bgAdjust(x.lumi, ...))
		} else if (method == 'bgAdjust.affy') {
			x.matrix <- apply(x.matrix, 2, bg.adjust, ...) 
		} else if (method == 'forcePositive') {
			offset <- apply(x.matrix, 2, min, na.rm=TRUE)
			offset[offset <= 0] <- offset[offset <= 0] - 1.01 	# to avoid higher fold-change when one value is less than 1
			offset[offset > 0] <- 0
			offset <- rep(1, nrow(x.matrix)) %*% t(offset)
			x.matrix <- x.matrix - offset
		} else {
			cat('The method is not supported!\n')
			return(x.lumi)
		}
	}

	if (is(x.lumi, 'ExpressionSet')) {
		if (!is.function(method)) exprs(x.lumi) <- x.matrix
		if (is(x.lumi, 'LumiBatch')) {
			# history tracking
			history.finished <- as.character(Sys.time())
			history.command <- capture.output(print(match.call(lumiB)))
        	
			if (is.null(x.lumi@history$lumiVersion)) x.lumi@history$lumiVersion <- rep(NA, nrow(x.lumi@history))
			lumiVersion <- packageDescription('lumi')$Version
			x.lumi@history<- rbind(x.lumi@history, data.frame(submitted=history.submitted, 
					finished=history.finished, command=history.command, lumiVersion=lumiVersion))
		}
	} else {
		x.lumi <- x.matrix
	}
	return(x.lumi)
}