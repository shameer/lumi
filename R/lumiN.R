`lumiN` <-
function(x.lumi, method=c('quantile', 'rsn', 'ssn', 'loess', 'vsn', 'rankinvariant'), verbose=TRUE, ...) {

	if (is(x.lumi, 'ExpressionSet')) {
	    # x.lumi is a lumi object
	    x.matrix <- exprs(x.lumi)		
	} else if (is.numeric(x.lumi)) {
		x.matrix <- as.matrix(x.lumi)
		if (method == 'rsn') x.lumi <- x.matrix
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}

	method <- match.arg(method)
    if (method == 'vsn') {
		if(!require(vsn)) stop('Package "vsn" should be installed for "vsn" method!')
		if (max(x.matrix, na.rm=TRUE) < 50) {
			warning('The data seems log2 transformed. VSN should be directly applied to the raw data!')
		}
	} 
	
	if (is(x.lumi, 'LumiBatch')) {
		history.submitted <- as.character(Sys.time())
	}
	if (verbose) cat(paste('Perform', method, 'normalization ...\n'))

	norm.matrix <- switch(method,
		rsn = rsn(x.lumi, ...),
		ssn = ssn(x.lumi, ...),
		loess = normalize.loess(x.matrix, ...),
		quantile = normalize.quantiles(x.matrix, ...),
		vsn = exprs(vsn::vsn2(x.matrix, ...)),
		rankinvariant = rankinvariant(x.lumi, ...) )  # add by Arno Velds
	
	## take the average as the parameters as the transformation parameters if there is no targetArray
	if (is.null(attr(norm.matrix, 'targetArray')) && !is.null(attr(x.lumi, 'vstParameter'))) {
		parameter <- attr(x.lumi, 'vstParameter')
		fun <- attr(x.lumi, 'transformFun')
		tt <- table(fun)
		selFun <- names(tt)[which.max(tt)]		# the most frequent function as the transformation function
		selInd <- which(fun == selFun)
		parameter <- colMeans(parameter)
		fun <- selFun
		attr(norm.matrix, 'vstParameter') <- parameter
		attr(norm.matrix, 'transformFun') <- fun
	} 

	if (is.matrix(norm.matrix)) {
		colnames(norm.matrix) <- colnames(x.matrix)
		rownames(norm.matrix) <- rownames(x.matrix)
	} 
	if (!is.null(attr(norm.matrix, 'vstParameter'))) {
		attr(x.lumi, 'vstParameter') <- attr(norm.matrix, 'vstParameter')
		attr(x.lumi, 'transformFun') <- attr(norm.matrix, 'transformFun')
		if (!is.null(attr(norm.matrix, 'targetArray')))
			attr(x.lumi, 'targetArray') <- attr(norm.matrix, 'targetArray')
	}

	if (is(x.lumi, 'ExpressionSet')) {
		new.lumi <- x.lumi
		exprs(new.lumi) <- if (is(norm.matrix, 'ExpressionSet')) exprs(norm.matrix) else norm.matrix
		
	    # history tracking
		if (is(x.lumi, 'LumiBatch')) {
		    history.finished <- as.character(Sys.time())
			history.command <- capture.output(print(match.call(lumiN)))
			if (is.null(new.lumi@history$lumiVersion)) new.lumi@history$lumiVersion <- rep(NA, nrow(new.lumi@history))
			lumiVersion <- packageDescription('lumi')$Version
			new.lumi@history<- rbind(new.lumi@history, data.frame(submitted=history.submitted, 
					finished=history.finished, command=history.command, lumiVersion=lumiVersion))
		}
	    return(new.lumi)
	} else {
		return(norm.matrix)
	}
}

