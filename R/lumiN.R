`lumiN` <-
function(x.lumi, method=c('rsn', 'loess', 'quantile', 'vsn'), ...) {

	if (is(x.lumi, 'ExpressionSet')) {
	    # x.lumi is a lumi object
	    x.matrix <- exprs(x.lumi)		
	} else if (is.numeric(x.lumi)) {
		x.matrix <- as.matrix(x.lumi)
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

	norm.matrix <- switch(method,
		rsn = rsn(x.matrix, ...),
		loess = normalize.loess(x.matrix, ...),
		quantile = normalize.quantiles(x=x.matrix, ...),
		vsn = exprs(vsn::vsn2(intensities=x.matrix, ...)) )

	colnames(norm.matrix) <- colnames(x.matrix)
	rownames(norm.matrix) <- rownames(x.matrix)

	if (is(x.lumi, 'ExpressionSet')) {
		new.lumi <- x.lumi
		exprs(new.lumi) <- norm.matrix
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

