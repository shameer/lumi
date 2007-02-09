`lumiN` <-
function(x.lumi, method=c('RSN', 'loess', 'quantile', 'VSN'), targetArray=NULL,
 		ifPlot=FALSE, ...) {

	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	method <- match.arg(method)
    history.submitted <- as.character(Sys.time())
    new.lumi <- x.lumi 
    if (method == 'VSN') {
		if(!require(vsn)) stop('Package "vsn" should be installed for "VSN" method!')
	}

	exprs(new.lumi) <- switch(method,
		RSN = exprs(lumiN.rsn(x.lumi, targetArray=targetArray, ifPlot=ifPlot, ...)),
		loess = normalize.loess(exprs(x.lumi), ...),
		quantile = normalize.quantiles(x=exprs(x.lumi)),
		VSN = exprs(vsn(intensities=exprs(x.lumi), ...)) )
    
	colnames(exprs(new.lumi)) <- colnames(exprs(x.lumi))
	rownames(exprs(new.lumi)) <- rownames(exprs(x.lumi))
    # history tracking
    history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiN)))
	new.lumi@history<- rbind(new.lumi@history,
	       data.frame(submitted=history.submitted, finished=history.finished, command=history.command))
    
    return(new.lumi)
}

