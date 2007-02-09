`lumiT` <-
function(x.lumi, method=c("vst", "vst.quadratic", 'log2', 'cubicRoot'), ifPlot=FALSE, ...) {
	if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')

	method <- match.arg(method)
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
				transPara <<- c(transPara, list(attr(x, 'parameter')))
				transFun <<- c(transFun, attr(x, 'transformFun'))
		        return(TRUE)
			})
	    exprs(new.lumi) <- transExpr
	}
	colnames(exprs(new.lumi)) <- colnames(exprs(x.lumi))
	rownames(exprs(new.lumi)) <- rownames(exprs(x.lumi))

	# history tracking
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiT)))
	new.lumi@history<- rbind(new.lumi@history,
	       data.frame(submitted=history.submitted, finished=history.finished, command=history.command))

	if (method == 'vst') {
		attr(new.lumi, 'parameter') <- transPara
		attr(new.lumi, 'transformFun') <- transFun
	}

    return(new.lumi)
}

