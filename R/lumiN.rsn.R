`lumiN.rsn` <-
function(x.lumi, targetArray=NULL, excludeFold=2, span=0.03, ifPlot=FALSE,...) {
	if (is(x.lumi, 'ExpressionSet')) {
	    # x.lumi is a lumi object
	    exprs <- exprs(x.lumi)		
	} else if (is.numeric(x.lumi)) {
		exprs <- as.matrix(x.lumi)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}
	
	## check whether the data was variance stabilized.
	if (max(exprs, na.rm=TRUE) > 100) {
		if (is(x.lumi, 'LumiBatch')) {
			warning('The data has not been variance stabilized!')
			print('Perform VST transform ...')
			x.lumi <- lumiT(x.lumi)
			exprs <- exprs(x.lumi)
		} else {
			warning('The data has not been variance stabilized!')
			print('Perform log2 transform ...')
			exprs <- log2(exprs)
		}
	}

    ## Define interal function 
   	pairwiseN <- function(ind, exprs, targetArray, method=c('RSN', 'loess'), 
				exprs0=NULL, ifPlot=FALSE) {
        cat(as.character(Sys.time()), ", processing array ", ind, "\n")
    
        # normal array ind against targetArray      
        u1 <- exprs[,ind]
        u2 <- exprs[, targetArray]
        u1.original <- u1
        u2.original <- u2
        
		## define window functions
		win <- function(x, sigma=1, type=c('gaussian', 'tricube', 'bisqure')) {
			type <- match.arg(type)
			ww <- switch(type,
				gaussian = exp(-x^2/(2*sigma^2)) / (sqrt(2*pi)*sigma),
				tricube = (1 - (abs(x)/(3 * sigma))^3)^3,
				bisquare = (1 - (x/(3 * sigma))^2)^2 )
			return(ww)			
		}

		if (ind != targetArray) {
			## calculate the weights based on the fold change 
			if (!is.null(exprs0)) {
				fd <- exprs0[,ind] - exprs0[, targetArray]
				if (!is.null(excludeFold)) {
					sigma <- log2(excludeFold)/3
				} else {
					sigma <- sd(fd)
				}
 				wt <- win(abs(fd), sigma)
				wt <- wt/max(wt)
			} else {
				wt <- rep(length(u1))
			}
		
			u1.normalized <- monoSmu(u1, u2, newX=u1.original, span=span, ifPlot=ifPlot, wt=wt, rotate=TRUE,
					xlab=paste("array", ind), ylab=paste("array", targetArray), ...)
        } else {
            u1.normalized <- u1
        }
        return(u1.normalized)
    }

    nArray <- ncol(exprs)
	
    if (ifPlot) par(mfrow=c(2,2))
    
	## do quantile normalization for the purpose of estimating fold change
	## Based on the estimated fold change, we can down-weight the differentiated genes.
    if (!is.null(excludeFold)) {
        exprs0 <- normalize.quantiles(exprs)
    } else {
        exprs0 <- NULL
    }

    if (is.null(targetArray)) {
        # find the sample which is the most similar to the mean profile of all samples,
 		meanProfile <- apply(exprs, 1, mean)
		targetArray <- which.min(abs(colSums(exprs - meanProfile)))
    }

    normalized <- lapply(1:nArray, FUN=pairwiseN, exprs=exprs, method=method, 
					exprs0=exprs0, targetArray=targetArray, ifPlot=ifPlot)
	if (is(x.lumi, 'ExpressionSet')) {
		exprs(x.lumi) <- matrix(unlist(normalized), ncol=nArray, byrow=FALSE)
	} else {
		x.lumi <- matrix(unlist(normalized), ncol=nArray, byrow=FALSE)
	}
    
    return(x.lumi)
}
