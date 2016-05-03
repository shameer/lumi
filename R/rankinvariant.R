`rankinvariant` <-
function(x.lumi, targetArray=NULL, rrc=.05, lowRank = seq(0.5, 0.25, -0.05), highRank=.9, minSize=.02, maxit=200) {
	if (is(x.lumi, 'ExpressionSet')) {
		expr  <- exprs(x.lumi)
	} else if (is.numeric(x.lumi)) {
		expr <- as.matrix(x.lumi)
	} else {
		stop('The object should be a matrix or class "ExpressionSet" inherited!')
	}

	stopifnot(min(lowRank)>0)
	stopifnot(max(lowRank)<1)
	stopifnot(length(highRank) == 1)
	stopifnot(highRank < 1)
	stopifnot(max(lowRank) < highRank)

	externalTarget <- FALSE
	if (!is.null(targetArray)) {
		## check the format of the targetArray
		if (is(targetArray, 'ExpressionSet')) {
			targetArray <- exprs(targetArray)[,1]
		} 
		if (length(targetArray) > 1) {
			if (length(targetArray) != nrow(expr)) stop('targetArray should be an index or a vector has the same length as other samples.')
			expr <- cbind(targetArray, expr)
			targetArray <- 1
			externalTarget <- TRUE
		}
		if (is.numeric(targetArray)) {
			if (targetArray < 1 || targetArray > nrow(expr)) {
				warning('The provided targetArray is invalid and will be set as NULL!')
				targetArray <- NULL
			}
		} else if (is.character(targetArray)) {
			if (!(targetArray %in% colnames(expr))) {
				warning('The provided targetArray is invalid and will be set as NULL!')
				targetArray <- NULL
			}
		} else {
			warning('The provided targetArray is invalid and will be set as NULL!')
			targetArray <- NULL
		}
	}

	## check whether the data was variance stabilized.
	if (max(expr, na.rm=TRUE) > 100) {
		log2Trans <- FALSE
	} else {
		log2Trans <- TRUE
		expr <- 2^(expr)
	}

	#add virtual sample if target is NULL
	if(is.null(targetArray)) {
		expr <- cbind(rowMeans(expr),expr)
		targetArray <- 1
		externalTarget <- TRUE
	}
	nrows <- nrow(expr)
	nArray <- ncol(expr)

	#prepare normalized data obj
	normalized <- expr

	#rank all data and virtual
	er <- apply(expr, 2, rank)

	## determine the rank invariant geneset
	ri <- rep(F, nrow(expr))
	rankscores <- apply(er, 2, function(r, v) abs( (r - v )/length(v) ), er[,targetArray])
	for(i in 1:nArray) {
		if(i == targetArray) next;
		for(lowrank in lowRank) {	
			if(getOption("verbose")) cat("Sample", i - ifelse(externalTarget,1,0), "using lowrank ", lowrank,"\n")
			ri <- rankscores[,i] < rrc & er[,i] > nrows*lowrank & er[,i] < nrows*highRank
			if(getOption("verbose")) cat("Nr rankinvariant probes=", sum(ri),"\n")
			if(sum(ri) >= minSize * nrows) {
				#normalize this sample
				coef <- rlm(expr[ri, targetArray] ~ expr[ri, i], psi=psi.bisquare, maxit=maxit)$coef
				normalized[,i] <- (expr[,i]*coef[2]) + coef[1]
				break;
			}
		}
		if(lowrank == tail(lowRank,1)) stop(paste("No rankinvariant set for sample", i - ifelse(externalTarget,1,0),"use broader lowRank/highRank or use another normalization method"))
	}

	#remove virtual from normalized data if target is external of rowMeans
	if(externalTarget) normalized <- normalized[,-1]

	## transformed as original scale in not log2transformed
	if (log2Trans) {
		if (min(normalized) <= 0) {
			normalized <- lumiB(normalized, method='forcePositive')
		}
		normalized <- log2(normalized)	
	} 
	
	if (is(x.lumi, 'ExpressionSet')) {
		exprs(x.lumi) <- normalized
	} else {
		x.lumi <- normalized
	}
	
	return(x.lumi)
}
