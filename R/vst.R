`vst` <-
function(u, std, nSupport=min(length(u), 500), method=c('iterate', 'quadratic'), ifPlot=FALSE) {
	# u is the mean of probe beads
	# std is the standard deviation of the probe beads
	
	method <- match.arg(method)
	ord <- order(u); u.bak <- u
	u <- u[ord]; std <- std[ord]
	
	## check for the negative values, which is not allowed in the vst function
	if (any(u < 0)) {
		stop('Negative expression value is not allowed!')
	}
	if (any(std < 0)) {
		stop('Negative expression standard deviation is not allowed!')
	}
	
 	## downsampling to speed up
	downSampledU <- 2^seq(from=min(log2(u), na.rm=TRUE), to=max(log2(u), na.rm=TRUE), length=nSupport)
	
	if (method == 'quadratic') {
		dd <- data.frame(y=std, x2=u^2, x1=u)
		##dd = data.frame(y=std^2, x2=u^2, x1=u)
		lm2 <- lm(y ~ x2 + x1, dd)
		smoothStd <- predict(lm2, data.frame(x2=downSampledU^2, x1=downSampledU))
	} else {
		minU <- max(log2(100), log(min(u)))
		maxU <- log2(max(u))
		# uCutoff <- 2^((maxU + minU)/2)
		uCutoffLow <- 2^(minU + (maxU - minU)/3)
		uCutoffHigh <- 2^(minU + (maxU - minU) * 3/4)
		selInd <- (u > uCutoffLow & u < uCutoffHigh)
		selLowInd <- (u < uCutoffLow)
		iterNum <- 0
		c3.i <- 0
		while(iterNum < 3) {
			selInd.i <- selInd & (std^2 > c3.i)
			dd <- data.frame(y=sqrt(std[selInd.i]^2 - c3.i), x1=u[selInd.i])
			if (nrow(dd) > 5000) dd <- dd[sample(1:nrow(dd), 5000),]
			lm.i <- lm(y ~ x1, dd)
			c1.i <- lm.i$coef[2]
			c2.i <- lm.i$coef[1]
			y <- std[selLowInd]
			x <- u[selLowInd]
			cc <- y^2 - (c1.i * x + c2.i)^2
			c3.i <- mean(cc, trim=0.05)
			if (c3.i < 0) {
				c3.i <- 0
				break
			}
			iterNum <- iterNum + 1
		}
		c1 <- c1.i; c2 <- c2.i; c3 <- c3.i
		if (c3 < 0) c3 <- 0
		smoothStd <- ((c1 * downSampledU + c2)^2 + c3)^(1/2)
	}

	if (ifPlot) {
		par(mfrow=c(1,2))
		len <- length(u)
		ind <- sample(1:len, min(5000, len))
		plot(u[ind], std[ind], pch='.', log='xy', xlab="mean", ylab="standard deviation")
		lines(downSampledU, smoothStd, col=2, lwd=1.5)
	}
	
	## calculate the integration (h function is the integral)
	if (method == 'iterate') {
		if (c3 == 0) {
			## Transform function h(x) = g * log(a + b * x)
			g <- 1/c1
			a <- c2
			b <- c1
			tmp <- a + b * u.bak
			if (any(tmp) < 0) {
				transformedU <- log(u.bak)
				g <- 1; a <- 0; b <- 1
			} else {
				transformedU <- g * log(a + b * u.bak)
			}
			if (ifPlot) hy <- g * log(a + b * downSampledU) 
			transFun <- 'log'
		} else {
			## Transform function h(x) = g * asinh(a + b * x)
			g <- 1/c1
			a <- c2/sqrt(c3)
			b <- c1/sqrt(c3)
			transformedU <- g * asinh(a + b * u.bak)
			if (ifPlot) hy <- g * asinh(a + b * downSampledU) 
			transFun <- 'asinh'
		}
		transform.parameter <- c(a, b, g, 0)
		names(transform.parameter) <- c('a', 'b', 'g', 'Intercept')
	} else {
		dd <- diff(downSampledU)	## the interval between samples
		dd <- c(dd[1], dd)
		hy <- cumsum(1/smoothStd * dd)		
		# get a smoothed version and interpolation of original data order
	    transformedU <- monoSpline(x=downSampledU, y=hy, newX=u.bak, nKnots=20, ifPlot=FALSE)
		transFun <- 'quadratic'
		transform.parameter <- NULL
	}	

    if (ifPlot) {
        plot(downSampledU, hy, col=2, xlab="original value", ylab="variance stabilization transformed value")
        ii <- order(u.bak)
        lines(u.bak[ii], transformedU[ii], col=3, type='l')
		par(mfrow=c(1,1))
    }

	## rescale to the similar range with log2
	if (method == 'iterate') {
		cutInd <- which.min(abs(u.bak - uCutoffLow))
		maxInd <- which.max(u.bak)
		y <- c(u.bak[cutInd], u.bak[maxInd])
		x <- c(transformedU[cutInd], transformedU[maxInd])
		m <- lm(log2(y) ~ x)
		transform.parameter <- c(a, b, g * m$coef[2], m$coef[1])
		names(transform.parameter) <- c('a', 'b', 'g', 'Intercept')
		## The transform parameter is in the transFun below
		## transFun <- g * asinh(a + b * x) * m$coef[2] + m$coef[1]
		transformedU <- predict(m, data.frame(x=transformedU))
	} else {
		minRequire <- length(which(u.bak < 1)) / length(u.bak)
		low <- ifelse(minRequire < 0.05, 0.05, minRequire)
		y <- c(quantile(u.bak, low), quantile(u.bak, 0.95))
		x <- c(quantile(transformedU, low), quantile(transformedU, 0.95))
		m <- lm(log2(y) ~ x)
		transformedU <- predict(m, data.frame(x=transformedU))
	}
	attr(transformedU, 'parameter') <- transform.parameter
	attr(transformedU, 'transformFun') <- transFun
	return(transformedU)
}

