`inverseVST` <-
function(x, fun=c('asinh', 'log'), parameter) {
	if (missing(x) | missing(fun) | missing(parameter)) {
		stop('Please provide all input parameters!')
	}
	fun <- match.arg(fun)
	## parameter is in the format of: "a, b, g, Intercept"
	a <- parameter[1]
	b <- parameter[2]
	g <- parameter[3]
	intercept <- parameter[4]
	if (fun == 'asinh') {
		## Transform function h(x) = g * asinh(a + b * x)
		## transFun <- g * asinh(a + b * x) * m$coef[2] + m$coef[1]
		inv <- (sinh((x - intercept)/g) - a)/b
	} else if (fun == 'log') {
		## Transform function h(x) = g * log(a + b * x)
		## transFun <- g * log(a + b * x) * m$coef[2] + m$coef[1]
		inv <- (exp((x - intercept)/g) - a)/b
	}
	return(inv)
}

