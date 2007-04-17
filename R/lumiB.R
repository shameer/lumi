lumiB <- function(lumiBatch, method = c('forcePositive', 'none', 'bg.adjust'), ...) 
{
	if (!(is.function(method)) & !(method %in% c('forcePositive', 'none', 'bg.adjust'))) {
		print('The method is not supported yet!')
		return(lumiBatch)
	} else if (method == 'none') {
		return(lumiBatch)
	} else if (method == 'forcePositive') {
		if (min(exprs(lumiBatch)) > 0)  return(lumiBatch)
	}
	
	history.submitted <- as.character(Sys.time())
	if (method == 'bg.adjust') {
		exprs(lumiBatch) <- apply(exprs(lumiBatch), 2, bg.adjust, ...) 
	} else if (method == 'forcePositive') {
		exprs(lumiBatch) <- exprs(lumiBatch) - min(exprs(lumiBatch)) + 1
	} else if (is.function(method)) {
		lumiBatch <- method(lumiBatch, ...)
	} else {
		print('The method is not supported!')
		return(lumiBatch)
	}

	# history tracking
	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(lumiB)))

	if (is.null(lumiBatch@history$lumiVersion)) lumiBatch@history$lumiVersion <- rep(NA, nrow(lumiBatch@history))
	lumiVersion <- packageDescription('lumi')$Version
	lumiBatch@history<- rbind(lumiBatch@history, data.frame(submitted=history.submitted, 
			finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	return(lumiBatch)
}