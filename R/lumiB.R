lumiB <- function(lumiBatch, method = c('none', 'bg.adjust'), ...) 
{
	if (!(is.function(method)) & !(method %in% c('none', 'bg.adjust'))) {
		print('The method is not supported yet!')
		return(lumiBatch)
	} else if (method == 'none') {
		return(lumiBatch)
	} 
	history.submitted <- as.character(Sys.time())
	if (method == 'bg.adjust') {
		exprs(lumiBatch) <- apply(exprs(lumiBatch), 2, bg.adjust, ...) 
	} else if (is.function(method)) {
		lumiBatch <- method(lumiBatch, ...)
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