\name{methylationCall}
\alias{methylationCall}

\title{ Estimated methylation call }

\description{
Estimated methylation call based on the fitting results of \code{\link{gammaFitEM}}
}

\usage{
methylationCall(x, threshold = 0.95, ...)
}

\arguments{
	\item{x}{ a vector of M-values covering the whole genome or a "gammaFit" class object returned by \code{\link{gammaFitEM}} }
	\item{threshold}{the probability threshold to make a methylation call. The threshold can be a vector of two: unmethylation threshold and methylation threshold}
	\item{\dots}{ other parameters used by \code{\link{gammaFitEM}} }
}
\details{
	Retrieve the probability element returned by \code{\link{gammaFitEM}}, and convert it as three status calls based on probability threshold
	
	NOTE: the methylation status modeling algorithm was developed based on 27K methylation array. It has not been tested for 450K array. Considering 450K array covers both promoter and gene body, the two component Gamma mixture model assumption may not be valid any more.

}
\value{
A vector of three methylation status: "Unmethy" (unmethylation posterior probability > unmethylation threshold), "Methy" (methylation posterior probability > methylation threshold), or "Margin". The sum of unmethylation posterior probability and methylation posterior probability equals one. The methylation probability is returned as an attribute of "probability".
}

\author{ Pan DU }

\seealso{
\code{\link{gammaFitEM}}
}

\examples{
data(example.lumiMethy)
M <- exprs(example.lumiMethy)
fittedGamma <- gammaFitEM(M[,1], initialFit=NULL, maxIteration=50, tol=0.0001, plotMode=TRUE, verbose=FALSE)
methyCall <- methylationCall(fittedGamma)
table(methyCall)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
