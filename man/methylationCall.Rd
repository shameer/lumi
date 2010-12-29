\name{methylationCall}
\alias{methylationCall}

\title{ retrieve the estimated methylation call }

\description{
retrieve the estimated methylation call from the fitting results of \code{\link{gammaFitEM}}
}

\usage{
methylationCall(x, k = NULL, theta = NULL, shift = NULL, proportion = NULL, threshold = 0.99, truncate = TRUE)
}

\arguments{
	\item{x}{ a vector of M-values covering the whole genome or a "gammaFit" class object returned by \code{\link{gammaFitEM}} }
	\item{k}{parameter k of gamma distribution}
	\item{theta}{parameter theta of gamma distribution}
	\item{shift}{parameter shift of gamma distribution}
	\item{proportion}{the proportion of two components (gamma distributions)}
	\item{threshold}{the probability threshold to make a binary methylation call}
	\item{truncate}{ determine whether truncating the density estimation at the extreme side }
}
\details{
	Retrieve the probability element returned by \code{\link{gammaFitEM}}, and convert it as binary calls based on probability threshold
}
\value{
Return a list which includes the following items:
	\item{methylationCall}{the binary call of methylation status}
	\item{probability}{the methylation call probability}

}
\references{

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
