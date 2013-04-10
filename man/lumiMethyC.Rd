\name{lumiMethyC}
\alias{lumiMethyC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Color bias adjust of Illumina Infinium methylation data
}
\description{
Color bias adjust of Illumina Infinium methylation data, which is an object in MethyLumiM class.
}
\usage{
lumiMethyC(methyLumiM, method = c("quantile", "ssn", "none"), verbose = TRUE, overwriteBigMatrix=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{method}{ color bias adjustment methods or user provided function, see "details" for more information of user defined function. }
  \item{verbose}{a boolean to decide whether to print out some messages }
  \item{overwriteBigMatrix}{whether to overwrite the result to the BigMatrix data, only valid when the input data is BigMatrix-based}
  \item{\dots}{ other parameters used by corresponding method }
}
\details{
	The first two arguments of the user defined function should be two intensity matrix (pool of methylated and unmethylated probe intensities) of red and green channel respectively. The return of the user defined function should be a list including color adjusted matrix of red and green channel. For example: return(list(red=redData, green=grnData)). "redData" and "grnData" are two color adjusted matrix.
}
\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after background level adjustment.
}

\author{
Pan DU
}

\seealso{
See Also \code{\link{adjColorBias.quantile}} and \code{\link{adjColorBias.ssn}}
}
\examples{
data(example.lumiMethy)
# before adjustment
plotColorBias1D(example.lumiMethy)
# plot in 2D plot of one selected sample
plotColorBias2D(example.lumiMethy, selSample = 1)
lumiMethy.adj = lumiMethyC(example.lumiMethy)
# after adjustment
plotColorBias1D(lumiMethy.adj)
# plot in 2D plot of one selected sample
plotColorBias2D(lumiMethy.adj, selSample = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
