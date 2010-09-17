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
lumiMethyC(methyLumiM, method = c("quantile", "ssn", "none"), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{method}{ color bias adjustment methods }
  \item{verbose}{a boolean to decide whether to print out some messages }
  \item{\dots}{ other parameters used by corresponding method }
}
\details{

}
\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after background level adjustment.
}
\references{

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
