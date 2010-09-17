\name{lumiMethyB}
\alias{lumiMethyB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjust background level of Illumina Infinium methylation data
}
\description{
Adjust background level of Illumina Infinium methylation data, which is an object in MethyLumiM class.
}
\usage{
lumiMethyB(methyLumiM, method = c("bgAdjust2C", "forcePositive", "none"), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{method}{ background adjustment methods }
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
See Also \code{\link{bgAdjustMethylation}} and \code{\link{estimateMethylationBG}}
}
\examples{
data(example.lumiMethy)
lumiMethy.bgAdj = lumiB(example.lumiMethy)
attr(lumiMethy.bgAdj, "EstimatedBG")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
