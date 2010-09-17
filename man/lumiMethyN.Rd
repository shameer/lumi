\name{lumiMethyN}
\alias{lumiMethyN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize the Illumina Infinium methylation data
}
\description{
Normalize the Illumina Infinium methylation data, which is an object in MethyLumiM class.
}
\usage{
lumiMethyN(methyLumiM, method = c("ssn", "quantile", "none"), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
 \item{method}{ supported normalization methods }
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
 See Also \code{\link{normalizeMethylation.ssn}} and \code{\link{normalizeMethylation.quantile}}
}
\examples{
data(example.lumiMethy)
lumiMethy.norm = lumiN(example.lumiMethy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
