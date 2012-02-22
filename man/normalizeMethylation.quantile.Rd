\name{normalizeMethylation.quantile}
\alias{normalizeMethylation.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile normalization of Illumina Infinium methylation data at probe level
}
\description{
Quantile normalization of Illumina Infinium methylation data at probe level. Input data is a MethyLumiM object
}
\usage{
normalizeMethylation.quantile(methyLumiM, separateColor = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{separateColor}{determine whether separately normalize two color channels}
  \item{\dots}{ other parameters used by normalize.quantiles.robust method }
}
\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after background level adjustment.
}
\author{
Pan DU
}
\seealso{
 See Also \code{\link{lumiMethyN}}, and \code{\link{normalizeMethylation.ssn}}
}
\examples{
data(example.lumiMethy)
lumiMethy.norm = normalizeMethylation.quantile(example.lumiMethy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
