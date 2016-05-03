\name{normalizeMethylation.ssn}
\alias{normalizeMethylation.ssn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Shift and scaling normalization of Illumina Infinium methylation data at probe level
}
\description{
Shift and scaling normalization of Illumina Infinium methylation data at probe level
}
\usage{
normalizeMethylation.ssn(methyLumiM, separateColor = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{separateColor}{determine whether separately normalize two color channels}
}

\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after background level adjustment.
}

\author{
Pan DU
}
\seealso{
 See Also \code{\link{lumiMethyN}}, and \code{\link{normalizeMethylation.quantile}}
}
\examples{
data(example.lumiMethy)
lumiMethy.norm = normalizeMethylation.ssn(example.lumiMethy)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
