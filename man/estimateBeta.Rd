\name{estimateBeta}
\alias{estimateBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate methylation Beta-value matrix 
}
\description{
Estimate methylation Beta-value matrix from MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities
}
\usage{
estimateBeta(methyLumiM, offset = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{offset}{ An offset value added to the denominator to avoid close to zero intensities  }
}

\details{
Beta-value is ratio between Illumina methylated probe intensity and total probe intensities (sum of methylated and unmethylated probe intensities, see \code{\link{estimateIntensity}}). An offset value added to the denominator to avoid close to zero intensities in the denominator. Beta-value is in the range of 0 and 1. If we assume the probe intensity follows Gamma distribution, then the Beta-value follows a Beta distribution.

}
\value{
A matrix of methylation Beta-value
}
\references{

}
\author{
Pan DU
}


\seealso{
See Also as \code{\link{estimateIntensity}} and \code{\link{estimateM}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
