\name{estimateM}
\alias{estimateM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate methylation M-value matrix 
}
\description{
Estimate methylation M-value matrix from MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities
}
\usage{
estimateM(methyLumiM, minValue = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{minValue}{ minimum value of Illumina methylated and unmethylated probe intensities }
}
\details{
M-value is the log2 ratio between Illumina methylated and unmethylated probe intensities. As variations of small intensities can cause big changes in the ratio estimation, the function add a parameter of "minValue". If the probe intensity is less than "minValue", it will be reset as "minValue".
}
\value{
 A matrix of methylation M-value
}
\references{

}
\author{
Pan DU
}

\seealso{
  See Also as \code{\link{estimateBeta}}, \code{\link{estimateIntensity}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
