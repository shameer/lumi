\name{estimateIntensity}
\alias{estimateIntensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate intensity of Illumina methylation data
}
\description{
Estimate intensity of Illumina methylation data, which is the sum of Illumina methylated and unmethylated probe intensities
}
\usage{
estimateIntensity(methyLumiM)
}

\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
}
\details{
	The Intensity basically is the sum of Illumina methylated and unmethylated probe intensities. 
}
\value{
A matrix of methylation Intensity-value
}
\references{

}
\author{
Pan DU
}
\seealso{
See Also as \code{\link{estimateBeta}} and \code{\link{estimateM}}

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
