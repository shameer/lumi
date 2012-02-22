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
estimateIntensity(methyLumiM, returnType=c("ExpressionSet", "matrix"))
}

\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{returnType}{determine whether return an ExpressionSet or matrix object}
}
\details{
	The Intensity basically is the sum of Illumina methylated and unmethylated probe intensities. 
}
\value{
An ExpressionSet or matrix object of methylation Intensity-value
}
\author{
Pan DU
}
\seealso{
See Also as \code{\link{estimateBeta}} and \code{\link{estimateM}}

}
\examples{
data(example.lumiMethy)
methyLumiIntensity = estimateIntensity(example.lumiMethy)
boxplot(methyLumiIntensity)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
