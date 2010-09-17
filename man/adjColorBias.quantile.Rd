\name{adjColorBias.quantile}
\alias{adjColorBias.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Color bias adjustment of Illumina Infinium methylaton microarrays using smooth quantile normalization
}
\description{
Color bias adjustment of Illumina Infinium methylaton microarrays using smooth quantile normalization \code{\link{smoothQuantileNormalization}}
}
\usage{
adjColorBias.quantile(methyLumiM, refChannel = c("green", "red"), logMode = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ a MethyLumiM object or any eSet object with "methylated" and "unmethylated" data matrix element in the assayData slot }
  \item{refChannel}{ the reference color channel for color bias adjustment }
  \item{logMode}{ whether perform the adjustment in log scale or not }
  \item{\dots}{ other parameters used by \code{\link{smoothQuantileNormalization}}}
}
\details{
Perform color bias adjustment of Illumina Infinium methylaton microarrays. It requires the input methyLumiM object includes the color channel information in the featureData. Basically, there should be a "COLOR_CHANNEL" column in the data.frame returned by pData(featureData(methyLumiM)).

The basic idea of color bias adjustment is to treat it as the normalization between two color channels. It uses smooth quantile normalization \code{\link{smoothQuantileNormalization}} to normalize two color channels.
}
\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after color bias adjustment.
}
\references{

}
\author{
Pan DU
}

\seealso{
 See Also  \code{\link{lumiMethyC}}, \code{\link{smoothQuantileNormalization}} and \code{\link{adjColorBias.ssn}}
}
\examples{
data(example.lumiMethy)
# before adjustment
plotColorBias1D(example.lumiMethy)
lumiMethy.adj = adjColorBias.quantile(example.lumiMethy)
# after adjustment
plotColorBias1D(lumiMethy.adj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
