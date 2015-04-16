\name{adjColorBias.ssn}
\alias{adjColorBias.ssn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Color bias adjustment of Illumina Infinium methylaton microarrays using simple shift and scaling normalization
}
\description{
Color bias adjustment of Illumina Infinium methylaton microarrays using simple shift and scaling normalization
}
\usage{
adjColorBias.ssn(methyLumiM, refChannel = c("green", "red", "mean"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{methyLumiM}{ a MethyLumiM object or any eSet object with "methylated" and "unmethylated" data matrix element in the assayData slot }
\item{refChannel}{ the reference color channel for color bias adjustment }
}
\details{
Perform color bias adjustment of Illumina Infinium methylaton microarrays. It requires the input methyLumiM object includes the color channel information in the featureData. Basically, there should be a "COLOR_CHANNEL" column in the data.frame returned by pData(featureData(methyLumiM)).

The basic idea of color bias adjustment is to treat it as the normalization between two color channels. It uses simple scaling normalization to normalize two color channels. The background levels are estimated using function \code{\link{estimateMethylationBG}}.
}
\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after color bias adjustment.
}
\author{
Pan DU
}

\seealso{
See Also  \code{\link{lumiMethyC}}, \code{\link{estimateMethylationBG}} and \code{\link{adjColorBias.quantile}}
}
\examples{
data(example.lumiMethy)
# before adjustment
plotColorBias1D(example.lumiMethy)
lumiMethy.adj = adjColorBias.ssn(example.lumiMethy)
# after adjustment
plotColorBias1D(lumiMethy.adj)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{method}

