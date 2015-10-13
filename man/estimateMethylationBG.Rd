\name{estimateMethylationBG}
\alias{estimateMethylationBG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the background levels of Illumina Infinium methylaton microarrays
}
\description{
Estimate the background levels of Illumina Infinium methylaton microarrays. It is called by function \code{\link{bgAdjustMethylation}}
}
\usage{
estimateMethylationBG(methyLumiM, separateColor = FALSE, nbin = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object or any eSet object with "methylated" and "unmethylated" data matrix element in the assayData slot }
  \item{separateColor}{determine whether to separately process two color channels}
  \item{nbin}{ the number of bins in the histogram used to estimate the mode position of the density }
}
\details{
When the controlData includes the negative control probe information, the background estimation will be the median of the negative control probes. Red and Green color channels will be estimated separately. 

In the case the negative control data is not available, the background will be estimated based on the mode positions of unmethylated or methylated distribution (the smaller one). The assumption is that the lots of CpG sites are unmethylated, which results in a density mode of the intensities measured by methylated probes. The position of this mode represents the background level.
}
\value{
a vector of estimated background levels for individual samples.
}
\author{
Pan DU
}

\seealso{
See Also \code{\link{lumiMethyB}} and \code{\link{bgAdjustMethylation}}
}
\examples{
data(example.lumiMethy)
estimatedBG = estimateMethylationBG(example.lumiMethy)
estimatedBG
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
