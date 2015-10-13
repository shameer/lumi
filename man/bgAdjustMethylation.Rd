\name{bgAdjustMethylation}
\alias{bgAdjustMethylation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate and adjust the background levels of Illumina Infinium methylaton microarrays
}
\description{
Estimate and adjust the background levels of Illumina Infinium methylaton microarrays
}
\usage{
bgAdjustMethylation(methyLumiM, separateColor = FALSE, targetBGLevel = 300, negPercTh = 0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object or any eSet object with "methylated" and "unmethylated" data matrix element in the assayData slot }
  \item{separateColor}{determine whether separately process two color channels}
  \item{targetBGLevel}{ adjust background level to a non-zero target background level }
  \item{negPercTh}{the threshold of the percentage of negative values after subtract estimated background levels. A warning will be given if too many probes having intensities lower than background levels.}
}
\details{
The estimation of background level of Infinium methylaton microarray is based on the assumption that the lots of CpG sites are unmethylated, which results in a density mode of the intensities measured by methylated probes. The position of this mode represents the background level.
}
\value{
Return an object (same class as input methyLumiM) with updated "methylated" and "unmethylated" data matrix after background level adjustment. The estimated background level was kept in the attribute, "EstimatedBG", of the returned methyLumiM object.
}
\author{
Pan DU
}

\seealso{
 See Also \code{\link{lumiMethyB}} and \code{\link{estimateMethylationBG}}
}
\examples{
data(example.lumiMethy)
lumiMethy.bgAdj = bgAdjustMethylation(example.lumiMethy)
attr(lumiMethy.bgAdj, "EstimatedBG")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
