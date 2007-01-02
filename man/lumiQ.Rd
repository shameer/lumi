\name{lumiQ}
\alias{lumiQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quality control estimation of the LumiBatch object }
\description{  }

\usage{
lumiQ(x.lumi, logMode = TRUE, sampleRelation.param = list(), detectOutLier.param = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{logMode}{ transform as log2 or not (the function can check whether it is already log transformed.)}
  \item{sampleRelation.param}{ the parameters of sampleRelation estimation, see \code{\link{getSampleRelation}} }
  \item{detectOutLier.param}{ the parameters of sampleRelation estimation, see \code{\link{detectOutlier}} }
}
\details{
  Quality control of a LumiBatch object includes estimating the mean and standard deviation of the chips, chip correlation, the coefficient of variance of each measurement, detectable probe ratio of each chip, density of each chip, sample (chip) relations, detecting outliers of samples (chips). See \code{\link{LumiQC}} for details.
}
\value{
a LumiQC object
}

\author{ Pan Du }

\seealso{\code{\link{LumiQC}}, \code{\link{plot.LumiQC}} }
\examples{
## load example data
data(example.lumi)

## Do quality control estimation
lumi.Q <- lumiQ(example.lumi)

## Plot the results
## plot the pairwise sample correlation
plot(lumi.Q, what='pair')

## see more examples in plot.LumiQC help documents

}
\keyword{methods}
