\name{addControlData2methyLumiM}
\alias{addControlData2methyLumiM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add methylation control data to a MethyLumiM object}
\description{
Add methylation control data to a MethyLumiM object
}
\usage{
addControlData2methyLumiM(controlData, methyLumiM, checkConsistency = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{a methylation control data file (output by GenomeStudio), or a MethyLumiQC object}
  \item{methyLumiM}{a MethyLumiM object to add control data}
  \item{checkConsistency}{whether to check the sample names consistency between methyLumiM and controlData}
  \item{\dots}{other parameters for reading controlData}
}
\details{
This function aims to add the controlData (MethyLumiQC object) to the controlData slot of a methyLumiM object
For control data, methylated data matrix in assayData slot corresponds to green channel, and unmethylated data matrix in assayData slot corresponds to red channel.

}
\value{
Return the methyLumiM object with the controlData added
}
\author{ Pan DU}
\seealso{
\code{\link{lumiMethyR}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
