\name{getControlData}
\alias{getControlData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get control probe information }
\description{
  Get control probe information from Bead Studio output or a LumiBatch object.
}
\usage{
getControlData(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the control data can be a LumiBatch object or the Control Probe Profile file outputted by BeadStudio }
}
\details{
  
}
\value{
   return a data.frame with first two columns as "controlType" and "ProbeID". The rest columns are the expression amplitudes for individual samples.
}
\references{  }
\author{ Pan Du }
\seealso{ \code{\link{addControlData2lumi}} }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
