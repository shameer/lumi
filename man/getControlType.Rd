\name{getControlType}
\alias{getControlType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get the types of the control probes }
\description{
  Get the types of the control probes, which is in the first column of the controlData data.frame for LumiBatch objects. For methylation data, it is the return of controlTypes function 
}
\usage{
getControlType(controlData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{ a LumiBatch object including control data, a control data data.frame, or a MethyLumiQC object for methylation data }
}
\value{
	return the unique type of control probe type.
}
\author{ Pan Du }

\seealso{  \code{\link{addControlData2lumi}}, \code{\link[methylumi]{controlTypes}} for methylation data }
\examples{
	controlFile <- system.file('doc', 'Control_Probe_Profile.txt', package='lumi')
	## return a data.frame
	controlData <- getControlData(controlFile)
	getControlType(controlData)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
