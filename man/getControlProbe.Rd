\name{getControlProbe}
\alias{getControlProbe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get the control probe Ids}
\description{
  Get the control probe Ids corresponding to the control probe type provided. The control probe ids are kept in the second column of controlData data.frame.
}
\usage{
getControlProbe(controlData, type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlData}{ a LumiBatch object including control data or a control data data.frame }
  \item{type}{ the type of control probe, which can be get by using \code{\link{getControlType}} function }
}
\details{
  
}
\value{
returns the corresponding probe Ids for the control type.
}
\references{  }
\author{ Pan Du }
\seealso{  \code{\link{addControlData2lumi}} }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
