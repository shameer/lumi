\name{lumiMethyStatus}
\alias{lumiMethyStatus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the methylation status of individual methylation sites
}
\description{
Estimate the methylation status of individual methylation sites by fitting a two component Gamma mixture model for each sample
}
\usage{
lumiMethyStatus(methyLumiM, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM class object}
  \item{\dots}{Other parameters used by \code{\link{methylationCall}} }
}
\details{
This function calls \code{\link{methylationCall}} and returns the methylation status of individual methylation sites. The methylation status includes: "Unmethy" (unmethylation probability > unmethylation threshold), "Methy" (methylation probability > methylation threshold), or "Margin". The methylation probability is returned as an attribute of "probability".
}
\value{ return a methylation status matrix with "probability" attribute }

\author{
Pan Du
}

\seealso{
  See Also \code{\link{methylationCall}} and \code{\link{gammaFitEM}}
}
\examples{
data(example.lumiMethy)
methyCall <- lumiMethyStatus(example.lumiMethy)
head(methyCall)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
