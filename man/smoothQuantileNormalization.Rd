\name{smoothQuantileNormalization}
\alias{smoothQuantileNormalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth quantile normalization
}
\description{
Smooth quantile normalization with a reference sample
}
\usage{
smoothQuantileNormalization(dataMatrix, ref = NULL, adjData=NULL, logMode = TRUE, bandwidth = NULL, degree = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrix}{ a matrix of microarray intensity data}
  \item{ref}{ a vector of reference sample intensity, which matches the dataMatrix }
  \item{adjData}{ data to be adjusted based on the ref and dataMatrix distribution }
  \item{logMode}{ whether perform the analysis in log2 scale}
  \item{bandwidth}{ a parameter used by \code{\link[KernSmooth]{locpoly}} }
  \item{degree}{ a parameter used by \code{\link[KernSmooth]{locpoly}} }
  \item{\dots}{other parameters used by \code{\link[KernSmooth]{locpoly}} }
}
\value{
a data matrix with intensity normalized.
}
\author{
Pan DU
}

\seealso{
 See Also \code{\link{adjColorBias.quantile}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
