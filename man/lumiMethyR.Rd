\name{lumiMethyR}
\alias{lumiMethyR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reading Illumina methylation microarray data
}
\description{
This function is a wrap of \code{\link[methylumi]{methylumiR}} function in methylumi package.
}
\usage{
lumiMethyR(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ parameters of \code{\link[methylumi]{methylumiR}} function in methylumi package }
}
\details{
This function is a wrap of \code{\link[methylumi]{methylumiR}} function in methylumi package. It will coerce the returned object as MethyLumiM class.
}
\value{ return a MethyLumiM object }

\author{
Pan Du
}

\seealso{
  See Also \code{\link[methylumi]{methylumiR}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
