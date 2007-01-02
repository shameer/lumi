\name{pairs-methods}
\docType{methods}
\alias{pairs-methods}
\alias{pairs.ExpressionSet}
\alias{pairs,ExpressionSet-method}

\title{Pair plot of an ExpressionSet object}
\description{ Creating \code{\link[graphics]{pairs}} plot of sample intensities in an ExpressionSet object }

\usage{
\S4method{pairs}{ExpressionSet}(x, ..., logMode = TRUE)
}

\arguments{
	\item{x}{ a \code{\link[Biobase]{ExpressionSet}} object}
	\item{\dots}{ optional arguments to \code{\link[graphics]{pairs}}.}
	\item{logMode}{ whether plot the data in log2 scale }
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[graphics]{pairs}}
}
\examples{
## load example data
data(example.lumi)

pairs(example.lumi)
}
\keyword{methods}
\keyword{hplot}
