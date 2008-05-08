\name{MAplot-methods}
\docType{methods}
\alias{MAplot-methods}
\alias{MAplot.ExpressionSet}
\alias{MAplot,ExpressionSet-method}

\title{MAplot of a ExpressionSet object}
\description{ Creating pairwise MAplot of sample intensities in a ExpressionSet object }

\usage{
\S4method{MAplot}{ExpressionSet}(object, ..., logMode=TRUE, subset = 5000, seed = 123)
}

\arguments{
\item{object}{ an \code{\link[Biobase]{ExpressionSet}} object}
\item{\dots}{ optional arguments to \code{\link[affy]{MAplot}}.}
\item{logMode}{ whether plot the data in log2 scale or not}	
\item{subset}{ subset of rows used to plot. It can be an index vector, or the length of a random subset}
\item{seed}{ the random seed for random subset}
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[affy]{MAplot}}
}
\examples{
## load example data
data(example.lumi)

MAplot(example.lumi)
}
\keyword{methods}
\keyword{hplot}
