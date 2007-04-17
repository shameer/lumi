\name{boxplot-methods}
\docType{methods}
\alias{boxplot-methods}
\alias{boxplot.ExpressionSet}
\alias{boxplot,ExpressionSet-method}

\title{boxplot of a ExpressionSet object}
\description{ Creating \code{\link[graphics]{boxplot}} of sample intensities in a ExpressionSet object }

\usage{
\S4method{boxplot}{ExpressionSet}(x, range = 0, main, logMode = TRUE, subset = 5000, ...)
}

\arguments{
\item{x}{ a \code{\link[Biobase]{ExpressionSet}} object}
\item{range}{ parameter of \code{\link[graphics]{boxplot}} }
\item{main}{ title of the boxplot }
\item{logMode}{ whether plot the data in log2 scale or not}
\item{subset}{subset of rows used to plot. It can be an index vector, or the length of a random subset}
\item{\dots}{ optional arguments to \code{\link[graphics]{boxplot}}.}
}

\seealso{
\code{\link{LumiBatch-class}}, \code{\link[graphics]{boxplot}}
}
\examples{
## load example data
data(example.lumi)

boxplot(example.lumi)
}
\keyword{methods}
\keyword{hplot}
