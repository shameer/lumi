\name{pairs-methods}
\docType{methods}
\alias{pairs-methods}
\alias{pairs.ExpressionSet}
\alias{pairs,ExpressionSet-method}

\title{Pair plot of an ExpressionSet object}
\description{ Creating \code{\link[graphics]{pairs}} plot of sample intensities in an ExpressionSet object }

\usage{
\S4method{pairs}{ExpressionSet}(x, ..., smoothScatter = FALSE, logMode = TRUE, subset = 5000, main = NULL)
}

\arguments{
	\item{x}{ a \code{\link[Biobase]{ExpressionSet}} object}
	\item{\dots}{ optional arguments to \code{\link[graphics]{pairs}}.}
	\item{smoothScatter}{ whether use \code{\link[geneplotter]{smoothScatter}} function to plot points}	
	\item{logMode}{ whether plot the data in log2 scale }
	\item{subset}{subset of rows used to plot. It can be an index vector, or the length of a random subset}
	\item{main}{ title of the plot }	
}

\details{
To increase the plot efficiency, by default, we only plot RANDOMLY selected subset of points (based on parameter "subset"). If users want to plot all the points, they can set the parameter "subset = NULL". When smoothScatter is set as TRUE, the subsetting will be suppressed because \code{\link[geneplotter]{smoothScatter}} function has good plot efficiency for large number of points.
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[graphics]{pairs}}
}
\examples{
## load example data
data(example.lumi)

pairs(example.lumi)

if (require(geneplotter))
	pairs(example.lumi, smoothScatter=TRUE)
}
\keyword{methods}
\keyword{hplot}
