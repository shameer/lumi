\name{MAplot-methods}
\docType{methods}
\alias{MAplot-methods}
\alias{MAplot.ExpressionSet}
\alias{MAplot,ExpressionSet-method}

\title{MAplot of a ExpressionSet object}
\description{ Creating pairwise MAplot of sample intensities in a ExpressionSet object }

\usage{
\S4method{MAplot}{ExpressionSet}(object, ..., smoothScatter = FALSE, logMode = TRUE, subset = 5000, main = NULL)
}

\arguments{
\item{object}{ an \code{\link[Biobase]{ExpressionSet}} object}
\item{\dots}{ optional arguments to \code{\link[affy]{MAplot}}.}
\item{smoothScatter}{ whether use \code{\link[geneplotter]{smoothScatter}} function to plot points}	
\item{logMode}{ whether plot the data in log2 scale or not}	
\item{subset}{ subset of rows used to plot. It can be an index vector, or the length of a random subset}
\item{main}{ title of the plot }	
}

\details{
To increase the plot efficiency, by default, we only plot subset of points (based on parameter "subset"). If users want to plot all the points, they can set the parameter "subset = NULL". When smoothScatter is set as TRUE, the subsetting will be suppressed because \code{\link[geneplotter]{smoothScatter}} function has good plot efficiency for large number of points.
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[affy]{MAplot}}
}
\examples{
## load example data
data(example.lumi)

MAplot(example.lumi)

if (require(geneplotter))
	MAplot(example.lumi, smoothScatter=TRUE)
}
\keyword{methods}
\keyword{hplot}
