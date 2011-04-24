\name{plotDensity}
\docType{methods}
\alias{plotDensity}
\alias{plotDensity.ExpressionSet}

\title{plot the density distribution}
\description{ plot the density distribution of a dataMatrix or ExpressionSet object }

\usage{
 plotDensity(dataMatrix, logMode=TRUE, addLegend=TRUE, legendPos="topright", subset = NULL, ...)
}

\arguments{
	\item{dataMatrix}{ a data matrix or \code{\link[Biobase]{ExpressionSet}} object }
	\item{logMode}{ determine whether the CDF plot is based on a log2 scale }
	\item{addLegend}{ whether add legend to the plot or not }
	\item{legendPos}{ the position of the legend. If its length is two, then it specifies the x and y location of legend. }
	\item{subset}{plot subset of randomly selected rows. All data will be plotted if it is NULL. }
	\item{\dots}{ additional parameters for \code{\link[graphics]{matplot}} function }
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[stats]{density}}
}
\examples{
## load example data
data(example.lumi)

plotDensity(example.lumi)
}
\keyword{methods}
\keyword{hplot}
