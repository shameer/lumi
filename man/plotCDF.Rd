\name{plotCDF}
\docType{methods}
\alias{plotCDF}
\alias{plotCDF.ExpressionSet}
\alias{plotCDF,ExpressionSet-method}

\title{plot the cumulative distribution function  of a ExpressionSet object}
\description{ plot the cumulative distribution function of a ExpressionSet object from high to low value or in reverse }

\usage{
 plotCDF(x, reverse=TRUE, logMode=TRUE, xlab = NULL, ylab = "Cumulative density", 
	col=1:dim(x)[2], lwd=1, xlim = NULL, index.highlight = NULL, color.highlight = 2, 
	addLegend = TRUE, main="", ...)
}

\arguments{
	\item{x}{ a numeric or \code{\link[Biobase]{ExpressionSet}} object }
	\item{reverse}{ determine whether plot the CDF from high to low value or not }
	\item{logMode}{ determine whether the CDF plot is based on a log2 scale }
	\item{xlab}{ xlab of the plotCDF plot }
	\item{ylab}{ ylab of the plotCDF plot }
	\item{col}{ line colors of the plotCDF plot }
	\item{lwd}{ line width of plot function }
	\item{xlim}{ parameter of the plot function}
	\item{index.highlight}{ the column index of the highlighted plotCDF curve }
	\item{color.highlight}{ color of highlighted plotCDF curve }
	\item{addLegend}{ whether add legend to the plot or not }
    \item{main}{ title for the plot}
	\item{\dots}{ additional parameters for \code{\link[stats]{plot.ecdf}} function }
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[stats]{ecdf}}
}
\examples{
## load example data
data(example.lumi)

plotCDF(example.lumi)
}
\keyword{methods}
\keyword{hplot}
