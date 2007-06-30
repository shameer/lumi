\name{density-methods}
\docType{methods}
\alias{density-methods}
\alias{density.ExpressionSet}
\alias{density,ExpressionSet-method}

\title{Density plot of a ExpressionSet object}
\description{ Creating density plot of sample intensities in a ExpressionSet object. It is equivalent to \link{hist-methods}. }

\usage{
\S4method{density}{ExpressionSet}(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l",
 	index.highlight = NULL, color.highlight = 2, symmetry = NULL, addLegend = TRUE, subset = 5000, ...)
}

\arguments{
	\item{x}{ a \code{\link[Biobase]{ExpressionSet}} object}
	\item{logMode}{ determine whether the density plot is based on a log2 scale }
	\item{xlab}{ xlab of the density plot }
	\item{ylab}{ ylab of the density plot }
	\item{type}{ parameter of plot function }
	\item{index.highlight}{ the column index of the highlighted density curve }
	\item{color.highlight}{ color of highlighted density curve }
	\item{symmetry}{ the boundary position suppose to be symmetric distributed  }
	\item{addLegend}{ whether add legend to the plot or not }
	\item{subset}{subset of rows used to plot. It can be an index vector, or the length of a random subset}
	\item{\dots}{ additional parameters for matplot function }
}

\seealso{
  \code{\link{LumiBatch-class}}, \link{hist-methods}, \code{\link[stats]{density}}
}
\examples{
## load example data
data(example.lumi)

density(example.lumi)
}
\keyword{methods}
\keyword{hplot}
