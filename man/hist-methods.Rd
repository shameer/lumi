\name{hist-methods}
\docType{methods}
\alias{hist-methods}
\alias{hist.ExpressionSet}
\alias{hist,ExpressionSet-method}

\title{Density plot of a ExpressionSet object}
\description{ Creating density plot of sample intensities in a ExpressionSet object }

\usage{
\S4method{hist}{ExpressionSet}(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l",
 	index.highlight = NULL, color.highlight = 2, symmetry = NULL, addLegend = TRUE, ...)
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
	\item{\dots}{ additional parameters for matplot function }
}

\seealso{
  \code{\link{LumiBatch-class}}, \code{\link[graphics]{hist}}
}
\examples{
## load example data
data(example.lumi)

hist(example.lumi)
}
\keyword{methods}
\keyword{hplot}
