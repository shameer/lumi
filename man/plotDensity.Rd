\name{plotDensity}
\alias{plotDensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density plot of the columns in a matrix }
\description{  }

\usage{
plotDensity(x, logMode=TRUE, xlab = NULL, ylab = "density", type = "l", index.highlight = NULL, color.highlight = 2, symmetry = NULL, addLegend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a matrix with columns as the profiles }
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
\value{
Invisibly return the x and y vectors of density plot
}
\author{Pan Du }
\seealso{ }

\note{ 
This function is modified based on the plotDensity function in affy package.
}

\examples{
## load example data
data(example.lumi)

## plot the density plot with the third chip highlighted
plotDensity(example.lumi, index.highlight=3)

}
\keyword{methods}
