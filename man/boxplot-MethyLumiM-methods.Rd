\name{boxplot,MethyLumiM-method}
\docType{methods}
\alias{boxplot,MethyLumiM-method}

\title{boxplot of a MethyLumiM object}
\description{ Creating a \code{\link[hdrcde]{hdr.boxplot}} of the M-value in a MethyLumiM object }

\usage{
\S4method{boxplot}{MethyLumiM}(x, main, logMode = TRUE, ...)
}

\arguments{
\item{x}{ a \code{\link{MethyLumiM-class}} object}
\item{main}{ title of the boxplot }
\item{logMode}{only works when the dataType of x is "Intensity"}
\item{\dots}{ optional arguments to \code{\link[hdrcde]{hdr.boxplot}}.}
}

\details{
 Because the density plot of M-values usually includes two modes, using the traditional boxplot cannot accurately represent the distribution of the data. Here we use 
violin plot to show the density of M-values by samples
}
\seealso{
\code{\link{MethyLumiM-class}}, \code{\link[lattice]{panel.violin}} and \code{\link{boxplot,ExpressionSet-method}}
}
\examples{
## load example data
data(example.lumiMethy)
boxplot(example.lumiMethy)
}
\keyword{methods}
\keyword{hplot}
