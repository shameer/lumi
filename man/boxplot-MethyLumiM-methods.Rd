\name{boxplot,MethyLumiM-method}
\docType{methods}
\alias{boxplot,MethyLumiM-method}

\title{boxplot of a MethyLumiM object}
\description{ Creating a \code{\link[hdrcde]{hdr.boxplot}} of the M-value in a MethyLumiM object }

\usage{
\S4method{boxplot}{MethyLumiM}(x, main, prob=c(seq(10,90, by=10), 95), col=gray(rev(seq(prob)/length(prob))), logMode = TRUE, ...)
}

\arguments{
\item{x}{ a \code{\link{MethyLumiM-class}} object}
\item{main}{ title of the boxplot }
\item{prob}{ Probability coverage required for HDRs (highest density regions), see \code{\link[hdrcde]{hdr.boxplot}}}
\item{col}{ the plot color to represent the data in different quantile levels}
\item{logMode}{only works when the dataType of x is "Intensity"}
\item{\dots}{ optional arguments to \code{\link[hdrcde]{hdr.boxplot}}.}
}

\details{
 Because the density plot of M-values usually includes two modes, using the traditional boxplot cannot accurately represent the distribution of the data. Here we use different color levels to represent the M-values in different quantile levels. The M-values locating outside the range of the maximum probability specified in "prob" are plotted as dots, which is similar with the outliers in the regular boxplot.
}
\seealso{
\code{\link{MethyLumiM-class}}, \code{\link[hdrcde]{hdr.boxplot}} and \code{\link{boxplot,ExpressionSet-method}}
}
\examples{
## load example data
data(example.lumiMethy)
boxplot(example.lumiMethy)
}
\keyword{methods}
\keyword{hplot}
