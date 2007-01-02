\name{monoSmu}
\alias{monoSmu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Monotonic smooth method }
\description{
  Smoothing the curve with monotonic constraint
}
\usage{
monoSmu(x, y, newX = NULL, nSupport = min(200, length(x)), nKnots = max(length(nSupport)/10, 5), rotate = FALSE, ifPlot = FALSE, xlab = 'x', ylab = 'y', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector represents x values }
  \item{y}{ a vector represents y values }
  \item{newX}{ the new values to be transformed. If not provided, "x" will be used. }
  \item{nSupport}{ downsampled data points }
  \item{nKnots}{ parameter used by \code{\link{monoSpline}} }
  \item{rotate}{ determine whether to rotate the axis with 45 degrees in clockwise }
  \item{ifPlot}{ determine whether to plot intermediate results }
  \item{xlab}{ the xlab of the plot }
  \item{ylab}{ the ylab of the plot }
  \item{\dots}{ parameters used by \code{\link{supsmu}} and \code{\link{plot}} }
}
\details{

}
\value{
Return the transformed "newX" based on the smoothed curve
}
\references{ }
\author{ Simon Lin, Pan Du }

\seealso{ \code{\link{monoSpline}}  }
\examples{

x <- sort(runif(200) * 100)
y <- x + rnorm(200) * 10
temp <- monoSmu(x, y, ifPlot=TRUE)

}
\keyword{ methods }
