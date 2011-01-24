\name{plotGammaFit}
\alias{plotGammaFit}

\title{ plot the fitting results of \code{\link{gammaFitEM}} }

\description{
plot the fitting results of \code{\link{gammaFitEM}}. It basically plot the histogram of M-values plus the estimated gamma density plots and their mixture.
}

\usage{
plotGammaFit(x, gammaFit = NULL, k = NULL, theta = NULL, shift = NULL, proportion = NULL, plotType = c("histogram", "density"), ...)
}

\arguments{
  \item{x}{ a vector of M-values covering the whole genome }
  \item{gammaFit}{ a "gammaFit" class object returned by \code{\link{gammaFitEM}} }
  \item{k}{parameter k of gamma distribution}
  \item{theta}{parameter theta of gamma distribution}
  \item{shift}{parameter shift of gamma distribution}
  \item{proportion}{the proportion of two components (gamma distributions)}
  \item{plotType}{ determine the way to show the distribution of the input data, either histogram or density plot }
  \item{\dots}{ Other parameters used by \code{\link[graphics]{hist}} or \code{\link[graphics]{plot}} (for "density" plotType) function. }
}

\details{
 This function is to visualize the fitting results, which helps us understand how well the fitting is.
}
\value{
Invisibly return TRUE if the plot is successful.
}
\references{

}
\author{ Pan DU }

\seealso{
\code{\link{gammaFitEM}}
}

\examples{

data(example.lumiMethy)
M <- exprs(example.lumiMethy)
fittedGamma <- gammaFitEM(M[,1], initialFit=NULL, maxIteration=50, tol=0.0001, plotMode=FALSE, verbose=FALSE)
plotGammaFit(M[,1], gammaFit=fittedGamma)

}

\keyword{methods}
