\name{gammaFitEM}
\alias{gammaFitEM}

\title{ Estimate the methylation status by fitting a Gamma mixture model using EM algorithm }

\description{
Estimate the methylation status by fitting a two component Gamma mixture model using EM algorithm based on the all M-values of a particular sample 
}
\usage{
gammaFitEM(M, initialFit = NULL, fix.k = NULL, weighted = TRUE, maxIteration = 50, tol = 1e-04, plotMode = FALSE, truncate=FALSE, verbose = FALSE)
}

\arguments{
  \item{M}{ a vector of M-values covering the whole genome }
  \item{initialFit}{ the initial estimation of the gamma parameters returned by .initialGammaEstimation function }
  \item{fix.k}{ the k parameter of  the gamma function which is fixed during estimation }
  \item{weighted}{ determine whether to down-weight the long tails of two component densities beyond their modes }
  \item{maxIteration}{ maximum iterations allowed before converging }
  \item{tol}{ the difference threshold used to determine convergence }
  \item{plotMode}{ determine whether plot the histogram and density plot estimation }
  \item{truncate}{ determine whether to truncate the tails beyond the modes during parameter estimation}
  \item{verbose}{ determine whether plot intermediate messages during iterations }
}
\details{

The assumption of this function is that the M-value distribution is composed of the mixture of two shifted gamma distributions, which are defined as:
dgamma(x-s[1], shape=k[1], scale=theta[1]) and dgamma(s[2]-x, shape=k[2], scale=theta[2]). Here s represents the shift.


NOTE: the methylation status modeling algorithm was developed based on 27K methylation array. It has not been tested for 450K array. Considering 450K array covers both promoter and gene body, the two component Gamma mixture model assumption may not be valid any more.

}
\value{

The return is a list with "gammaFit" class attribute, which includes the following items:
 \item{logLikelihood}{the log-likelihood of the fitting model}
 \item{k}{parameter k of gamma distribution}
 \item{theta}{parameter theta of gamma distribution}
 \item{shift}{parameter shift of gamma distribution}
 \item{proportion}{the proportion of two components (gamma distributions)}
 \item{mode}{the mode positions of the gamma distributions}
 \item{probability}{the estimated methylation status posterior probability of each CpG site}
}

\author{ Pan Du }

\seealso{
  \code{\link{methylationCall}} and \code{\link{plotGammaFit}}
}

\examples{

data(example.lumiMethy)
M <- exprs(example.lumiMethy)
fittedGamma <- gammaFitEM(M[,1], initialFit=NULL, maxIteration=50, tol=0.0001, plotMode=TRUE, verbose=FALSE)

}

\keyword{methods}
