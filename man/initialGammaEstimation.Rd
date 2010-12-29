\name{initialGammaEstimation}
\alias{initialGammaEstimation}

\title{Initial gamma parameter estimation for function \code{\link{gammaFitEM}} }

\description{
Initial gamma parameter estimation for function \code{\link{gammaFitEM}} 
}

\usage{
initialGammaEstimation(x, initialFit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of M-values covering the whole genome or a "gammaFit" class object returned by \code{\link{gammaFitEM}} }
  \item{initialFit}{ previous initial parameter estimation}
}
\details{
Estimate the initial parameters based on histogram distribution of the data
}
\value{
The return is a list with "gammaFit" class attribute, which includes the following items:
 \item{k}{parameter k of gamma distribution}
 \item{theta}{parameter theta of gamma distribution}
 \item{shift}{parameter shift of gamma distribution}
 \item{proportion}{the proportion of two components (gamma distributions)}
 \item{mode}{the mode positions of the gamma distributions}

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
initialEstimate <- initialGammaEstimation(M[,1])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

