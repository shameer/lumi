\name{lumiB}
\alias{lumiB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Background correction of Illumina data }
\description{
  Background correction of Illumina data
}
\usage{
lumiB(lumiBatch, method = c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiBatch}{ a LumiBatch Object, which can be the return of \code{\link{lumiR}} }
  \item{method}{ the background correction method, it can be any function with a LumiBatch Object as the first argument and return a LumiBatch Object }
  \item{\dots}{ other parameters used by the user provided background correction method  }
}
\details{
  We assume the BeadStudio output data is background corrected. So by default, it will do nothing. The 'bgAdjust' method will estimate the background based on the control probe information, which is kept in the controlData slot of LumiBatch object. The 'forcePositive' method will force all expression values to be positive by adding an offset (minus minimum value plus one), it does nothing if all expression values are positive. The purpose of this is to avoid NA when do logarithm transformation. 'none' does not but return the LumiBatch object. 'bgAdjust.affy' will call the \code{\link[affy]{bg.adjust}} function in affy package. User can also provide their own function with a LumiBatch Object as the first argument and return a LumiBatch Object with background corrected.

Thanks Kevin Coombes (M.D. Anderson Cancer Center) suggested adding this function.
}
\value{
 return a LumiBatch Object with background corrected.
}

\author{ Pan Du, Kevin Coombes }
\seealso{ \code{\link{bgAdjust}}, \code{\link{lumiExpresso}} }
\examples{
## load example data
data(example.lumi)

## Do the default background correction method
lumi.B <- lumiB(example.lumi, method='bgAdjust', probs=0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
