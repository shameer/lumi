\name{lumiB}
\alias{lumiB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Background correction of Illumina data }
\description{
  Background correction of Illumina data
}
\usage{
lumiB(lumiBatch, method = c('none', 'bg.adjust'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiBatch}{ a LumiBatch Object, which can be the return of \code{\link{lumiR}} }
  \item{method}{ the background correction method, it can be any function with a LumiBatch Object as the first argument and return a LumiBatch Object }
  \item{\dots}{ other parameters used by the user provided background correction method  }
}
\details{
  The default background correction method does nothing but return the input LumiBatch object. 'bg.adjust' method is based on the \code{\link[affy]{bg.adjust}} function in affy package. User can also provide their own function with a LumiBatch Object as the first argument and return a LumiBatch Object with background corrected.

Thanks Kevin Coombes (M.D. Anderson Cancer Center) suggested adding this function.
}
\value{
 return a LumiBatch Object with background corrected.
}

\author{ Pan Du, Kevin Coombes }
\seealso{ \code{\link{lumiExpresso}} }
\examples{
## load example data
data(example.lumi)

## Do the default background correction method
lumi.B <- lumiB(example.lumi)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
