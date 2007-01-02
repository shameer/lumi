\name{detectionCall}
\alias{detectionCall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the detectable probe ratio in each sample }
\description{  }

\usage{
detectionCall(x.lumi, Th = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{Th}{ the threshold. By default, if the detection value is larger than 0.99, we suppose it is detectable.  }
}

\value{
Return the detectable probe ratio of each sample
}
\references{  }
\author{ Pan Du }

\seealso{  \code{\link{lumiQ}} }
\examples{
## load example data
data(example.lumi)
## load example data
data(example.lumi)

## estimate the detect call (percentage of expressed genes)
temp <- detectionCall(example.lumi)
print(temp)

}
\keyword{ methods }
