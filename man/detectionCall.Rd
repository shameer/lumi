\name{detectionCall}
\alias{detectionCall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the detectable probe ratio in each sample }
\description{  }

\usage{
detectionCall(x.lumi, Th = 0.01, type = c('probe', 'sample'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{Th}{ the threshold. By default, when the detection p-value is less than 0.01, we suppose it is detectable. For the old version of BeadStudio output (version 2 or earlier), the threshold will automatically transferred as 1 - Th, because in the old format, value close to 1 is suppose to be detectable. }
  \item{type} { determine to calculate the detection count by probe or by sample }
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

## estimate the detect call (percentage of expressed genes) of each sample
temp <- detectionCall(example.lumi, type='sample')
print(temp)

## estimate the present count of each gene (probe)
temp <- detectionCall(example.lumi, type='probe')
hist(temp)

}
\keyword{ methods }
