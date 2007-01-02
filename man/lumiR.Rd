\name{lumiR}
\alias{lumiR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in Illumina expression data }
\description{
  Read in Illumina expression data. We assume the data was saved in a comma or tab separated text file.  
}
\usage{
lumiR(fileName, sep = NULL, detectionTh = 0.99, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{ fileName of the data file }
  \item{sep}{ the separation character used in the text file.  }
  \item{detectionTh}{ the threshold of determining detectability of the expression. The detection value larger than the detectionTh is claimed as detectable. }
  \item{na.rm}{ determine whether to remove NA }
}
\details{
The function can automatically determine the separation character if it is Tab or comma. Otherwise, the user should specify the separator manually.
}
\value{
return a LumiBatch object
}

\author{ Simon Lin, Pan Du }
\seealso{ \code{\link{LumiBatch}}, \code{\link{addNuId2lumi}} }

\examples{
## specify the file name
# fileName <- 'Barnes_gene_profile.txt' # Not Run
## load the data
# x.lumi <- lumiR(fileName)

}
\keyword{ IO }

