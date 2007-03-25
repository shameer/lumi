\name{lumiR}
\alias{lumiR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in Illumina expression data }
\description{
  Read in Illumina expression data. We assume the data was saved in a comma or tab separated text file.  
}
\usage{
lumiR(fileName, sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{ fileName of the data file }
  \item{sep}{ the separation character used in the text file.  }
  \item{detectionTh}{ the p-value threshold of determining detectability of the expression. The detection value less than the detectionTh is claimed as detectable. }
  \item{na.rm}{ determine whether to remove NA }
  \item{lib}{ the annotation library used by \code{\link{addNuId2lumi}}  }
}
\details{
The function can automatically determine the separation character if it is Tab or comma. Otherwise, the user should specify the separator manually. If the annotation library is provided, the Illumina Id will be replaced with nuID, which is used as the index Id for the lumi annotation packages.

Current version of lumiR can adaptively read the output of BeadStudio Verson 2 and 3. The format Version 3 made quite a few changes comparing with previous versions. One change is the detection value. It was called detectable when the detection value is close to one for Version 2 format. However, the detection value became a p-value in the Version 3. As a result, the detectionTh is automatically changed based on the version. The detectionTh 0.01 for the Version 3 will be changed as the detectionTh 0.99 for Version 2. Another big change is that Version 3 separately output the control probe (gene) information and a "Samples Table". As a result, the controlData slot in LumiBatch class was added to keep the control probe (gene) information, and a QC slot to keep the quality control information, including the "Sample Table" output by BeadStudio version 3.
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

