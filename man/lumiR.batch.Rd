\name{lumiR.batch}
\alias{lumiR.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read BeadStudio output files in batch}
\description{
  Read BeadStudio output files in batch and combine them as a single LumiBatch object
}
\usage{
lumiR.batch(fileList, convertNuID = TRUE, lib.mapping = NULL, detectionTh = 0.01, QC = TRUE, transform = c('none', 'vst', 'log2', 'cubicRoot'), sampleInfoFile = NULL, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileList}{ a vector of file names or a directory keeping the data files in the format of .csv }
  \item{convertNuID}{ determine whether convert the probe identifier as nuID }
  \item{lib.mapping}{ same as \code{\link{lumiR}} parameter lib.mapping (optional) }
  \item{detectionTh}{ the p-value threshold of determining detectability of the expression. See more details in \code{\link{lumiQ}}  }
  \item{QC}{ determine whether to do quality control assessment after read in the data.}
  \item{transform}{ determine whether to do transform after input each file }
  \item{sampleInfoFile}{ a Tab-separated text file or a data.frame keeping the sample information (optional) }
  \item{verbose}{ a boolean to decide whether to print out some messages }
  \item{\dots}{ other parameters used by \code{\link{lumiR}} }
}
\details{
The function basically call lumiR for individual files and then combine the returns. The sampleInfoFile parameter is optional. It provides the sample information (for phenoData slot in LumiBatch object), it is a Tab-separated text file. ID column is required. It represents sample ID, which is defined based on the column names of BeadStudio output file. For example, sample ID of column "1881436070\_A\_STA.AVG\_Signal" is "1881436070\_A\_STA". The sample ID column can also be found in the "Samples Table.txt" file output by BeadStudio. Another "Label" column (if provided) will be used as the sampleNames of LumiBatch object. All information of sampleInfoFile will be directly added in the phenoData slot in LumiBatch object.

To save memory space in the case of reading large data set, we can do transformation using lumiT function right after input the data, and the information like se.exprs, beadNum will be removed from the LumiBatch object after transformation.
}
\value{
A LumiBatch object which combines the individual LumiBatch object corresponding to each file
}
\author{ Pan Du }

\seealso{ \code{\link{lumiR}} }
\examples{
## fileList <- c('file1.csv', 'file2.cvs')
## x.lumi <- lumiR.batch(fileList, sampleInfoFile='sampleInfo.txt')

}
\keyword{ methods }
