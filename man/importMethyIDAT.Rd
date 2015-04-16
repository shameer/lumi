\name{importMethyIDAT}
\alias{importMethyIDAT}
\title{
Import Illumina methylation .idat files as an MethyLumiM object
}
\description{
Import Illumina methylation .idat files as an MethyLumiM object. An extension of \code{\link[methylumi]{lumIDAT}} function
}
\usage{
importMethyIDAT(sampleInfo, dataPath = getwd(), lib = NULL, bigMatrix=FALSE, dir.bigMatrix='.', savePrefix.bigMatrix, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampleInfo}{A data.frame of sample information or a character vector of barcodes.}
  \item{dataPath}{The path of .idat files}
  \item{lib}{Annotation library}
  \item{bigMatrix}{whether to save the data as BigMatrix (designed for very large dataset)}
  \item{dir.bigMatrix}{the parent directory to save the BigMatrix data files}
  \item{savePrefix.bigMatrix}{the folder name prefix of the directory to save the BigMatrix data files. The fold name will be like this: paste(savePrefix.bigMatrix, '_bigmat', sep='') }
  \item{\dots}{other parameters used by \code{\link[methylumi]{lumIDAT}} function}
}
\details{
This function is an extension of \code{\link[methylumi]{lumIDAT}}. It adds sample information and probe annotation information to the data. As Illumina organizes the output .idat files by barcodes, the function will automatically check the sub-folders in the names of barcodes for .idat files. The "sampleInfo" parameter can be either a barcode vector, e.g., "7310440039_R04C02" "7310440039_R05C02". Or a data.frame with required columns of 'Sentrix_Barcode' and 'Sentrix_Position'. If "sampleInfo" is a data.frame, it will be added as the pData of the output MethyLumiM object.
}

\value{
A MethyLumiM object
}

\author{
Pan Du, Tim Triche
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[methylumi]{lumIDAT}}, \code{\link{lumiMethyR}}, \code{\link{addAnnotationInfo}}
}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

