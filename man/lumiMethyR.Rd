\name{lumiMethyR}
\alias{lumiMethyR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reading Illumina methylation microarray data
}
\description{
This function is a wrap of \code{\link[methylumi]{methylumiR}} function in methylumi package.
}
\usage{
lumiMethyR(filename, lib=NULL, controlData=NULL, qcfile=NULL, sampleDescriptions=NULL, sep = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ file name output by GenomeStudio}
  \item{lib}{Annotation library of Illumina methylation microarray}
  \item{controlData}{the controlData file name or a MethyLumiQC object to be added to the "controlData" slot of the MethyLumiM object}
  \item{qcfile} { parameter of \code{\link[methylumi]{methylumiR}} function}
  \item{sampleDescriptions} { parameter of \code{\link[methylumi]{methylumiR}} function}
  \item{sep} { parameter of \code{\link[methylumi]{methylumiR}} function}
}
\details{
This function is a wrap of \code{\link[methylumi]{methylumiR}} function in methylumi package. It will coerce the returned object as MethyLumiM class. The methylated and unmethylated probe intensity information is required for color-bias adjustment and normalization. If users have the Illumina IDAT files, we suggest use  \code{\link{importMethyIDAT}} function to import the data. The  \code{\link{importMethyIDAT}} function will automatically retrieve the required information and return a MethyLumiM object. 
}
\value{ return a MethyLumiM object }

\author{
Pan Du
}

\seealso{
  See Also \code{\link{importMethyIDAT}}, \code{\link[methylumi]{methylumiR}} and \code{\link{addControlData2methyLumiM}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
