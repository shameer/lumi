\name{nuID2IlluminaID}
\alias{nuID2IlluminaID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matching nuIDs to Illumina IDs based on Illumina ID mapping library}
\description{
  Matching nuIDs to Illumina IDs based on Illumina ID mapping library
}
\usage{
nuID2IlluminaID(nuID, lib.mapping=NULL, species = c("Human", "Mouse", "Rat", "Unknown"), chipVersion = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nuID}{ a vector of nuIDs }
  \item{lib.mapping}{ the ID mapping library. If it is provided, the parameter "species" will be ignored. }
  \item{species}{ the species of the chip designed for. If users do not know it, it can be set as "Unknown".}
  \item{chipVersion}{ chipVersion information returned by function \code{\link{getChipInfo}} }
  \item{\dots}{ other parameters of \code{\link{getChipInfo}} }
}
\details{
When the parameter "chipVersion" is not provided, this function basically returned the "idMapping" item returned by function \code{\link{getChipInfo}}. 
}
\value{
The mapping information from nuID to Illumina ID. It will be a matrix with each column corresponding to one matched manifest file when parameter "returnAllMatches" is TRUE. In this case, the columns are sorted from the best match to worst.
}
\author{Pan Du}
\seealso{ See Also \code{\link{getChipInfo}}, \code{\link{IlluminaID2nuID}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{utilities}
