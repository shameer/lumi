\name{IlluminaID2nuID}
\alias{IlluminaID2nuID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Matching Illumina IDs to nuID based on Illumina ID mapping library }
\description{
  Matching Illumina IDs to nuID based on Illumina ID mapping libraries.
}
\usage{
IlluminaID2nuID(IlluminaID, lib.mapping=NULL, species = c("Human", "Mouse", "Rat", "Unknown"), chipVersion = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IlluminaID}{ a vector of Illumina IDs }
  \item{lib.mapping}{ the ID mapping library. If it is provided, the parameter "species" will be ignored. }
  \item{species}{ the species of the chip designed for. If users do not know it, it can be set as "Unknown". }
  \item{chipVersion}{ chipVersion information returned by function \code{\link{getChipInfo}}  }
  \item{\dots}{ other parameters of \code{\link{getChipInfo}} }
}
\details{
When the parameter "chipVersion" is not provided, this function basically returned the "idMapping" item returned by function \code{\link{getChipInfo}}. 
}
\value{
The mapping information from Illumina ID to nuID. It will be a matrix with each column corresponding to one matched manifest file when parameter "returnAllMatches" is TRUE. In this case, the columns are sorted from the best match to worst. If IlluminaID is NULL and chipVersion is provided, it will return all mapping information of the chip. 

}
\author{Pan Du}
\seealso{ \code{\link{getChipInfo}}, \code{\link{nuID2IlluminaID}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{utilities}
