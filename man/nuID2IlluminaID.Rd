\name{nuID2IlluminaID}
\alias{nuID2IlluminaID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matching nuIDs to Illumina IDs based on Illumina ID mapping library}
\description{
  Matching nuIDs to Illumina IDs based on Illumina ID mapping library
}
\usage{
nuID2IlluminaID(nuID, lib.mapping=NULL, species = c("Human", "Mouse", "Rat", "Unknown"), idType=c('All', 'Probe', 'Gene', 'Accession', 'Search_key', 'Symbol'), chipVersion = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nuID}{ a vector of nuIDs }
  \item{lib.mapping}{ the ID mapping library. If it is provided, the parameter "species" will be ignored. }
  \item{species}{ the species of the chip designed for. If users do not know it, it can be set as "Unknown".}
  \item{idType}{ the Illumina ID type}
  \item{chipVersion}{ chipVersion information returned by function \code{\link{getChipInfo}} }
  \item{\dots}{ other parameters of \code{\link{getChipInfo}} }
}
\details{
The parameter "idType" represents different types of Illumina IDs. It returns the entire table when idType = "All". When idType = 'Probe', it returns "ProbeId" or "Probe\_Id". When idType = 'Gene', it returns "Target" or "ILMN\_Gene" IDs.
 
This function basically returned the "idMapping" item returned by function \code{\link{getChipInfo}}.
If nuID is NULL and chipVersion is provided, it will return all mapping information of the chip. 
}
\value{
The mapping information from nuID to Illumina ID. It will be a matrix with each column corresponding to one matched manifest file when parameter "returnAllMatches" is TRUE. In this case, the columns are sorted from the best match to worst.
}
\author{Pan Du}
\seealso{ \code{\link{getChipInfo}}, \code{\link{IlluminaID2nuID}} }
\examples{
	## load example data
	data(example.lumi)
	nuIDs <- featureNames(example.lumi)
	if (require(lumiHumanIDMapping.db)) {
		illuminaID <- nuID2IlluminaID(nuIDs[1:5], lib='lumiHumanIDMapping.db')
		illuminaID
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{utilities}
