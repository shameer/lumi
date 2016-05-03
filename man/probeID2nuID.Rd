\name{probeID2nuID}
\alias{probeID2nuID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mapping Illumina ProbeID as nuID }
\description{
Mapping Illumina ProbeID as nuID. 
}
\usage{
probeID2nuID(probeID, lib.mapping = "lumiHumanIDMapping", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probeID}{ a vector of Illumina ProbeID }
  \item{lib.mapping}{ an Illumina ID mapping library }
  \item{\dots}{ other parameters of \code{\link{IlluminaID2nuID}} }
}
\details{
The function will call \code{\link{IlluminaID2nuID}} when ID mapping library were provided.
}
\value{
see function \code{\link{IlluminaID2nuID}}
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }
\seealso{ \code{\link{nuID2probeID}}, \code{\link{IlluminaID2nuID}} }
\examples{
if (require(lumiHumanIDMapping)) {
    probeID2nuID('0001240020', lib='lumiHumanIDMapping')
}

}
\keyword{ methods }
