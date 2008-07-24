\name{targetID2nuID}
\alias{targetID2nuID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mapping Illumina TargetID (GeneID) into nuID}
\description{
Mapping Illumina TargetID (GeneID) into nuID.
}
\usage{
targetID2nuID(targetID, lib.mapping = "lumiHumanIDMapping.db", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targetID}{ a vector of Illumina TargetID (GeneID) }
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
\author{Pan Du}
\seealso{ \code{\link{nuID2targetID}}, \code{\link{IlluminaID2nuID}} }
\examples{
if (require(lumiHumanIDMapping.db)) {
    targetID2nuID('GI_21389350-S', lib='lumiHumanIDMapping.db')
}
}
% R documentation directory.
\keyword{ methods }
