\name{nuID2targetID}
\alias{nuID2targetID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mapping nuID into Illumina TargetID }
\description{
Mapping nuID into Illumina TargetID or GeneID.
}
\usage{
nuID2targetID(nuID, lib.mapping = "lumiHumanIDMapping.db", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nuID}{ a vector of nuID }
  \item{lib.mapping}{ an Illumina ID mapping library }
  \item{\dots}{ other parameters of \code{\link{nuID2IlluminaID}} }
}
\details{
The function will call \code{\link{nuID2IlluminaID}} when ID mapping library were provided.
}
\value{
see function \code{\link{nuID2IlluminaID}}
}
\references{
Du, P., Kibbe, W.A. and Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", Biology Direct 2007, 2:16 (31May2007).
}
\author{ Pan Du }
\seealso{ \code{\link{targetID2nuID}}, \code{\link{nuID2IlluminaID}} }
\examples{
if (require(lumiHumanIDMapping.db)) {
    nuID2targetID("B2J6WGhV.RevOJYff4", lib.mapping = "lumiHumanIDMapping.db")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
