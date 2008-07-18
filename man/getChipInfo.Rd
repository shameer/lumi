\name{getChipInfo}
\alias{getChipInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get Illumina Chip Information based on probe identifiers}
\description{
	Retrieve the matched Illumina chip information by searching the provided probe identifiers through the Illumina identifiers in all manifest files.
}
\usage{
getChipInfo(x, lib.mapping = NULL, species = c("Human", "Mouse", "Rat", "Unknown"), idMapping = FALSE, returnAllMatches = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of probe identifiers, ExpressionSet object or a matrix with probe identifiers as row names }
  \item{lib.mapping}{ the ID mapping library. If it is provided, the parameter "species" will be ignored. }
  \item{species}{ species of the chip designed for. If users do not know it, it can be set as "Unknown". }
  \item{idMapping}{ determine whether return the idMapping information (between Illumina ID and nuID) }
  \item{returnAllMatches}{ determine whether return all matches or just the best match }
  \item{verbose}{ determine whether print some warning information }
}
\details{
  The function searches the provided probe Identifiers (Illumina IDs or nuIDs) through all the manifest file ID information kept in the IDMapping libraries (lumiHumanIDMapping.db, lumiMouseIDMapping.db, lumiRatIDMapping.db). The Illumina IDs kept in the library include "Search\_key" ("Search\_Key"), "Target" ("ILMN\_Gene"), "Accession", "Symbol", "ProbeId" ("Probe\_Id"). To determine the best match, the function calculate the number of matched probes. The higher "matchProbeNumber" is claimed as better. When the "matchProbeNumber" is the same, the manifest file will fewer probes is claimed as better.
}
\value{
  The function returns a list with following items:
  \item{chipVersion}{the file name of the manifest file for the corresponding version and release}
  \item{species}{the species of the chip designed for}
  \item{IDType}{the type of probe identifier}
  \item{chipProbeNumber}{the number of probes in the manifest file}
  \item{matchProbeNumber}{the number of input probes matching the manifest file}
  \item{idMapping}{id mapping information between Illumina ID and nuID}
  When parameter "returnAllMatches" is TRUE, the items of "chipVersion", "IDType", "chipProbeNumber", "matchProbeNumber" will be a vector corresponding to the matched manifest files, whose "matchProbeNumber" is larger than zero, and the "idMapping" will be a matrix with each column corresponding to one matched manifest file. All of the items are sorted from the best match to worst (The higher "matchProbeNumber" is claimed as better. When the "matchProbeNumber" is the same, the manifest file with fewer probes is claimed as better.). 
}
\references{  }
\author{Pan Du}
\seealso{ See Also \code{\link{help}} }
\examples{
	## load example data
	data(example.lumi)
	if (require(lumiHumanIDMapping.db)) {
		chipInfo <- getChipInfo(example.lumi, species='Human')
		chipInfo
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{utilities}
