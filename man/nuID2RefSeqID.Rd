\name{nuID2RefSeqID}
\alias{nuID2RefSeqID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Map nuID to RefSeq ID }
\description{
  Map nuID to RefSeq ID based on IDMapping libraries.
}
\usage{
nuID2RefSeqID(nuID = NULL, lib.mapping, filterTh = c(Strength1 = 95, Uniqueness = 95), returnAllInfo = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nuID}{ a vector of nuIDs. If it is NULL, all mappings will be returned. }
  \item{lib.mapping}{ the ID mapping library }
  \item{filterTh}{ the mapping quality filtering threshold used to filter the ID mapping. }
  \item{returnAllInfo}{ determine to return the detailed mapping information or just the matched RefSeq IDs }
}
\details{
  This function is based on the return of \code{\link{getNuIDMappingInfo}} function. It uses mapping quality information to filter out the bad mappings. The names of "filterTh" are basically the field names of "nuID\_MappingInfo" table, which include 'Strength1', 'Strength2', 'Uniqueness' and 'Total hits'. For the definition of these metrics, please refer to the IDMapping library or see the reference website.
}
\value{
	returns the matched RefSeq IDs or a data.frame with each row corresponding to an input nuID (when "returnAllInfo" is TRUE).	
}
\references{ 
 https://prod.bioinformatics.northwestern.edu/nuID/
}
\author{Warren Kibbe, Pan Du, Simon Lin}
\seealso{ See Also \code{\link{getNuIDMappingInfo}} }
\examples{
	## load example data
	data(example.lumi)
	if (require(lumiHumanIDMapping)) {
		nuIDs <- featureNames(example.lumi)
		mappingInfo <- nuID2RefSeqID(nuIDs, lib.mapping='lumiHumanIDMapping')
		head(mappingInfo)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{utilities}
