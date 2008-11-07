\name{getNuIDMappingInfo}
\alias{getNuIDMappingInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ get the mapping information from nuID to RefSeq ID }
\description{
  Get the mapping information (including mapping quality information) of nuIDs to the most recent RefSeq release. These information was kept in the IDMapping libraries.
}
\usage{
getNuIDMappingInfo(nuID = NULL, lib.mapping)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nuID}{ a vector of nuIDs. If it is NULL, all mappings will be returned. }
  \item{lib.mapping}{ the ID mapping library }
}
\details{
 The function basically return the nuID mapping information kept in the "nuID\_MappingInfo" table of IDMapping libraries (lumiHumanIDMapping, lumiMouseIDMapping, lumiRatIDMapping). For more details of nuID mapping, please refer to the help of corresponding IDMapping library.
}
\value{
	It returns a data.frame with each row corresponding to an input nuID.
}
\author{Warren Kibbe, Pan Du, Simon Lin}
\examples{
	## load example data
	data(example.lumi)
	if (require(lumiHumanIDMapping)) {
		nuIDs <- featureNames(example.lumi)
		mappingInfo <- getNuIDMappingInfo(nuIDs, lib.mapping='lumiHumanIDMapping')
		head(mappingInfo)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{utilities}
