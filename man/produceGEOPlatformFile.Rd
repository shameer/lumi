\name{produceGEOPlatformFile}
\alias{produceGEOPlatformFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce GEO Platform Submission File in SOFT format }
\description{
	Produce GEO Sample Submission File in SOFT format based on the provided LumiBatch object and Illumina ID Mapping library
}
\usage{
produceGEOPlatformFile(x.lumi, lib.mapping = NULL, nuIDMode = FALSE, fileName = "GEOPlatformFile.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ The LumiBatch object keeping all probes }
  \item{lib.mapping}{ The Illumina ID Mapping library, e.g., "lumiHumanIDMapping.db" }
  \item{nuIDMode}{ Determine whether producing the platform indexed by nuID }
  \item{fileName}{ Filename of the GEO Platform File name }
}
\details{
 The function produces the GEO platform submission file based on the chip information kept in the Illumina ID Mapping library (specified by lib.mapping parameter). The determination of chip type will be automatically done by selecting the best matching of the probe IDs with individual chips.
}
\value{
Save the result as a text file in SOFT platform submission format.
}

\references{ http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html }
\author{ Pan Du }

\seealso{ \code{\link{produceGEOSubmissionFile}} }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
