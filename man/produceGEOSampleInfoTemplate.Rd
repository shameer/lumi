\name{produceGEOSampleInfoTemplate}
\alias{produceGEOSampleInfoTemplate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce the template of GEO sample information }
\description{
  Produce the template of GEO sample information, which is used for function \code{\link{produceGEOSubmissionFile}}.
}
\usage{
produceGEOSampleInfoTemplate(lumiNormalized, lib.mapping = NULL, fileName = "GEOsampleInfo.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiNormalized}{ The normalized data (LumiBatch object) }
  \item{lib.mapping}{ The Illumina ID Mapping library, e.g., "lumiHumanIDMapping" }
  \item{fileName}{ The file name of Tab separated sample information file }
}
\details{
  This function just produces a template of sample information with some default fillings. Users need to fill in the detailed sample descriptions, especially the Sample\_title, Sample\_description and some protocols. No blank fields are allowed. Function \code{\link{produceGEOSubmissionFile}} will produce the file GEO submission file based on this sample information. The users should not use "\#" in the description as it is a reserved character.
}
\value{
	Save the result as a Tab separated text file or return a data.frame if the fileName is NULL.
}
\references{ http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html }
\author{ Pan Du }

\seealso{ \code{\link{produceGEOSubmissionFile}} }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
