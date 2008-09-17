\name{produceGEOSubmissionFile}
\alias{produceGEOSubmissionFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce GEO Sample Submission File in SOFT format }
\description{
  Produce GEO Sample Submission File in SOFT format based on the provided LumiBatch object and sample information
}
\usage{
produceGEOSubmissionFile(lumiNormalized, lumiRaw, lib.mapping, sampleInfo = NULL, fileName = "GEOSubmissionFile.txt", supplementaryRdata = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lumiNormalized}{ The normalized data (LumiBatch object) }
  \item{lumiRaw}{  The raw data (LumiBatch object), e.g., returned by \code{\link{lumiR}} }
  \item{lib.mapping}{ The Illumina ID Mapping library, e.g., "lumiHumanIDMapping.db" }
  \item{sampleInfo}{ The sample information filename or data.frame, which is returned by \code{\link{produceGEOSampleInfoTemplate}} }
  \item{fileName}{ The file name of GEO Submission file  }
  \item{supplementaryRdata}{ determine whether produce the Rdata supplement data, which include both lumiNormalized and lumiRaw R objects. }
}
\details{
The function produces the GEO sample submission file including both normalized and raw data information. The sample information should be provided by the user as a data.frame or Tab separated text file. The template of the sample information file can be produced by function \code{\link{produceGEOSampleInfoTemplate}}. Users need to fill in the detailed sample descriptions in the template, especially the Sample\_title, Sample\_description and some protocols. 
}
\value{
Save the result as a text file in SOFT sample submission format. The supplementary Rdata will be saved in a file 'supplementaryData.Rdata'.
}
\references{ http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html }
\author{ Pan Du }

\seealso{ \code{\link{produceGEOSampleInfoTemplate}}, \code{\link{produceGEOPlatformFile}} }
\examples{
## Not run
## 'GEOsampleInfo.txt' is file with filled sample information
# produceGEOSubmissionFile(lumiNormalized, lumiRaw, lib='lumiHumanIDMapping.db', sampleInfo='GEOsampleInfo.txt')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
