\name{produceMethylationGEOSubmissionFile}
\alias{produceMethylationGEOSubmissionFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce GEO Sample Submission File of Illumina methylation microarray data in SOFT format}
\description{
Produce GEO Sample Submission File in the SOFT format based on the provided MethyLumiM object and sample information
}
\usage{
produceMethylationGEOSubmissionFile(methyLumiM, methyLumiM.raw = NULL, lib.mapping = NULL, idType = "Probe", sampleInfo = NULL, fileName = "GEOSubmissionFile.txt", supplementaryRdata = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{The normalized data in MethyLumiM class}
  \item{methyLumiM.raw}{The raw data in MethyLumiM class}
  \item{lib.mapping}{Currently not used for Illumina methylation data}
  \item{idType}{Currently no other options for Illumina methylation data}
  \item{sampleInfo}{The sample information filename or data.frame, which is returned by \code{\link{produceGEOSampleInfoTemplate}}}
  \item{fileName}{The file name of GEO Submission file}
  \item{supplementaryRdata}{determine whether produce the Rdata supplement data, which include both methyLumiM and methyLumiM.raw R objects.}
  \item{\dots}{other parameters used by function \code{\link{nuID2IlluminaID}}, but not implemented for methylation data}
}
\details{
The function produces the GEO sample submission file including both normalized and raw data information  in the SOFT format. The sample information should be provided by the user as a data.frame or Tab separated text file following the format of the template, which can be produced by function \code{\link{produceGEOSampleInfoTemplate}}. Users need to fill in the detailed sample descriptions in the template, especially the Sample\_title, Sample\_description and some protocols. Users are also required to fill in the "Sample\_platform\_id" by checking information of the GEO Illumina platform.

When the parameter "supplementaryRdata" is TRUE, the R objects, methyLumiM, methyLumiM.raw and sampleInfo, will be saved in a file named 'supplementaryData.Rdata'.
}
\value{
Save the result as a text file in SOFT sample submission format. The supplementary Rdata will be saved in a file 'supplementaryData.Rdata'.
}
\references{ http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html }
\author{ Pan Du }

\seealso{ \code{\link{produceGEOSampleInfoTemplate}}, \code{\link{produceGEOPlatformFile}} }
\examples{
## Not run
## Produce the sample information template
# produceGEOSampleInfoTemplate(methyLumiM, fileName = "GEOsampleInfo.txt")
## After editing the 'GEOsampleInfo.txt' by filling in sample information
# produceMethylationGEOSubmissionFile(methyLumiM, methyLumiM.raw, sampleInfo='GEOsampleInfo.txt')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
