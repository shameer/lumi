\name{addAnnotationInfo}
\alias{addAnnotationInfo}
\title{Add probe color channel and basic annotation information based on the annotation library of Illumina methylation microarray }
\description{
Add probe color channel and basic annotation information based on the annotation library of Illumina methylation microarray
}
\usage{
addAnnotationInfo(methyLumiM, lib = NULL, annotationColumn=c('COLOR_CHANNEL', 'CHROMOSOME', 'POSITION'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{lib}{Annotation library of Illumina methylation microarray. For example, "IlluminaHumanMethylation27k.db" for 27k array, "IlluminaHumanMethylation450k.db" for 450k infinium array }
  \item{annotationColumn}{only include 'COLOR_CHANNEL', 'CHROMOSOME' and 'POSITION' information}
}
\value{
return the MethyLumiM object with COLOR_CHANNEL, CHROMOSOME and chromome POSITION information added to the featureData.
}

\author{Pan DU}

\seealso{
\code{\link{lumiMethyR}}
}
\examples{
	data(example.lumiMethy)
	head(pData(featureData(example.lumiMethy)))
	## removing color channel information
	# testData = example.lumiMethy
	# pData(featureData(testData))$COLOR_CHANNEL = NULL
	# testData = addAnnotationInfo(testData, lib="IlluminaHumanMethylation27k.db")
	## check whether the color channel information is added
	# head(pData(featureData(testData)))
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
