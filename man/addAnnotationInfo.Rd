\name{addAnnotationInfo}
\alias{addAnnotationInfo}
\title{Add probe color channel and basic annotation information based on the annotation library of Illumina methylation microarray }
\description{
Add probe color channel and basic annotation information based on the annotation library of Illumina methylation microarray
}
\usage{
addAnnotationInfo(methyLumiM, lib = 'FDb.InfiniumMethylation.hg19', annotationColumn=c('COLOR_CHANNEL', 'CHROMOSOME', 'POSITION'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{lib}{Annotation library of Illumina methylation microarray. }
  \item{annotationColumn}{only include 'COLOR_CHANNEL', 'CHROMOSOME' and 'POSITION' information}
}
\value{
return the MethyLumiM object with COLOR_CHANNEL, CHROMOSOME and chromome POSITION information added to the featureData.
}
\details{
The "lib" parameter supports both FeatureDb based annotation libraries and old array-based annotation libraries.
'FDb.InfiniumMethylation.hg19' is the FeatureDb based annotation library, which includes both 450k and 27k data. 
"IlluminaHumanMethylation27k.db" (for 27k array) and "IlluminaHumanMethylation450k.db" (450k infinium array) are old types of annotation libraries.
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
