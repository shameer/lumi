\name{addColorChannelInfo}
\alias{addColorChannelInfo}
\title{Add probe color channel based on the annotation library of Illumina methylation microarray }
\description{
Add probe color channel based on the annotation library of Illumina methylation microarray
}
\usage{
addColorChannelInfo(methyLumiM, lib = "IlluminaHumanMethylation27k.db")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{a MethyLumiM object includes Illumina Infinium methylation data}
  \item{lib}{Annotation library of Illumina methylation microarray}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return the MethyLumiM object with "COLOR_CHANNEL" column added to the featureData.
}

\author{Pan DU}

\seealso{
\code{\link{lumiMethyR}}
}
\examples{
	data(example.lumiMethy)
	head(pData(featureData(example.lumiMethy)))
	# removing color channel information
	testData = example.lumiMethy
	pData(featureData(testData))$COLOR_CHANNEL = NULL
	testData = addColorChannelInfo(testData, lib="IlluminaHumanMethylation27k.db")
	# check whether the color channel information is added
	head(pData(featureData(testData)))
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
