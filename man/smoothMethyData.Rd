\name{smoothMethyData}
\alias{smoothMethyData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth the methylation data
}
\description{
Smooth the methylation data by a sliding window with fixed width in bp unit
}
\usage{
smoothMethyData(methyData, winSize = 250, lib = "IlluminaHumanMethylation450k.db", asDataFrame = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyData}{
A MethyLumiM or MethyGenoSet object
}
  \item{winSize}{
Sliding window size in bp unit
}
  \item{lib}{
Methylation annotation library
}
  \item{asDataFrame}{
Whether return the smoothed value as a data.frame.
}
  \item{\dots}{
other parameters
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object with the methylation values smoothed
}
\references{
}
\author{
Pan Du
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
