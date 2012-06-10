\name{getChrInfo}
\alias{getChrInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get the chromosome location information of methylation probes
}
\description{
get the chromosome location information of methylation probes
}
\usage{
getChrInfo(methyData, lib = "IlluminaHumanMethylation450k.db", as.GRanges = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyData}{
a MethyLumiM object
}
  \item{lib}{
Methylation annotation library
}
  \item{as.GRanges}{
whether convert as GRanges or data.frame.
}
}
\details{

}
\value{
a GRanges or data.frame object 
}
\references{

}
\author{
Pan Du
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
