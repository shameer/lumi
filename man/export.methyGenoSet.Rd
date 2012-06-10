\name{export.methyGenoSet}
\alias{export.methyGenoSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Export a MethyGenoSet object to be visualized in IGV or IGB 
}
\description{
Export a MethyGenoSet object to be visualized in IGV or IGB. Current version only supports "gct" or "bw" formats.
}
\usage{
export.methyGenoSet(methyGenoSet, file.format = c("gct", "bw"), exportValue = c("beta", "M"), savePrefix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyGenoSet}{
A MethyGenoSet object.
}
  \item{file.format}{
Export file format
}
  \item{exportValue}{
Export methylation values
}
  \item{savePrefix}{
The prefix used in the output filename.
}
}
\details{

}
\value{


}
\references{

}
\author{
Pan Du
}

\seealso{

}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
