\name{ssn}
\alias{ssn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simple Scaling Normalization }
\description{
  This function basically adjusts the samples to the same background level and scales based on the background level and foreground level.
}
\usage{
ssn(x.lumi, targetArray = NULL, scaling = TRUE, bgMethod=c('density', 'mean', 'median'), fgMethod=c('mean', 'density', 'median'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes }
  \item{targetArray}{ A target chip is the model for other chips to normalize. It can be a column index, a vector or a LumiBatch object with one sample. }
  \item{scaling}{ determine whether do scaling or just background shift }
  \item{bgMethod}{ optional methods of determining the background level }
  \item{fgMethod}{ optional methods of determining the foreground level }
  \item{\dots}{ other parameters used by \code{\link[stats]{density}} function }
}
\details{
	This function basically adjusts the samples to the same background level and scales based on the background level and foreground level. The adjustment is based on the raw scale data.
}
\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
}
\references{ }
\author{ Pan Du, Simon Lin }
\seealso{ \code{\link{lumiN}} }
\examples{

}
\keyword{ methods }
