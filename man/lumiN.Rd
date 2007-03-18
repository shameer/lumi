\name{lumiN}
\alias{lumiN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Between chip normalization of a LumiBatch object }
\description{
  A main function of between chip normalization of a LumiBatch object. Currently, four methods ("RSN", "loess", "quantile", "VSN") are supported.
}
\usage{
lumiN(x.lumi, method = c("RSN", "loess", "quantile", "VSN"), targetArray=NULL, ifPlot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes}
  	\item{method}{ four different between chips normalization methods ("RSN", "loess", "quantile", "VSN") are supported }
	\item{targetArray}{ the target array to which others are aligned. If it is null, the one having the most similar profile with the mean profile is selected. }
	\item{ifPlot}{ determine whether to plot intermediate results }
  	\item{\dots}{ other parameters used by \code{\link{lumiN.rsn}} }
}
\details{

}

\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
If it is a LumiBatch object, it also includes the VST transform function and its parameters as attributes: "transformFun", "parameter". See \code{\link{inverseVST}} for details.
}
\references{  }
\author{Pan Du, Simon Lin}
\note{ 

}
\seealso{  \code{\link{lumiN.rsn}} }
\examples{
## load example data
data(example.lumi)

## Do lumi transform
lumi.T <- lumiT(example.lumi)

## Do lumi between chip normaliazation
lumi.N <- lumiN(lumi.T, ifPlot=TRUE)

}
\keyword{ methods }
