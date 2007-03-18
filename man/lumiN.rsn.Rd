\name{lumiN.rsn}
\alias{lumiN.rsn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Robust Spline Normalization between chips  }
\description{
  Robust spline normalization (monotonic curves) between chips
}
\usage{
lumiN.rsn(x.lumi, targetArray = NULL, excludeFold = 2, span = 0.03, ifPlot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes }
  \item{targetArray}{ a target chip which is the model for other chips for normalization }
  \item{excludeFold}{ exclude the genes with fold change larger than "excludeFold" during fitting the curve in normalization }
  \item{span}{ the span parameter used by \code{\link{monoSmu}} }
  \item{ifPlot}{ determine whether to plot intermediate results }
  \item{\dots} { other parameters used by \code{\link{monoSmu}} }
}
\details{

}
\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
If it is a LumiBatch object, it also includes the VST transform function and its parameters as attributes: "transformFun", "parameter". See \code{\link{inverseVST}} for details.
}
\references{ 
Lin, S.M., Du, P., Kibbe, W.A.,  "Model-based Variance-stabilizing Transformation for Illumina Mi-croarray Data", submitted
}
\author{Pan Du, Simon Lin}
\note{ 

}
\seealso{ \code{\link{lumiN}}, \code{\link{monoSmu}} }
\examples{

}
\keyword{ methods }
