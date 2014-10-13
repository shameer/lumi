\name{lumiT}
\alias{lumiT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Transfer the Illumina data to stabilize the variance }
\description{
  Transfer the Illumina data to stabilize the variance.
}
\usage{
lumiT(x.lumi, method = c("vst", 'log2', 'cubicRoot'), ifPlot = FALSE, simpleOutput = TRUE, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ LumiBatch object }
  \item{method}{ four methods are supported: "vst", "log2", "cubicRoot" }
  \item{ifPlot}{ determine whether to plot the intermediate results }
  \item{simpleOutput}{ determine whether to simplify the output LumiBatch object, which will set the se.exprs, detection and beadNum slots as NULL. }
  \item{verbose}{ a boolean to decide whether to print out some messages }
  \item{\dots}{ other parameters used by \code{\link{vst}} }
}
\details{
lumiT is an interface of difference variance stabilizing transformation. See \code{\link{vst}} for details of VST (Variance Stabilizing Transform) of Illumina data.

NOTE: This correction of STDERR as STDEV was moved to the lumiR function.
}

\value{
Return a LumiBatch object with transformed expression values. It also includes the VST transform function and its parameters as attributes: "transformFun", "parameter". See \code{\link{inverseVST}} for details.
	
}
\references{
 Lin, S.M., Du, P., Kibbe, W.A., (2008) 'Model-based  Variance-stabilizing Transformation for Illumina Microarray Data', Nucleic Acids Res. 36, e11
 }
\author{ Pan Du, Simon Lin }

\seealso{ \code{\link{vst}} }
\examples{

## load example data
data(example.lumi)

## Do default VST variance stabilizing transform
lumi.T <- lumiT(example.lumi, ifPlot=TRUE)


}
\keyword{ methods }
