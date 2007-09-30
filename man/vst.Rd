\name{vst}
\alias{vst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Variance Stabilizing Transformation }
\description{
  Stabilizing the expression variance based on the bead level expression variance and mean relations
}
\usage{
vst(u, std, nSupport = min(length(u), 500), backgroundIndex=NULL, fitMethod = c('linear', 'quadratic'), ifPlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ mean expression of the beads with same sequence }
  \item{std}{ expression standard deviation of the beads with same sequence }
  \item{nSupport}{ the number of down-sampling to speed processing }
  \item{backgroundIndex}{ the index of background probes, which are determined based on the detection p-values.}
  \item{fitMethod}{ methods of fitting the relations between expression variance and mean relations }
  \item{ifPlot}{ plot intermediate results or not }
}
\details{
 The variance-stabilizing transformation (VST) takes the advantage of larger number of technical replicates available on the Illumina microarray. It models the mean-variance relationship of the within-array technical replicates at the bead level of Illumina microarray. An arcsinh transform is then applied to stabilize the variance. See reference for more details.

For the methods of fitting the relations between expression variance and mean relations, the 'linear' method is more robust and provides detailed parameters for inverseVST. 
}
\value{
Return the transformed (variance stabilized) expression values.
}
\references{
 Lin, S.M., Du, P., Kibbe, W.A.,  "Model-based Variance-stabilizing Transformation for Illumina Mi-croarray Data", submitted
 }
\author{ Pan Du, Simon Lin }

\seealso{ \code{\link{lumiT}}, \code{\link{inverseVST}} }
\examples{
## load example data
data(example.lumi)

## get the gene expression mean for one chip
u <- exprs(example.lumi)[,1]
## get the gene standard deviation for one chip
std <- se.exprs(example.lumi)[,1]

## do variance stabilizing transform
transformedU <- vst(u, std)

## do variance stabilizing transform with plotting intermediate result 
transformedU <- vst(u, std, ifPlot=TRUE)


}
\keyword{ methods }
