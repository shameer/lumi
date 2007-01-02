\name{getSampleRelation}
\alias{getSampleRelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate the sample relations }
\description{
  Estimate the sample relations based on MDS or hierarchical clustering
}
\usage{
getSampleRelation(x, selProbe = NULL, Th = 0.1, standardize = TRUE, method = c("mds", "cluster"), dimension = c(1, 2), col = NULL, ifPlot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an exprSet object or a matrix with each column as a chip}
  \item{selProbe}{ the selected probes used to determine the sample relations. If not provide, all the probes will be used. }
  \item{Th}{ the threshold of the coefficient of variance of probes used to select probes to estimate sample relations }
  \item{standardize}{ standardize the expression profiles or not }
  \item{method}{ "MDS" or "hierarchical clustering" }
  \item{dimension}{ the principle components to visualize the MDS plot }
  \item{col}{ the color for each sample during plot. Only support the "mds" method }
  \item{ifPlot}{ plot the results or not }
}
\details{
  Estimate the sample relations based on selected probes (based on large coefficient of variance (mean / standard variance)). Two methods can be used: MDS (Multi-Dimensional Scaling) or hierarchical clustering methods.
}
\value{
Plot the results or return the distance matrix
}

\author{ Pan Du }

\seealso{  \code{\link{lumiQ}}, \code{\link{LumiQC}}, , \code{\link{plot.LumiQC}} }
\examples{

## load example data
data(example.lumi)

## estimate the sample distance and plot the relations with MDS
## the color of sample is automatically set based on the sample type
temp <- getSampleRelation(example.lumi, col=c('100US', '95US:5P', '100US', '95US:5P'), ifPlot=TRUE)

## estimate the sample distance and plot the relations with hierarchical clustering
temp <- getSampleRelation(example.lumi, method='cluster', ifPlot=TRUE)

}
\keyword{ methods }
