\name{boxplotColorBias}
\alias{boxplotColorBias}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the Illumina Infinium methylation color bias in terms of boxplot
}
\description{
Plot the Illumina Infinium methylation color bias in terms of boxplot. boxplot of red and green color channel will be plotted side by side
}
\usage{
boxplotColorBias(methyLumiM, logMode = TRUE, channel = c("both", "unmethy", "methy", "sum"), grid = TRUE, main = NULL, mar = NULL, verbose = F, subset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{logMode}{ whether plot the intensities in log-scale }
  \item{channel}{ estimate the intensity in different methods }
  \item{grid}{ whether to add grid on the plot }
  \item{main}{ title of the plot }
  \item{mar}{ margin of the plot }
  \item{verbose}{ whether print verbose information during plot }
  \item{subset}{plot subset of randomly selected rows. All data will be plotted if it is NULL. }
  \item{\dots}{ other parameters of \code{\link[graphics]{boxplot}} }
}
\details{
Plot the Illumina Infinium methylation color bias in terms of boxplot. boxplot of red and green color channel will be plotted side by side
}
\value{
Invisibly return TRUE if plot successfully.
}
\author{
Pan DU
}

\seealso{
 See Also as \code{\link[graphics]{boxplot}} and \code{\link{plotColorBias1D}} 
}
\examples{
data(example.lumiMethy)
boxplotColorBias(example.lumiMethy)
}

\keyword{methods}
\keyword{hplot}
