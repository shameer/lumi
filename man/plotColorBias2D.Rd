\name{plotColorBias2D}
\alias{plotColorBias2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the color bias of Illumina Infinium Methylation data in two dimensions
}
\description{
Plot the color bias (red and green channel) of Illumina Infinium Methylation data of one selected sample in two dimensions (methylated and unmethylated probe intensities)
}
\usage{
plotColorBias2D(methyLumiM, selSample = 1, combineMode = F, layoutRatioWidth = c(0.75, 0.25), layoutRatioHeight = c(0.25, 0.75), margins = c(5, 5, 2, 2), cex = 1.25, logMode = TRUE, main = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{selSample}{ The index of sample name of the selected sample to plot color bias }
  \item{combineMode}{ Whether combine two color channels together and plot as one color}
  \item{layoutRatioWidth}{ the plot figure ratio between scatter plot and density plot}
  \item{layoutRatioHeight}{the plot figure ratio between scatter plot and density plot}
  \item{margins}{margin of the plot}
  \item{cex}{A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. See \code{\link[graphics]{par}} }
  \item{logMode}{ Whether plot the intensities in log-scale}
  \item{main}{ the title of the plot }
}
\details{
The function basically plots the probe intensities in 2-dimension (methylated vs unmethylated), and colors the dots in Red and Green based on their color channel information. The related density plot will also be plotted at the right and top of the scatter plot.
}
\value{
Invisibly return TRUE if plot successfully.
}
\references{

}
\author{
Pan DU
}

\seealso{
 See Also as \code{\link{plotColorBias1D}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}
