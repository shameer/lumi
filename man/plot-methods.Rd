\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot.LumiQC}
\alias{plot,LumiQC-method}
\alias{plot,LumiQC,missing-method}

\title{Plot of a LumiQC object}
\description{
Creating quality control plots of a LumiQC object
}

\usage{
\S4method{plot}{LumiQC,missing}(x, what = c("density", "boxplot", "pair", "MAplot", "sampleRelation", "outlier", "cv"), 
parameterList = NULL, main, ...)
}

\arguments{
\item{x}{ a LumiQC object returned by \code{\link{lumiQ}} }
\item{what}{ one of the six kinds of QC plots }
\item{parameterList}{ the parameters for the corresponding QC plots }
\item{main}{ the title of the QC plot }
\item{\dots}{ additional parameters for the corresponding QC plots }
}

\details{
The parameter "what" of \code{plot} function controls the type of QC plots, which includes:
\item{\bold{density}:}{ the density plot of the chips, see \code{\link{hist-methods}}}
\item{\bold{boxplot}:}{ box plot of the chip intensities, see \code{\link{boxplot-methods}}}
\item{\bold{pair}:}{ the correlation among chips, plot as a hierarchical tree, see \code{\link{pairs-methods}}}
\item{\bold{MAplot}:}{ the MAplot between chips, see \code{\link{MAplot-methods}}}
\item{\bold{sampleRelation}:}{ plot the sample relations. See \code{\link{getSampleRelation}}}
\item{\bold{outlier}:}{ detect the outliers based on the sample distance to the center. See 
	\code{\link{detectOutlier}}}
\item{\bold{cv}:}{ the density plot of the coefficients of variance of the chips. See \code{\link{estimateLumiCV}}}
}

\seealso{ 
\code{\link{LumiQC-class}}, \code{\link{hist-methods}}, \code{\link{boxplot-methods}}, \code{\link{MAplot-methods}},
 	\code{\link{pairs-methods}}, \code{\link{getSampleRelation}}, \code{\link{estimateLumiCV}}, \code{\link{detectOutlier}} 
}

\examples{

## load example data
data(example.lumi)

## Quality control estimation
lumi.Q <- lumiQ(example.lumi)

## summary
summary(lumi.Q)

## plot the density
plot(lumi.Q, what='density')

## plot the pairwise sample correlation
plot(lumi.Q, what='pair')

## plot the pairwise MAplot
plot(lumi.Q, what='MAplot')

## sample relations
plot(lumi.Q, what='sampleRelation', parameterList=list(method='mds', col=c('100US', '95US:5P', '100US', '95US:5P')))

## detect outlier based on the distance to the mean profile
plot(lumi.Q, what='outlier')

## Density plot of coefficient of variance
plot(lumi.Q, what='cv')

}
\keyword{methods}
\keyword{hplot}
