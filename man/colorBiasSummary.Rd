\name{colorBiasSummary}
\alias{colorBiasSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A summary of colorBias information
}
\description{
A summary of colorBias information, which is a data.frame summarizing the intensities of individual samples
}
\usage{
colorBiasSummary(methyLumiM, logMode = TRUE, channel = c("both", "unmethy", "methy", "sum"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{logMode}{ Whether plot the intensities in log-scale }
  \item{channel}{estimate the intensity in different methods}
}
\details{
A summary of colorBias information. There are four options using "channel" parameter to plot the density plot. "both": estimate the density by pooling together methylated and unmethylated probe intensities. "unmethy" and "methy": plot either unmethylated or methylated probe density. "sum" plot the density of the sum of methylatled and unmethylated probe intensitys.}
\value{
A data.frame summarizing the intensities of individual samples
}

\author{
Pan DU
}

\examples{

}
\keyword{methods}
