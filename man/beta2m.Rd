\name{beta2m}
\alias{beta2m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert methylation Beta-value to M-value
}
\description{
Convert methylation Beta-value to M-value through a logistic transformation
}
\usage{
beta2m(beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{
	a matrix or vector of methylation Beta-value
}
}
\details{
Convert methylation Beta-value to M-value through a logistic transformation
}
\value{
return methylation M-value with the same size of input Beta-value
}
\references{
Du, P., Zhang, X, Huang, C.C., Jafari, N., Kibbe, W.A., Hou, L., and Lin, S.M.,  (2010) 'Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis', (under review)

}
\author{
Pan Du
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See Also as \code{\link{m2beta}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}

