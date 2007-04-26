\name{targetID2nuID}
\alias{targetID2nuID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mapping Illumina TargetID into nuID}
\description{
}
\usage{
targetID2nuID(targetID, lib = "lumiHumanV1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{targetID}{ a vector of Illumina TargetID }
  \item{lib}{ a nuID indexed Illumina annotation library }
}
\details{
}
\value{
return a list with the same length as the input TargetID vector.
}

\author{Pan Du}
\seealso{See also \code{\link{nuID2targetID}} }
\examples{
if (require(lumiHumanV1)) {
    targetID2nuID('GI_21389350-S', lib='lumiHumanV1')
}
}
% R documentation directory.
\keyword{ methods }