\name{example.lumi}
\alias{example.lumi}
\docType{data}
\title{Example lumiBatch object includes example data}
\description{
  
}
\usage{data(example.lumi)}
\format{
Formal class 'LumiBatch' with 7 slots
  \cr \code{history}     		: 'data.frame':	1 obs. of  3 variables:
  \cr \code{assayData}			: length 4 <environment> 
  \cr \code{phenoData}			: Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  \cr \code{featureData}		: Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  \cr \code{experimentData}		: Formal class 'MIAME' [package "Biobase"] with 13 slots
  \cr \code{annotation}			: chr(0) 
  \cr \code{.\_\_classVersion\_\_}: Formal class 'Versions' [package "Biobase"] with 1 slots
}

\details{
The data is from (Chudin, et al., 2006). It used Sentrix Human-6 Expression BeadChip. Two samples "100US" and "95US.5P" (each has two technique replicates) were selected. In order to save space, 8000 genes were randomly selected. As a result, the example data includes 8000 genes, each has 4 measurements. 
}
\references{
  Chudin, E., Kruglyak, S., Baker, S.C., Oeser, S., Barker, D. and McDaniel, T.K. (2006) A model of technical variation of microar-ray signals, J Comput Biol, 13, 996-1003.
}
\examples{
## load the data
data(example.lumi)

## summary of the data
example.lumi

}
\keyword{datasets}
