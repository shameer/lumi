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
The data is from (Barnes, 2005). It used Sentrix HumanRef-8 Expression BeadChip. Two samples "100US" and "95US:5P" (each has two technique replicates) were selected. In order to save space, 8000 genes were randomly selected. As a result, the example data includes 8000 genes, each has 4 measurements. The full data set was included in the Bioconductor Experiment data package lumiBarnes.
}
\references{
Barnes, M., Freudenberg, J., Thompson, S., Aronow, B. and Pavlidis, P. (2005) Ex-perimental comparison and cross-validation of the Affymetrix and Illumina gene expression analysis platforms, Nucleic Acids Res, 33, 5914-5923.
}
\examples{
## load the data
data(example.lumi)

## summary of the data
example.lumi

}
\keyword{datasets}
