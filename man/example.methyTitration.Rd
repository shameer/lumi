\name{example.methyTitration}
\alias{example.methyTitration}
\docType{data}
\title{Example Illumina Infinium Methylation titration data in MethyLumiM class}
\description{
An example Illumina Infinium Methylation27k dataset, which includes a titration dataset
}
\usage{data(example.methyTitration)}
\details{
The example dataset is a titration dataset. It includes 8 samples: "A1", "A2", "B1", "B2", "C1", "C2", "D" and "E". They are mixtures of Sample A (a B-lymphocyte sample) and Sample B (is (a colon cancer sample from a female donor) at five different titration ratios: 100:0 (A), 90:10 (C), 75:25 (D), 50:50 (E) and 0:100 (B). 

To save storage space, we randomly subset 10000 CpG sites among about 27000 measured CpG sites.

}

\examples{
data(example.methyTitration)
sampleNames(example.methyTitration)
}
\keyword{datasets}
