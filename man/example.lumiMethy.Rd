\name{example.lumiMethy}
\alias{example.lumiMethy}
\docType{data}
\title{Example Illumina Infinium Methylation data in MethyLumiM class}
\description{
An example Illumina Infinium Methylation27k dataset, which includes two sub datasets: a titration dataset, and a control and treatment dataset
}
\usage{data(example.lumiMethy)}
\details{
The example dataset includes 12 samples composed by two sub datasets: a titration dataset, and a control and treatment dataset. The titration dataset includes 8 samples: "A1", "A2", "B1", "B2", "C1", "C2", "D" and "E". They are mixtures of Sample A (a B-lymphocyte sample) and Sample B (is (a colon cancer sample from a female donor) at five different titration ratios: 100:0 (A), 90:10 (C), 75:25 (D), 50:50 (E) and 0:100 (B). 

Another sub dataset is composed of two control samples and two treatment samples: "Ctrl1", "Ctrl2", "Treat1" and "Treat2". The two technique replicates were measured in different batches. We can see clear batch effects in colors bias. 

To save storage space, we randomly subset 10000 CpG sites among about 27000 measured CpG sites.

}

\references{

}
\examples{
data(example.lumiMethy)
sampleNames(example.lumiMethy)
}
\keyword{datasets}
