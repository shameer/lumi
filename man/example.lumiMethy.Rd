\name{example.lumiMethy}
\alias{example.lumiMethy}
\docType{data}
\title{Example Illumina Infinium Methylation data in MethyLumiM class}
\description{
An example Illumina Infinium Methylation27k dataset, which includes a control and treatment dataset with both technique and biological replicates
}
\usage{data(example.lumiMethy)}
\details{
The example dataset includes four control and four treatment samples together with their technique replicates. The original samples and technique replicates were measured in two batches. Here are the names of sixteen samples: Treat1, Treat2, Treat3, Treat4, Ctrl1, Ctrl2, Ctrl3, Ctrl4, Treat1.rep, Treat2.rep, Treat3.rep, Treat4.rep, Ctrl1.rep, Ctrl2.rep, Ctrl3.rep, Ctrl4.rep.

To save storage space, we randomly subset 5000 CpG sites among about 27000 measured CpG sites.

}

\examples{
data(example.lumiMethy)
sampleNames(example.lumiMethy)
}
\keyword{datasets}
