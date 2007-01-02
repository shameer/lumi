\name{addNuId2lumi}
\alias{addNuId2lumi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add the nuID information to the LumiBatch object }
\description{ Replace the Illumina Id (Target ID) as nuID (nucleotide universal identifier) for indexing genes in the LumiBatch object }

\usage{
addNuId2lumi(x.lumi, annotationFile=NULL, sep = NULL, lib = NULL,  annotationColName = c(sequence = "Probe_Sequence", target = "Target"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{annotationFile}{ a annotation file, which includes the TargetID and probe sequence information}
  \item{sep}{ the separation used in the annotation file. Automatically detect the separator if it is "," or "\t". }
  \item{lib}{ a lumi annotation package }
  \item{annotationColName}{ the annotation column name used for the probe sequence and TargetID  }
}
\details{
Since the default TargetID may not be consistent between different arrays and batches, we invented a nuID, which is one-to-one matching with the probe sequence. This function is to replace the TargetID with the nuID.
}
\value{
a LumiBatch object with TargetID replaced by nuID. 
}
\references{ Du, P., Kibbe, W.A., Lin, S.M., "nuID: A Novel Identifier for Oligos, Ideal for Oligonucleotide-based  Microarrays", submitted.}
\author{ Pan Du }

\seealso{  \code{\link{lumiR}} }
\examples{
## load example data
# data(example.lumi)

## specify the annotation file for the Illumina chip
# annotationFile <- 'Human_RefSeq-8.csv'
## Replace the Target ID with nuID
# lumi.nuId <- addNuId2lumi(example.lumi, annotationFile)

## An alternative way is to load the Annotation library and match the targetID with nuID
## This is not recommended because some targetIDs of Illumina microarray are recycled for different chips 
##		and could cause mismatch between the library and the annotation file. Also the lumi annotation library 
##      may be outdated and does not include your annotation information.
# lumi.nuId <- addNuId2lumi(example.lumi, lib='Human.lumi')

}
\keyword{ methods }
