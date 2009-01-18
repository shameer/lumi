\name{addNuID2lumi}
\alias{addNuID2lumi}
\alias{addNuId2lumi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add the nuID information to the LumiBatch object }
\description{ Replace the Illumina Id (Target ID or Probe Id) as nuID (nucleotide universal identifier) for indexing genes in the LumiBatch object }

\usage{
addNuID2lumi(x.lumi, annotationFile=NULL, sep = NULL, lib.mapping = NULL,  annotationColName = c(sequence = "Probe_Sequence", target = "ILMN_Gene", probe='Probe_Id'), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ a LumiBatch object }
  \item{annotationFile}{ a annotation file, which includes the Illumina ID (target or probe ids) and probe sequence information}
  \item{sep}{ the separation used in the annotation file. Automatically detect the separator if it is "," or "\t". }
  \item{lib.mapping}{ a Illumina ID mapping package, e.g, lumiHumanIDMapping }
  \item{annotationColName}{ the annotation column name in the annotation file used for the probe sequence and TargetID and ProbeID }
  \item{verbose}{ a boolean to decide whether to print out some messages }
}
\details{
Since the default Illumina IDs (TargetID (ILMN\_Gene ID) and ProbeId (Probe\_Id)) are not consistent between different arrays and batches, we invented a nuID, which is one-to-one matching with the probe sequence. This function is to replace the Illumina ID with the nuID. If the annotation library (the unzipped manifest file (.bgx)) is provided, the function will automatically check whether the Illumina ID is provided for the microarray data. We recommend output the data using ProbeID when using Illumina BeadStudio software, because the TargetID (ILMN\_Gene ID) are not unique. 
}
\value{
a LumiBatch object with Illumina ID replaced by nuID. 
}
\references{ Du, P., Kibbe, W.A., Lin, S.M., "nuID: A universal naming schema of oligonucleotides for Illumina, Affymetrix, and other microarrays", submitted.}
\author{ Pan Du }

\seealso{  \code{\link{IlluminaID2nuID}}, \code{\link{lumiR}} }
\examples{
## load example data
# data(example.lumi)

## specify the annotation file for the Illumina chip
# annotationFile <- 'Human_RefSeq-8.csv'
## Replace the Target ID with nuID
# lumi.nuID <- addNuID2lumi(example.lumi, annotationFile)

## An alternative way is to load the Annotation library and match the targetID (or Probe Id) with nuID
# lumi.nuId <- addNuID2lumi(example.lumi, lib.mapping='lumiHumanIDMapping')

}
\keyword{ methods }
