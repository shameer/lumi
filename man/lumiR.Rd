\name{lumiR}
\alias{lumiR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read in Illumina expression data }
\description{
  Read in Illumina expression data. We assume the data was saved in a comma or tab separated text file.  
}
\usage{
lumiR(fileName, sep = NULL, detectionTh = 0.01, na.rm = TRUE, convertNuID = TRUE, lib.mapping = NULL, dec = '.', parseColumnName = FALSE, checkDupId = TRUE,
QC = TRUE, columnNameGrepPattern = list(exprs='AVG_SIGNAL', se.exprs='BEAD_STD', detection='DETECTION', beadNum='Avg_NBEADS'),
inputAnnotation=TRUE, annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION'), verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{ fileName of the data file }
  \item{sep}{ the separation character used in the text file.  }
  \item{detectionTh}{ the p-value threshold of determining detectability of the expression. See more details in \code{\link{lumiQ}} }
  \item{na.rm}{ determine whether to remove NA }
  \item{convertNuID}{ determine whether convert the probe identifier as nuID }
  \item{lib.mapping}{ a Illumina ID mapping package, e.g, lumiHumanIDMapping, used by \code{\link{addNuID2lumi}} }
  \item{dec}{ the character used in the file for decimal points.}
  \item{parseColumnName}{ determine whether to parse the column names and retrieve the sample information (Assume the sample information is separated by "\_".) }
  \item{checkDupId}{ determine whether to check duplicated TargetIDs or ProbeIds. The duplicated ones will be averaged.}
  \item{QC}{ determine whether to do quality control assessment after read in the data.}
  \item{columnNameGrepPattern}{ the string grep patterns used to determine the slot corresponding columns.}
  \item{inputAnnotation}{ determine whether input the annotation information outputted by BeadStudio if exists.}
  \item{annotationColumn}{ the column names of the annotation information outputted by BeadStudio}
  \item{verbose}{ a boolean to decide whether to print out some messages }
  \item{\dots}{ other parameters used by \code{\link[utils]{read.table}} function }
}
\details{
The function can automatically determine the separation character if it is Tab or comma. Otherwise, the user should specify the separator manually. If the annotation library is provided, the Illumina Id will be replaced with nuID, which is used as the index Id for the lumi annotation packages. If the annotation library is not provided, it will try to directly convert the probe sequence (if provided in the BeadStudio output file) as nuIDs.

The parameter "columnNameGrepPattern" is designed for some advanced users. It defines the string grep patterns used to determine the slot corresponding columns. For example, for the "exprs" slot in LumiBatch object, it is composed of the columns whose name includes "AVG\_SIGNAL". In some cases, the user may not want to read the "detection" and "beadNum" related columns to save memory. The user can set the "detection" and "beadNum" as NA in "columnNameGrepPattern". If the 'se.exprs' is set as NA or the corresponding columns are not available, then lumiR will create a ExpressionSet object instead of LumiBatch object.

The parameter "parseColumnName" is designed to parse the column names and retrieve the sample information. We assume the sample information is separated by "\_" and the last element after "\_" is the sample label (sample names of the LumiBatch object).  If the parsed sample labels are not unique, then the entire string will be used as the sample label. For example: "1881436055\_A\_STA 27aR" is included in one of the column names of BeadStudio output file. Here, the program will first treat "STA 27aR" as the sample label. If it is not unique across the samples, "1881436055\_A\_STA 27aR" will be the sample label. If it is still not unique, the program will report warning messages. All the parsed information is kept in the phenoData slot. By default, "parseColumnName" is FALSE. We suggest the users use it only when they know what they are doing.

Current version of lumiR can adaptively read the output of BeadStudio Verson 1 and 3. The format Version 3 made quite a few changes comparing with previous versions. One change is the detection value. It was called detectable when the detection value is close to one for Version 1 format. However, the detection value became a p-value in the Version 3. As a result, the detectionTh is automatically changed based on the version. The detectionTh 0.01 for the Version 3 will be changed as the detectionTh 0.99 for Version 1. Another big change is that Version 3 separately output the control probe (gene) information and a "Samples Table". As a result, the controlData slot in LumiBatch class was added to keep the control probe (gene) information, and a QC slot to keep the quality control information, including the "Sample Table" output by BeadStudio version 3.

The recent version of BeadStudio can also output the annotation information together with the expression data. In the users also want to input the annotation information, they can set the parameter "inputAnnotation" as TRUE. At the same time, they can also specify which columns to be inputted by setting parameter "annotationColumn". The BeadStudio annotation columns include: SPECIES, TRANSCRIPT, ILMN\_GENE, UNIGENE\_ID, GI, ACCESSION, SYMBOL, PROBE\_ID, ARRAY\_ADDRESS\_ID, PROBE\_TYPE, PROBE\_START, PROBE\_SEQUENCE, CHROMOSOME, PROBE\_CHR\_ORIENTATION, PROBE\_COORDINATES, DEFINITION, ONTOLOGY\_COMPONENT, ONTOLOGY\_PROCESS, ONTOLOGY\_FUNCTION, SYNONYMS, OBSOLETE\_PROBE\_ID. As the annotation data is huge, by default, we only input: ACCESSION, SYMBOL, PROBE\_START, CHROMOSOME, PROBE\_CHR\_ORIENTATION, PROBE\_COORDINATES, DEFINITION. This annotation information is kept in the featureData slot of ExpressionSet, which can be retrieved using pData(featureData(x.lumi)), suppose x.lumi is the LumiBatch object. As some annotation information may be outdated. We recommend using Bioconductor annotation packages to retrieve the annotation information.
}
\value{
return a LumiBatch object
}

\author{ Simon Lin, Pan Du }
\seealso{ \code{\link{LumiBatch}}, \code{\link{addNuID2lumi}} }

\examples{
## specify the file name
# fileName <- 'Barnes_gene_profile.txt' # Not Run
## load the data
# x.lumi <- lumiR(fileName)

## load the data with empty detection and beadNum slots
# x.lumi <- lumiR(fileName, columnNameGrepPattern=list(detection=NA, beadNum=NA))
}
\keyword{ IO }

