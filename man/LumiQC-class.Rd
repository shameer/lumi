\name{LumiQC-class}
\docType{class}
\alias{LumiQC}
\alias{LumiQC-class}
\alias{[,LumiQC-method}
\alias{summary,LumiQC-method}
\alias{show,LumiQC-method}
\alias{initialize,LumiQC-method}
\alias{getHistory,LumiQC-method}


\title{Class LumiQC: quality control information of the samples (Illumina microarrays) }
\description{This is a class for the quality control of Illumina microarray. It extends \code{\link[Biobase]{ExpressionSet-class}}.}
\section{Extends}{
   Directly extends class \code{\link[Biobase]{ExpressionSet}}.
}

\section{Creating Objects}{
  The \code{LumiQC} object is created using the function \code{\link{lumiQ}} with \code{LumiBatch} object as input.}
}
\section{Slots}{
   Slot specific to \code{LumiQC}:
   \describe{
		\item{\code{mean}:}{ the mean of each sample}
		\item{\code{std}:}{ the standard deviation of each sample}
		\item{\code{sampleCor}:}{ the correlation bewteen the samples }
		\item{\code{detectionCount}:}{ the percentage of detected genes }
		\item{\code{sampleRelation}:}{ the sample relations (distance among samples) }
		\item{\code{outlier}:}{ identification of outlier samples. See \code{\link{detectOutlier}}. }
		\item{\code{history}:}{ the operation records of LumiQC objects. }
   }

   Slots inherited from \code{ExpressionSet}:
   \describe{
      \item{\code{assayData}:}{ contains equal dimensional matrices:
 		\code{exprs} (contains gene expression level, which is the mean of its bead replicates.),
		\code{cv} (coefficients of variance of each measurement).
 		For more details of \code{assayData}, please see \code{\link[Biobase]{ExpressionSet}} } 
  	}
}

\section{Methods}{

  \bold{Class-specific methods:}
  \describe{
		\item{\code{summary(LumiQC)}:}{ the summary of quality control information of the Illumina data }
		\item{\code{getHistory(LumiQC)}:}{ get the operation history of the \code{LumiQC} object }
		\item{\code{initialize(LumiQC)}:}{ object instantiation, used
		  by \code{new}; not to be called directly by the user.}    
    }

   \bold{Derived from \code{\link[Biobase]{ExpressionSet}}} (For the directly inherited methods, 
		please see \code{\link[Biobase]{ExpressionSet}} and \code{\link[Biobase]{eSet}}):
   \describe{
   		\item{\code{object[(,j)}:}{ conduct subsetting of data matrices and quality control information of samples in a LumiQC object.}
   }

   \bold{Standard generic methods} (For the directly inherited methods, 
	 please see \code{\link[Biobase]{ExpressionSet}} and \code{\link[Biobase]{eSet}}):
  \describe{
    \item{\code{initialize(LumiQC)}:}{ object instantiation, used
      	by \code{new}; not to be called directly by the user.}
    \item{\code{show(LumiQC)}:}{ the summary of quality control information of the Illumina data.}
   }
}

\details{
This class organize the quality control information of the Illumina data saved in a LumiBatch object. 
The LumiQC object can be produced by the function \code{\link{lumiQ}}.
}


\note{This class is better described in the vignette.}
\author{ Pan Du, Simon Lin }
\seealso{\code{\link{LumiBatch}}, \code{\link{lumiQ}}, \code{\link{plot.LumiQC}} }
\examples{
## load example data
data(example.lumi)

## Quality control estimation. lumi.Q is a LumiQC object
lumi.Q <- lumiQ(example.lumi)

## show the summary of the object
# summary(lumi.Q)
lumi.Q	

## get operation history
getHistory(lumi.Q)

## plot LumiQC object. 
## plot the density
plot(lumi.Q, what='density')

## plot the pairwise sample correlation
plot(lumi.Q, what='pair')

## sample relations
plot(lumi.Q, what='sampleRelation')

## See more examples in the help of plot.LumiQC 

}

\keyword{classes}
