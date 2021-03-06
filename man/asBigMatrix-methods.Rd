\name{asBigMatrix-methods}
\docType{methods}
\alias{asBigMatrix-methods}
\alias{asBigMatrix}
\alias{asBigMatrix,ExpressionSet-method}
\title{convert the data matrix in the assayData of a ExpressionSet as BigMatrix}
\description{
 convert the data matrix in the assayData of a ExpressionSet as BigMatrix
}
\usage{
\S4method{asBigMatrix}{ExpressionSet}(object, rowInd=NULL, colInd=NULL, nCol=NULL, dimNames=NULL, saveDir='.', savePrefix=NULL, ...)
}

\arguments{
\item{object}{ an object of \code{\link[Biobase]{ExpressionSet}} or its inherited class }
\item{rowInd}{ the subset of row index }
\item{colInd}{ the subset of column index }
\item{nCol}{ the number of columns of the data, which can be larger than the real data dimension. It is designed for adding future data.}
\item{dimNames}{ the dimension names, which is a list of two character vectors (rownames and colnames)}
\item{saveDir}{the parent directory to save the BigMatrix data files}
\item{savePrefix}{the folder name prefix of the directory to save the BigMatrix data files. The fold name will be like this: paste(savePrefix, '_bigmat', sep='') }
\item{\dots}{ optional arguments to \code{\link[bigmemoryExtras]{BigMatrix}} }
}
\details{
This function does not work in Windows because the dependent package bigmemoryExtras does not support it. 
In order to make lumi package still compilation under Windows, I deliberately remove the dependency of bigmemoryExtras package.
As a result, users need to manually load the bigmemoryExtras function before using this function.

The BigMatrix data files will be save in the directory file.path(saveDir, paste(savePrefix, '_bigmat', sep=''))
}
\seealso{
\code{\link[bigmemoryExtras]{BigMatrix}}
}

\keyword{methods}
