\name{rankinvariant}
\alias{rankinvariant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Rank Invariant Normalization }
\description{
  This function basically adjusts the samples to the same background level and then optionally scales to the same foreground level.
}
\usage{
rankinvariant(x.lumi, targetArray = NULL, rrc = .05, lowRank = seq(.5, .25, -.05), highRank = .9, minSize = .02, maxit=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{ an ExpressionSet inherited object or a data matrix with columns as samples and rows as genes }
  \item{targetArray}{ A target chip is the model for other chips to normalize. It can be a column index, a vector or a LumiBatch object with one sample. }
  \item{rrc}{ The relative rank change allowed for a gene to be selected as rank invariant }
  \item{lowRank}{ A vector with, in decreasing order, the minimum ranks where candidate genes can be selected as rank invariant }
  \item{highRank}{ The maximum rank where candidate genes can be selected as rank invariant }
  \item{minSize}{ Fraction of genes required to be selected as rank invariant }
  \item{maxit}{ Maximum number of iterations for \code{\link[MASS]{rlm}} to reach convergence }
}
\details{
Rank invariant normalization uses a set of genes that are rank invariant between a given sample and a target sample. The target sample can be predefined by setting the targetArray argument. If targetArray is NULL the average expression of all samples will be the target. Rank invariant genes are found for each sample seperately by calculation the relative rank change for each gene. Furthermore, only genes with ranks between the lowRank and highRank are considered. If the number of probes is less than minSize multiplies by the number of genes the next lowRank value tried. If no rank invariant set can be found an error is thrown.

The default settings of this function are the same as used Genomstudio (Illumina). The results produced by this method are similar, but not identical to Genomestudio.
}
\value{
Return an object with expression values normalized. The class of the return object is the same as the input object x.lumi.
}
\author{ Arno Velds (contact: a.velds (at) nki.nl) }
\seealso{ \code{\link{lumiN}} }
\keyword{ methods }
