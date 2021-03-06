% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/GeneBreakMethods.R
\name{getBreakpoints}
\alias{getBreakpoints}
\title{getBreakpoints}
\usage{
getBreakpoints(data, data2 = NULL, first.rm = TRUE)
}
\arguments{
\item{data}{An object of class \pkg{cghCall} or an object of class QDNAseqCopyNumbers or a data.frame containing feature annotations ("Chromosome", "Start", "End", "FeatureName") followed by copy number segment values (rows are features, columns are subjects).}

\item{data2}{A "data.frame" containing copy number calls following feature annotations with the four columns ("Chromosome", "Start", "End", "FeatureName", ...). This is optional and allows CNA-associated breakpoint filtering. (see ?bpFilter)}

\item{first.rm}{Remove the first 'artificial' breakpoint of the first DNA segment for each chromosome (default: first.rm=TRUE)}
}
\value{
Returns an object of class \linkS4class{CopyNumberBreakPoints}.
}
\description{
Builds the \linkS4class{CopyNumberBreakPoints} object from copynumber data and detects breakpoint locations.
}
\details{
The accuracy of chromosomal breakpoint locations depends on the quality and genomic resolution of processed copy number data. For CNA input data, we recommend to use established computational methods for CNA detection such as 'CGHcall' (Van De Wiel et al., 2007) for array-CGH or 'QDNAseq' (Scheinin et al., 2014) for MPS data, which both use the implemented Circular Binary Segmentation algorithm (Olshen et al. 2004).
}
\examples{
data( copynumber.data.chr20 )
breakpoints <- getBreakpoints( data = copynumber.data.chr20 )

## or alternatively
library(CGHcall)
cgh <- copynumber.data.chr20
segmented <- data.frame( Chromosome=chromosomes(cgh), Start=bpstart(cgh),
  End=bpend(cgh), FeatureName=rownames(cgh), segmented(cgh))
called <- data.frame( Chromosome=chromosomes(cgh), Start=bpstart(cgh),
  End=bpend(cgh), FeatureName=rownames(cgh), calls(cgh))
breakpoints <- getBreakpoints( data = segmented, data2 = called )

## options to inspect the data
breakpoints
accessOptions( breakpoints )
}
\references{
Van De Wiel, M.A. et al. (2007) CGHcall: calling aberrations for array CGH tumor profiles. Bioinformatics, 23, 892-894.

Scheinin, I. et al. (2014) DNA copy number analysis of fresh and formalin-fixed specimens by shallow whole-genome sequencing with identification and exclusion of problematic regions in the genome assembly. Genome Research, 24, 2022-2032.

Olshen, A.B. et al. (2004) Circular binary segmentation for the analysis of array-based DNA copy number data. 5, 557-572.
}

