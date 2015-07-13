#' GeneBreak: A package for gene breakpoint detection on copy number abberation data
#'
#' The GeneBreak package performs cohort based recurrent 
#' gene breakpoint detection on copynumber data. It is possible
#' to use the output of the function \code{\link[CGHcall]{CGHcall}} from 
#' the package \code{CGHcall} or the function \code{\link[QDNAseq]{callBins}} from the
#' package \code{QDNAseq} as the input for this analysis.
#' 
#' @section GeneBreak functions:
#' Analysis starts with the function \code{\link{getBreakpoints}} and continues with:\cr
#' \code{\link{bpFilter}} to exclude certain breakpoints from the analysis\cr
#' \code{\link{addGeneAnnotation}} to add gene location information\cr
#' \code{\link{bpGenes}} to determine which features (probes/bins) are related to which genes\cr
#' \code{\link{bpStats}} to determine which gene breaks are recurrent in the cohort\cr
#'
#' @docType package
#' @name GeneBreak
NULL
#> NULL