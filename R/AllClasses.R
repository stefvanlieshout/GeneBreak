#' An S4 class to represent a CopyNumberBreakPoints object.
#'
#' @slot segmDiff A matrix with breakpoints based on segment values
#' @slot callDiff A matrix with breakpoints based on call values
#' @slot segments A matrix with segmented copy number values
#' @slot calls A matrix with copy number calls
#' @slot featureAnnotation A dataframe with predefined information about the features (usually probes or bins)
#' @slot featureData A dataframe with calculated information about the features (usually probes or bins)
CopyNumberBreakPoints <- setClass( 
    'CopyNumberBreakPoints',
    
    slots = c(
        segmDiff     = "matrix",
        callDiff    = "matrix",
        calls       = "matrix",
        segments    = "matrix",
        breakpoints = "matrix",
        featureAnnotation  = "data.frame",
        featureData = "data.frame"
    ),
    prototype=list(
        segmDiff     = matrix(),
        callDiff    = matrix(),
        calls       = matrix(),
        segments    = matrix(),
        breakpoints = matrix(),
        featureAnnotation = data.frame(),
        featureData = data.frame()
    )
)

#' An S4 class to represent a CopyNumberBreakPointGenes object
#'
#' @slot geneAnnotation A data.frame with original gene annotation input
#' @slot geneData A data.frame with gene information added by package methods
#' @slot featuresPerGene A list with the associated features per gene
#' @slot breakpointsPerGene A matrix with breakage status per gene
CopyNumberBreakPointGenes <- setClass( 
    'CopyNumberBreakPointGenes',
    contains  = 'CopyNumberBreakPoints',
    slots = c(
        geneAnnotation = "data.frame",
        geneData = "data.frame",
        featuresPerGene = "list",
        breakpointsPerGene = "matrix"
    ),
    prototype=list(
        geneAnnotation = data.frame(),
        geneData = data.frame(),
        featuresPerGene = list(),
        breakpointsPerGene = matrix()
    )
)
# EOF