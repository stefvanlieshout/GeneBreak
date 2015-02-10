#' An S4 class to represent a CopyNumberBreakPoints object.
#'
#' @slot segmDiff A matrix
#' @slot callDiff A matrix
#' @slot segments A matrix
#' @slot calls A matrix
#' @slot featureAnnotation A dataframe
#' @slot featureData A dataframe
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