#' An S4 class to represent a CopyNumberBreakPoints object.
#'
#' @slot segmDiff A matrix 
#' @slot callDiff A matrix
#' @slot segments A matrix
#' @slot calls A matrix
#' @slot annotation A dataframe
CopyNumberBreakPoints <- setClass( 
	'CopyNumberBreakPoints', # name of classs
	
	slots = c( # contents
		segmDiff     = "matrix",
		callDiff    = "matrix",
		calls       = "matrix",
		segments    = "matrix",
		breakpoints = "matrix",
		featureAnnotation  = "data.frame",
		featureData = "data.frame"
	),
	prototype=list( # set defaults
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
#' @slot geneAnnotation A data.frame with originally given gene information
#' @slot geneData A data.frame with gene information added by package methods
#' @slot featuresPerGene A list with per gene the associated features
#' @slot breakpointsPerGene A matrix with breakage status per gene
CopyNumberBreakPointGenes <- setClass( 
	'CopyNumberBreakPointGenes', # name of classs
	contains  = 'CopyNumberBreakPoints',
	slots = c( # contents
		geneAnnotation = "data.frame",
		geneData = "data.frame",
		featuresPerGene = "list",
		breakpointsPerGene = "matrix"
	),
	prototype=list( # set defaults
		geneAnnotation = data.frame(),
		geneData = data.frame(),
		featuresPerGene = list(),
		breakpointsPerGene = matrix()
	)
)
# EOF