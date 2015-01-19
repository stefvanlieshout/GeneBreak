## ---------------
## eSet way
## ---------------
# CNBP <-
# 	function(
# 		phenoData=AnnotatedDataFrame(), 
# 		experimentData=MIAME(),
# 		annotation=character(), 
# 		segDiff=new("matrix"), 
# 		callDiff=new("matrix"),
# 		calls=new("matrix"), 
# 		segments=new("matrix"), 
# 		breakpoints=new("matrix"), 
# 		...
# 	)
# {
# 	.CNBP( 
# 		phenoData=phenoData, 
# 		experimentData=experimentData,
# 		annotation=annotation, R=R, G=G, Rb=Rb, Gb=Gb, 
# 		...
# 	)
# }

# setValidity( "CNBP", function(object) {
# 	assayDataValidMembers( assayData(object), c("R", "G", "Rb", "Gb"))
# })


## ---------------
## Own way of things
## ---------------
#' An S4 class to represent a CopyNumberBreakPoints object.
#'
#' @slot segDiff A matrix
#' @slot callDiff A matrix
#' @slot segments A matrix
#' @slot calls A matrix
#' @slot annotation A dataframe
CopyNumberBreakPoints <- setClass( 
	'CopyNumberBreakPoints', # name of classs
	
	slots = c( # contents
		segDiff     = "matrix",
		callDiff    = "matrix",
		calls       = "matrix",
		segments    = "matrix",
		breakpoints = "matrix",
		featureAnnotation  = "data.frame"
	),
	prototype=list( # set defaults
		segDiff     = matrix(),
		callDiff    = matrix(),
		calls       = matrix(),
		segments    = matrix(),
		breakpoints = matrix(),
		featureAnnotation = data.frame()
	)
)

#' An S4 class to represent a CopyNumberBreakPointGenes object
#'
#' @slot geneAnnotation A data.frame with gene information
#' @slot featuresPerGene A list with per gene the associated features
#' @slot breakpointsPerGene A matrix with breakage status per gene
CopyNumberBreakPointGenes <- setClass( 
	'CopyNumberBreakPointGenes', # name of classs
	contains  = 'CopyNumberBreakPoints',
	slots = c( # contents
		geneAnnotation = "data.frame",
		featuresPerGene = "list",
		#geneBreakStatus = "matrix",
		breakpointsPerGene = "matrix"
	),
	prototype=list( # set defaults
		geneAnnotation = data.frame(),
		featuresPerGene = list(),
		#geneBreakStatus = matrix(),
		breakpointsPerGene = matrix()
	)
)

# EOF