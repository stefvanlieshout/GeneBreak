BreakPointGenes <- setClass( 
	# name of class
	'BreakPointGenes',
	# contents
	slots = c(
		counts = "matrix",
		meta = "list",
		sampleNames = "character",
		geneNames = "character"
	),
	# set defaults
	prototype=list(
		counts = matrix(),
		meta = list(),
		sampleNames = c(),
		geneNames = c()
	)
)

BreakPointGenes2 <- setClass( # subclass
	# name of class
	'BreakPointGenes2',
	# is subclass of
	contains  = 'BreakPointGenes',
	# contents
	slots = c(
		mCounts = "matrix"
	),
	# set defaults
	prototype=list(
		mCounts = matrix()
	)
)

#' An S4 class to represent a CopyNumberBreakPointGenes object.
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

CopyNumberBreakPointGenes <- setClass( 
	'CopyNumberBreakPointGenes', # name of classs
	contains  = 'CopyNumberBreakPoints',
	slots = c( # contents
		geneAnnotation = "data.frame",
		geneProbes = "list",
		geneBreakpoints = "matrix"
	),
	prototype=list( # set defaults
		geneAnnotation = data.frame(),
		geneProbes = list(),
		geneBreakpoints = matrix()
	)
)

# EOF