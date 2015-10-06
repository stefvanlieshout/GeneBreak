#' An S4 class to represent a CopyNumberBreakPoints object.
#'
#' @slot segmDiff A matrix with breakpoints based on segment values
#' @slot callDiff A matrix with breakpoints based on call values
#' @slot segments A matrix with segmented copy number values
#' @slot calls A matrix with copy number calls
#' @slot featureAnnotation A dataframe with predefined information about the features (usually probes or bins)
#' @slot featureData A dataframe with calculated information about the features (usually probes or bins)
#' 
#' @section Accessors:
#' \itemize{
#'   \item callData( object ) Returns feature call values
#'   \item segmentData( object ) Returns feature segment values
#'   \item breakpointData( object ) Returns feature breakpoint values
#'   \item sampleNames( object ) Returns vector with sample names
#'   \item namesFeatures( object ) Returns vector with feature names
#'   \item featureChromosomes( object ) Returns vector of feature chromosomes
#'   \item featureInfo( object ) Returns feature data/information
#' }
#'
#' @section Methods:
#' \itemize{
#'   \item getBreakpoints Builds the \linkS4class{CopyNumberBreakPoints} object from copynumber data and detects breakpoint locations
#'   \item bpFilter Selects breakpoints by filter criteria options
#'   \item bpStats Applies cohort-based statistics to identify chromosomal locations that are recurrently affected by breakpoints
#'   \item bpPlot Plots breakpoint frequencies per chromosome
#' }
#' 
#' @author E. van den Broek and S. van Lieshout
#' @examples
#' data( copynumber.data.chr20 )
#' data( ens.gene.ann.hg18 )
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp )
#' bp <- bpStats( bp , level = 'feature' , method = 'BH' )
#' bpPlot( bp, c(20) )
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
#' 
#' @section Accessors:
#' \itemize{
#'    \item \code{callData( object )} Returns feature call values:
#'    \item \code{segmentData( object )} Returns feature segment values
#'    \item \code{breakpointData( object )} Returns feature breakpoint values
#'    \item \code{sampleNames( object )} Returns vector with sample names
#'    \item \code{namesFeatures( object )} Returns vector with feature names
#'    \item \code{featureChromosomes( object )} Returns vector of feature chromosomes
#'    \item \code{featureInfo( object )} Returns feature data/information
#'    \item \code{geneChromosomes( object )} Returns vector of gene chromosomes
#'    \item \code{geneInfo( object )} Returns gene data/information
#'    \item \code{featuresPerGene( object )} Returns a list of genes with coupled features
#'    \item \code{breakpointsPerGene( object )} Returns gene break status
#'    \item \code{recurrentGenes( object )} Returns recurrently broken genes
#' }
#'
#' @section Methods:
#' \itemize{
#'    \item getBreakpoints Builds the \linkS4class{CopyNumberBreakPoints} object from copynumber data and detects breakpoint locations
#'    \item bpFilter Selects breakpoints by filter criteria options
#'    \item addGeneAnnotation Maps features to gene locations
#'    \item bpGenes Indentifies genes affected by breakpoint locations
#'    \item bpStats Applies cohort-based statistics to identify genes and/or chromosomal locations that are recurrently affected by breakpoints
#'    \item bpPlot Plots breakpoint frequencies per chromosome
#' }
#' 
#' @author E. van den Broek and S. van Lieshout
#' @examples
#' data( copynumber.data.chr20 )
#' data( ens.gene.ann.hg18 )
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp )
#' bp <- addGeneAnnotation( bp, ens.gene.ann.hg18 )
#' bp <- bpGenes( bp )
#' bp <- bpStats( bp )
#' bpPlot( bp, c(20) )
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