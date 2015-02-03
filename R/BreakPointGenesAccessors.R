
## ---------------
## define access methods information
## ---------------
.dataAccessOptions = list(
    callData = "feature call values",
    segmentData = "feature segment values",
    breakpointData = "feature breakpoint values",
    sampleNames = "vector with sample names",
    featureNames = "vector with feature names",
    featureAnnotation = "feature annotation",
    featureChromosomes = "vector of feature chromosomes",
    featureData = "feature data"
)
.dataAccessOptionsGene = list(
    geneAnnotation = "gene annotation",
    geneChromosomes = "vector of gene chromosomes",
    geneData = "gene data",
    featuresPerGene = "a list of genes with coupled features",
    breakpointsPerGene = "gene break status",
    recurrentGenes = "recurrently broken genes"
)

setMethod( "show",
    signature = "CopyNumberBreakPoints",
    definition = function(object) {
        
        cat( "\n --- Object Info ---\n", sep="")
        cat( " This is an object of class \"", class(object), "\"\n", sep = "" )
        cat( " ", nrow( object@segmDiff ), " features by ", ncol( object@segmDiff ), " samples\n", sep = "")
        cat( " A total of ", sum( object@breakpoints ), " breakpoints\n", sep = "" )
        if ( class(object) == "CopyNumberBreakPointGenes" && nrow(object@breakpointsPerGene) > 0){
            
            geneBreaksTotal <- sum( object@breakpointsPerGene )
            genesBrokenTotal <- length( which( rowSums( object@breakpointsPerGene ) > 0 ) )

            if ( !is.na(object@breakpointsPerGene)[1] ){
                cat( " A total of ", geneBreaksTotal, " gene breaks in ", genesBrokenTotal, " genes\n", sep = "" )
            }
            else{
                cat( " See ?bpGenes for how to determine gene breakpoints\n", sep = "" )
            }
            
            ## recurrent breakpoints and recurent genes
            ## no statictics run if no FDR data present at all
            geneFDRcount <- length( object@geneData$FDR )
            featFDRcount <- length( object@featureData$FDR )

            if ( geneFDRcount < 1 & featFDRcount < 1 ){    
                cat( " See ?bpStats for how to determine breakpoint statistics\n", sep = "" )
            }
            else{
                if ( length( object@geneData$FDR ) ){
                    signGenes <- length( which( object@geneData$FDR < 0.1 ) )
                    cat( " A total of ", signGenes, " recurrent breakpoint genes (FDR < 0.1)\n", sep = "" )
                }
                if ( length( object@featureData$FDR ) ){
                    signGenes <- length( which( object@featureData$FDR < 0.1 ) )
                    cat( " A total of ", signGenes, " recurrent breakpoints (FDR < 0.1)\n", sep = "" )
                }    
            }

        }
        cat( " See accessOptions(object) for how to access data in this object\n", sep = "" )

        cat( "\n" )
        invisible(NULL)
    }
)

#' Access Object Data
#' @param object An object of class \code{CopyNumberBreakPoints}
#' @examples
#' accessOptions( bp )
#' accessOptions( bp_genes )
#' accessOptions( bp_stats )
#' # --- Object data access ---
#' # This is an object of class "CopyNumberBreakPointGenes"
#' # callData( obj ) => returns feature call values
#' # segmentData( obj ) => returns feature segment values
#' # breakpointData( obj ) => returns feature breakpoint values
#' # sampleNames( obj ) => returns vector with sample names
#' # etc...
#' @aliases accessOptions
setMethod( "accessOptions",
    signature = "CopyNumberBreakPoints",
    definition = function(object) {
        cat( "\n --- Object data access ---\n", sep="")
        cat( " This is an object of class \"", class(object), "\"\n", sep = "" )
        for ( i in names(.dataAccessOptions) ){
            cat( " ", i, "( obj ) => returns ", .dataAccessOptions[[ i ]],"\n", sep="")    
        }
        ## extra options if class is CopyNumberBreakPointGenes
        if ( class(object) == "CopyNumberBreakPointGenes" ){
            for ( i in names(.dataAccessOptionsGene) ){
                cat( " ", i, "( obj ) => returns ", .dataAccessOptionsGene[[ i ]],"\n", sep="")    
            }   
        }
        cat( "\n" )
        invisible(NULL)
    }
)


## ---------------
## CopyNumberBreakPoints specific slot access
## ---------------
setMethod( "callData", "CopyNumberBreakPoints",
    function(object) object@calls
)
setMethod( "segmentData", "CopyNumberBreakPoints",
    function(object) object@segments
)
setMethod( "featureData", "CopyNumberBreakPoints",
    function(object) object@featureData
)
setMethod( "breakpointData", "CopyNumberBreakPoints",
    function(object) object@breakpoints
)
setMethod( "featureAnnotation", "CopyNumberBreakPoints",
    function(object) object@featureAnnotation
)
setMethod( "featureNames", "CopyNumberBreakPoints",
    function(object) rownames(object@breakpoints)
)
setMethod( "sampleNames", "CopyNumberBreakPoints",
    function(object) colnames(object@breakpoints)
)
setMethod( "featureChromosomes", "CopyNumberBreakPoints",
    function(object) object@featureAnnotation$Chromosome
)

## ---------------
## CopyNumberBreakPointGenes specific slot access
## ---------------
setMethod( "geneAnnotation", "CopyNumberBreakPointGenes",
    function(object) object@geneAnnotation
)
setMethod( "featuresPerGene", "CopyNumberBreakPointGenes",
    function(object, geneName=NULL){
        if ( !is.null(geneName) ){
            cat( "Gene chosen: ", geneName, "\n", sep="")    
            idx <- which( object@geneAnnotation$Gene == geneName)
            if( length(idx) == 0 ){
                stop( "Sorry, no record found for gene ", geneName, sep="")
            }
            object@featuresPerGene[[ idx ]]
        }
    } 
)
setMethod( "breakpointsPerGene", "CopyNumberBreakPointGenes",
    function(object) object@breakpointsPerGene
)
setMethod( "geneChromosomes", "CopyNumberBreakPointGenes",
    function(object) object@geneAnnotation$Chromosome
)
setMethod( "geneData", "CopyNumberBreakPoints",
    function(object) object@geneData
)

#' Show recurrent genes
#' @param object Output of bpStats(): a \code{CopyNumberBreakPointGenes} object
#' @param fdr.threshold Genes with lower FDR are returned
#' @param summarize If TRUE only certain columns are returned
#' @param order.column Name of the column to sort output on
#' @return data.frame with recurrent genes
#' @examples
#' recurrentGenes( bp_stats )
#' @aliases recurrentGenes
setMethod( "recurrentGenes", "CopyNumberBreakPointGenes",
    function(object, fdr.threshold=0.1, summarize=TRUE, order.column="FDR"){
        summaryColumns <- c("Gene","geneBreaks", "samplesWithGeneBreaks", "featureTotal", "pvalue", "FDR")
        if ( length( object@geneData$FDR ) ){
            idx <- which( object@geneData$FDR < fdr.threshold ) 
            cat( " A total of ", length(idx), " recurrent breakpoint genes (at FDR < ", fdr.threshold,")\n", sep = "" )
            output <- cbind( object@geneAnnotation[idx,], object@geneData[idx,] )
            output <- output[ order( output[ ,order.column], decreasing=FALSE ), ]
            if( summarize == TRUE ) output <- output[, summaryColumns]
            return( output )
        }else{
            cat( " No statistics information available in object, see ?bpStats\n", sep = "" )
        }
    } 
)

