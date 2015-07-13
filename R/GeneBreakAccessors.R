
## ---------------
## define access methods information
## ---------------
.dataAccessOptions = list(
    callData = "feature call values",
    segmentData = "feature segment values",
    breakpointData = "feature breakpoint values",
    sampleNames = "vector with sample names",
    featureNames = "vector with feature names",
    featureChromosomes = "vector of feature chromosomes",
    featureInfo = "feature data/information"
)
.dataAccessOptionsGene = list(
    geneChromosomes = "vector of gene chromosomes",
    geneInfo = "gene data/information",
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
        ## --- this does not work here while it does in terminal R ---
        #inputVariableName <- deparse(substitute(object))
        #cat( " See accessOptions(",inputVariableName,") for how to access data in this object\n", sep = "" )
        ## --- 
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
#' @aliases accessOptions
setMethod( "accessOptions",
    signature = "CopyNumberBreakPoints",
    definition = function(object) {
        cat( "\n --- Object data access ---\n", sep="")
        cat( " This is an object of class \"", class(object), "\"\n", sep = "" )
        for ( i in names(.dataAccessOptions) ){
            cat( " ", i, "( object ) => returns ", .dataAccessOptions[[ i ]],"\n", sep="")    
        }
        ## extra options if class is CopyNumberBreakPointGenes
        if ( class(object) == "CopyNumberBreakPointGenes" ){
            for ( i in names(.dataAccessOptionsGene) ){
                cat( " ", i, "( object ) => returns ", .dataAccessOptionsGene[[ i ]],"\n", sep="")    
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
setMethod( "breakpointData", "CopyNumberBreakPoints",
    function(object) object@breakpoints
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

#' Show recurrent genes
#' @param object Output of bpStats(): a \code{CopyNumberBreakPointGenes} object
#' @param fdr.threshold Genes with lower FDR are returned
#' @param summarize If TRUE only certain columns are returned
#' @param order.column Name of the column to sort output on
#' @return data.frame with recurrent genes
#' @examples
#' recurrentGenes( breakpointStatistics )
#' @aliases recurrentGenes
setMethod( "recurrentGenes", "CopyNumberBreakPointGenes",
    function(object, fdr.threshold=0.1, summarize=TRUE, order.column="FDR"){
        summaryColumns <- c("Gene","sampleCount", "featureTotal", "pvalue", "FDR")
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

#' Gather Gene information
#' @param object of class \code{CopyNumberBreakPointGenes}
#' @return data.frame
#' @examples
#' geneInfo( bp_genes )
#' @aliases geneInfo
setMethod( "geneInfo", "CopyNumberBreakPointGenes",
    function( object ){
        ## the gene info is located in two data.frames
        annCount <- nrow(object@geneAnnotation)
        datCount <- nrow(object@geneData)
        if( annCount != datCount ){
            stop( "Somehow different amount of rows in annotations [", annCount, "] and data [", datCount, "]\n" )
        }
        output <- cbind( object@geneAnnotation, object@geneData )
        return( output )
    } 
)
#' Gather Feature information
#' @param object of class \code{CopyNumberBreakPoints}
#' @return data.frame
#' @examples
#' featureInfo( bp )
#' @aliases featureInfo
setMethod( "featureInfo", "CopyNumberBreakPoints",
    function( object ){
        ## the gene info is located in two data.frames of same nrow
        annCount <- nrow(object@featureAnnotation)
        datCount <- nrow(object@featureData)
        if( annCount != datCount ){
            stop( "Somehow different amount of rows in annotations [", annCount, "] and data [", datCount, "]\n" )
        }
        output <- cbind( object@featureAnnotation, object@featureData )
        return( output )
    } 
)