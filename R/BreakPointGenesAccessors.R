
## --- define what is shown when object is called --- ##
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
        #systemUser <- system( "whoami", T )
        #cat( " Hi ", systemUser, "\n", sep = "")
        cat( "\n --- Object Info ---\n", sep="")
        cat( " This is an object of class ", class(object), "\n", sep = "" )
        cat( " ", nrow( object@segmDiff ), " features by ", ncol( object@segmDiff ), " samples.\n", sep = "")
        cat( " A total of ", sum( object@breakpoints ), " breakpoints\n", sep = "" )
        if ( class(object) == "CopyNumberBreakPointGenes" && nrow(object@breakpointsPerGene) > 0){
            geneBreaksTotal <- sum( object@breakpointsPerGene )
            genesBrokenTotal <- length( which( rowSums( object@breakpointsPerGene ) > 0 ) )

            if ( !is.na(object@breakpointsPerGene)[1] ){
                cat( " A total of ", geneBreaksTotal, " gene breaks in ", genesBrokenTotal, " genes\n", sep = "" )
            }
            else{
                cat( " Run bpGenes() to determine gene breakpoints\n", sep = "" )
            }
            
            ## recurrent breakpoints and recurent genes
            if ( length( object@geneData$FDR ) ){
                signGenes <- length( which( object@geneData$FDR < 0.1 ) )
                cat( " A total of ", signGenes, " recurrent breakpoint genes (FDR < 0.1)\n", sep = "" )
            }
            else{
                cat( " Run bpStats() to determine breakpoint statistics\n", sep = "" )
            }
            if ( length( object@featureData$FDR ) ){
                signGenes <- length( which( object@featureData$FDR < 0.1 ) )
                cat( " A total of ", signGenes, " recurrent breakpoints (FDR < 0.1)\n", sep = "" )
            }
        }
        
        cat( "\n --- Object Data Access Options ---\n", sep="")
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

## --- access to specific slots/data --- ##

## CopyNumberBreakPoints
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

## CopyNumberBreakPointGenes
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
#' @param object
#' @param fdr.threshold Genes with lower FDR are returned
#' @return data.frame with gene annotation and data
#' @examples
#' recurrentGenes( bpStats )
setMethod( "recurrentGenes", "CopyNumberBreakPointGenes",
    function(object, fdr.threshold=0.1){
        if ( length( object@geneData$FDR ) ){
            idx <- which( object@geneData$FDR < fdr.threshold ) 
            output <- cbind( object@geneAnnotation[idx,], object@geneData[idx,] )
            cat( " A total of ", length(idx), " recurrent breakpoint genes (at FDR < ", fdr.threshold,")\n", sep = "" )
            return( output )
        }else{
            cat( " No statistics information available in object, see ?bpStats\n", sep = "" )
        }
    } 
)

