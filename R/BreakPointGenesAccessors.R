
## --- define what is shown when object is called --- ##
.dataAccessOptions = list(
    segDiff  = "segmentDiff",
    callDiff = "callDiff",
    callData = "feature call values",
    segmentData = "feature segment values",
    breakpointData = "feature breakpoint values",
    sampleNames = "vector with sample names",
    featureNames = "vector with feature names",
    featureAnnotation = "feature annotation"
)
.dataAccessOptionsGene = list(
    geneAnnotation = "gene annotation",
    featuresPerGene = "a list of genes with coupled features",
    breakpointsPerGene = "gene break status"
)

setMethod( "show",
    signature = "CopyNumberBreakPoints",
    definition = function(object) {
        #systemUser <- system( "whoami", T )
        #cat( " Hi ", systemUser, "\n", sep = "")
        cat( "\n --- Object Info ---\n", sep="")
        cat( " This is an object of class ", class(object), "\n", sep = "" )
        cat( " ", nrow( object@segDiff ), " features by ", ncol( object@segDiff ), " samples.\n", sep = "")
        cat( " A total of ", sum( object@breakpoints ), " breakpoints\n", sep = "" )
        if ( class(object) == "CopyNumberBreakPointGenes" && nrow(object@breakpointsPerGene) > 0){
            geneBreaksTotal <- sum( object@breakpointsPerGene )
            genesBrokenTotal <- length( which( rowSums( object@breakpointsPerGene ) > 0 ) )
            cat( " A total of ", geneBreaksTotal, " gene breaks in ", genesBrokenTotal, " genes\n", sep = "" )
        }
        
        cat( "\n --- Object Data Access ---\n", sep="")
        for ( i in names(.dataAccessOptions) ){
            cat( " ", i, "(obj) => returns ", .dataAccessOptions[[ i ]],"\n", sep="")    
        }
        
        ## extra options if class is CopyNumberBreakPointGenes
        if ( class(object) == "CopyNumberBreakPointGenes" ){
            for ( i in names(.dataAccessOptionsGene) ){
                cat( " ", i, "(obj) => returns ", .dataAccessOptionsGene[[ i ]],"\n", sep="")    
            }   
        }

        cat( "\n" )
        invisible(NULL)
    }
)

## --- access to specific slots/data --- ##

## CopyNumberBreakPoints
setMethod( "segDiff", "CopyNumberBreakPoints",
	function(object) object@segDiff
)
setMethod( "callDiff", "CopyNumberBreakPoints",
	function(object) object@callDiff
)
setMethod( "callData", "CopyNumberBreakPoints",
	function(object) object@calls
)
setMethod( "segmentData", "CopyNumberBreakPoints",
	function(object) object@segments
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

## CopyNumberBreakPointGenes
setMethod( "geneAnnotation", "CopyNumberBreakPointGenes",
    function(object) object@geneAnnotation
)
setMethod( "featuresPerGene", "CopyNumberBreakPointGenes",
    function(object, geneName=NULL){
        if ( !is.null(geneName) ){
            cat( "Gene chosen:", geneName, "\n", sep="")    
            idx <- which( object@geneAnnotation$Gene == geneName)
            if( length(idx) == 0 ){
                stop( "Sorry, no record found for gene ", geneName, sep="")
            }
            object@featuresPerGene[idx]
        }
        object@featuresPerGene
    } 
)
setMethod( "breakpointsPerGene", "CopyNumberBreakPointGenes",
    function(object) object@breakpointsPerGene
)

