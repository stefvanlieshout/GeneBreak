.runWorkflow <- function( object, geneAnnotation ) {
    startTime <- Sys.time()
    cat( "GeneBreak workflow started at: ", format(startTime),"\n", sep="" )
    
    ## input checks
    if ( (class( object ) != 'cghCall') && (class( object ) != 'QDNAseqSignals') )
        stop( '[ERR] input data not a cghCall or QDNAseqSignals object...' )

    ## perform all workflow steps
    bp <- getBreakpoints( data = object )
    bp <- bpFilter( object = bp )
    bp <- addGeneAnnotation( object = bp, geneAnnotation )
    bp <- bpGenes( bp )
    bp <- bpStats( bp )

    endTime <- Sys.time()
    timeDiff <- format( round( endTime - startTime, 3 ) )
    cat( "[Workflow runtime: ", timeDiff, "]\n", sep='')

    return(bp)
}