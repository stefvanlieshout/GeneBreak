## contains the functions that create the BP object and show

#' Creates a BreakPointGenes object from a CGHcall object
#' 
#' @param someObj A "BreakPointGenes" object.
#' @param y Empty.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' plotDaShizzle( bpGenes )
readCGHcall <- function( cghCallObj, samples=c('A', 'B'), fillInt=9 ) 
{
    maxSamples <- 3

    if ( length( cghCallObj ) == 0L )
        stop( '[err in doSomething] someFile is empty...' )
    
    #if ( class( bins ) != 'numeric' )
    #    stop( '[err in doSomething] bins not of type numeric...' )

    if ( length( samples ) > maxSamples )
        stop( '[err in doSomething] too many samples, max is ', maxSamples )

    #samples = c('sampleA', 'sampleB')
    genes = c('geneA', 'geneB')
    counts <- matrix( fillInt, nrow=length(genes), ncol=length(samples), dimnames=list( genes, samples) )


    info <- data.frame( name=samples, row.names=samples )
    info$samples     <- samples
    info$total.reads <- rep( 1000, length(samples) )
    info$used.reads  <- rep( 990, length(samples) )

    output <- new( 'BreakPointGenes', counts=counts, meta=info, sampleNames=samples, geneNames=genes )
    return( output )
}

#' Does something, nice huh!
#' 
#' @param someObj A "BreakPointGenes" object.
#' @param y Empty.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' doSomething( bpGenes )
#' doSomething( bpGenes, nounou )
doSomething <- function( someObject, someFiles=NULL, someString='bam' ) 
{

    if ( length(someFiles) == 0L )
        stop('[err in doSomething] someFiles is empty...')
    
    #if (class(bins) == 'data.frame')
    #    bins <- AnnotatedDataFrame(bins)

    colN <- c('colA', 'colB')
    rowN <- c('rowA', 'rowB')
    counts <- matrix(NA_integer_, nrow=2, ncol=2, dimnames=list( rowN, colN) )


    data <- data.frame(name=colN, row.names=colN, stringsAsFactors=FALSE)
    data$total.reads <- 101
    data$used.reads <- 100

    output <- list()
    output$files <- someFiles
    output$data <- data
    output$counts <- counts
    # free mem
    #gc(FALSE)
    return( output )
    #new('ClassName', counts=counts, data=data)
}

divideNumbers <- function( number1, number2=5 ) 
{
    if (length(number1) == 0L)
        stop('[err in doSomething] input number1 is empty...')

    output <- number1 / number2
    return( output )
    #new('ClassName', counts=counts, data=data)
}

#' Sing a long!
#' 
#' @param msg Wha-ever
#' @param msg2 Sure why not
#' @return Returns NOTHING
#' @examples
#' sayHello()
#' sayHello( msg="First msg", msg2="First msg" )
sayHello <- function( msg="Is it me you're looking fooooor?", msg2="I can see it in your eyes..")
{
    cat( msg, "\n" )
    cat( msg2, "\n" )
}


#' Plots da SHIZZLE....
#' 
#' @param x A "BreakPointGenes" object.
#' @param y Empty.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' plotDaShizzle( bpGenes )
plotDaShizzle <- function( data, main="Title", col="red"){
    data <- log2( 1:10 )
    plot( data, main=main, col=col)
    lines( data, col="blue" )
}

createObject <- function(){
    return( 1:10 )
}

#' Multiplies the data
#' 
#' @param data    A "BreakPointGenes" object.
#' @param mFactor The multiply factor
#' @examples
#' multiplyCounts( bpGenes, 10 )
#' @aliases multiplyCounts
setMethod( "multiplyCounts", "BreakPointGenes",
    function(object, mFactor) {
    	if ( mFactor < 0 ){
    		stop("Sorry, multiplying with a negative number not allowed")
    	}
        object@counts <- object@counts * mFactor
        return(object)
    }
)

setMethod( "addMultipliedCounts", "BreakPointGenes",
    function(object, mFactor) {
        object <- as( object, "BreakPointGenes2" )
        object@mCounts <- object@counts * mFactor
        return(object)
    }
)

function( data, main="Title", col="red"){
    data <- log2( 1:10 )
    plot( data, main=main, col=col)
    lines( data, col="blue" )
}

# EOF
