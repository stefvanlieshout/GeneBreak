.SEP_CHAR <- ','

#' getBreakpoints
#' 
#' @description
#' Builds the \linkS4class{CopyNumberBreakPoints} object from copynumber data and detects breakpoint locations.
#' @param data An object of class \pkg{cghCall} or an object of class QDNAseqCopyNumbers or a data.frame containing feature annotations ("Chromosome", "Start", "End", "FeatureName") followed by copy number segment values (rows are features, columns are subjects).
#' @param data2 A "data.frame" containing copy number calls following feature annotations with the four columns ("Chromosome", "Start", "End", "FeatureName", ...). This is optional and allows CNA-associated breakpoint filtering. (see ?bpFilter)
#' @param first.rm Remove the first 'artificial' breakpoint of the first DNA segment for each chromosome (default: first.rm=TRUE)
#' @return Returns an object of class \linkS4class{CopyNumberBreakPoints}.
#' @details The accuracy of chromosomal breakpoint locations depends on the quality and genomic resolution of processed copy number data. For CNA input data, we recommend to use established computational methods for CNA detection such as 'CGHcall' (Van De Wiel et al., 2007) for array-CGH or 'QDNAseq' (Scheinin et al., 2014) for MPS data, which both use the implemented Circular Binary Segmentation algorithm (Olshen et al. 2004).
#' @references Van De Wiel, M.A. et al. (2007) CGHcall: calling aberrations for array CGH tumor profiles. Bioinformatics, 23, 892-894.
#' @references Scheinin, I. et al. (2014) DNA copy number analysis of fresh and formalin-fixed specimens by shallow whole-genome sequencing with identification and exclusion of problematic regions in the genome assembly. Genome Research, 24, 2022-2032.
#' @references Olshen, A.B. et al. (2004) Circular binary segmentation for the analysis of array-based DNA copy number data. 5, 557-572.
#' @examples
#' data( copynumber.data.chr20 )
#' breakpoints <- getBreakpoints( data = copynumber.data.chr20 )
#'
#' ## or alternatively
#' library(CGHcall)
#' cgh <- copynumber.data.chr20
#' segmented <- data.frame( Chromosome=chromosomes(cgh), Start=bpstart(cgh), 
#'   End=bpend(cgh), FeatureName=rownames(cgh), segmented(cgh))
#' called <- data.frame( Chromosome=chromosomes(cgh), Start=bpstart(cgh), 
#'   End=bpend(cgh), FeatureName=rownames(cgh), calls(cgh))
#' breakpoints <- getBreakpoints( data = segmented, data2 = called )
#' 
#' ## options to inspect the data
#' breakpoints
#' accessOptions( breakpoints )
getBreakpoints <- function( data, data2=NULL, first.rm=TRUE ) {
    
    ## input checks
    if ( (class( data ) == 'cghCall') || (class( data ) == 'QDNAseqCopyNumbers') ){
	    ## slots are the same for CGHcall and QDNAseq objects
	    #segmData <- data@assayData$segmented
	    #bpChrs   <- data@featureData@data$Chromosome
	    #bpStart  <- data@featureData@data$Start
	    #bpEnd    <- data@featureData@data$End
	    segmData <- segmented(data)
	    bpChrs   <- chromosomes(data)
	    bpStart  <- bpstart(data)
	    bpEnd    <- bpend(data)
	    featureNames <- rownames( segmData )
	    callData <- matrix( data=NA, ncol=ncol(segmData), nrow=nrow(segmData) )

	    ## check slot
	    if ( exists( "calls", data@assayData ) ) callData <- calls(data)

    }
    else if ( class(data) == "data.frame" ){
    	requiredColnames <- c( "Chromosome", "Start", "End", "FeatureName" )
    	## check if all required columns are present in data
    	for (col in requiredColnames){
    		if ( ! col %in% colnames(data) ){
    			stop( '[ERR] input data is a data.frame but missing column "', col, '"' )		
    		}
    	}

    	## all ok: prepare segmented data
    	bpChrs <- data$Chromosome
    	bpStart <- data$Start
    	bpEnd <- data$End
    	featureNames <- data$FeatureName
    	segmData <- as.matrix( data[, 5:ncol(data) ] )
    	callData <- matrix( data=NA, ncol=ncol(segmData), nrow=nrow(segmData) )

    	## now check if data2 is set and setup call data
    	if ( !is.null(data2) ){
    		if ( class(data2) != "data.frame" ) stop( '[ERR] param data2 can only be a data.frame')
    		if ( ! identical( data[,1:4], data2[, 1:4] ) ) stop( '[ERR] data and data2 differ in annotation columns (column 1 to 4)')
    		callData <- as.matrix( data2[, 5:ncol(data2) ] )
    	}
    	rownames( segmData ) <- featureNames
    	rownames( callData ) <- featureNames
    }
    else {
    	stop( '[ERR] input data not a cghCall/QDNAseqSignals object nor a data.frame...' )
    }

    cat( "Breakpoint detection started...\n", sep="" )


    breakpoints <- NULL; segmDiff <- NULL; callDiff <- NULL; startChr <- c()
    featureAnnotation <- NULL;  featureInterval <- c()
    
    segmDiff <- rbind( segmData[ 1, ], apply( segmData, 2, diff ) ) 
    callDiff <- rbind( callData[ 1, ], apply( callData, 2, diff ) )
    rownames(segmDiff) <- featureNames
    rownames(callDiff) <- featureNames
    
    ## remove first chromosomal BP if requested: set bp value to 0
    if ( first.rm == TRUE ) {
        startChr <- c( 1, which( diff( bpChrs ) == 1 ) + 1 )
        segmDiff[ startChr, ] <- 0
        callDiff[ startChr, ] <- 0
    }
    
    featureAnnotation <- data.frame( 
        Chromosome = bpChrs, 
        Start = bpStart,
        End = bpEnd,
        row.names = featureNames 
    )

    ## probe distance is needed for statistics
    featureInterval <- c( 0, diff( bpStart ) )
    featureInterval[ startChr ] <- 0
    featureData <- data.frame( 
        featureInterval = featureInterval, 
        row.names = featureNames 
    )
    breakpoints <- ifelse( segmDiff != 0, 1, 0 )

    ## NOTE: still to change name!! and remove in bpStats
    featureData$sampleCount <- apply( breakpoints, 1, sum )

    samplesListed <- apply( breakpoints, 1, function(x){ (names(x)[ x>0 ]) })
    featureData$sampleNamesWithBreakpoints <- sapply( samplesListed, function( x ) { 
        ifelse( any( is.na(x) ), NA, paste( x, collapse = .SEP_CHAR ) )
    } )

    ## create object with all required slots
    output <- new( 'CopyNumberBreakPoints', 
        segmDiff = segmDiff,
        callDiff = callDiff,
        featureAnnotation = featureAnnotation,
        featureData = featureData,
        calls = callData,
        segments = segmData,
        breakpoints = breakpoints
    )
    return( output )
}

#' bpFilter
#'
#' @description
#' Selects breakpoints by filter criteria options.
#' @param object An object of class \linkS4class{CopyNumberBreakPoints}
#' @param filter Type of filter. This can be either "CNA-ass", "deltaSeg" or "deltaCall".
#' \itemize{
#'   \item CNA-ass: filter out breakpoints that are flanked by copy number neutral segments to obtain CNA-associated breakpoint locations
#'   \item deltaSeg: selects for breakpoints where the log2 ratio transition of the copy number segments exceeds the user-defined threshold
#'   \item deltaCall: selects only breakpoints of discrete copy number states (amplification, gain, neutral, loss)
#' }
#' @param threshold Set the minimal log2 ratio difference between segments. This parameter is required for the "deltaSeg" filter option
#' @return Returns an object of class \linkS4class{CopyNumberBreakPoints} with breakpoint matrix replaced by filtered breakpoints.
#' @details Filter options "CNA-ass" and "deltaCall" require calls in addition to segmented copynumber data (see input for \code{getBreakpoints()} )
#' @examples
#' data( copynumber.data.chr20 )
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp, filter = "CNA-ass" )
#' bp <- bpFilter( bp, filter = "deltaSeg", threshold = 0.2 )
#' 
#' ## options to inspect the data
#' bp
#' accessOptions( bp )
#' @aliases bpFilter
setMethod( "bpFilter", "CopyNumberBreakPoints",
    function( object, filter="CNA-ass", threshold=NULL) {

        cat("Applying BP selection...\n")
        allowed.filters <- c( "deltaSeg", "deltaCall", "CNA-ass" )

        ## Check input params
        if ( ! filter %in% allowed.filters ){
            stop( "Parameter \"filter\" incorrect, options are: ", paste( allowed.filters, collapse=', ' ), "\n" )
        }
        if ( !is.null(threshold) & !is.numeric( threshold ) ){
            stop( "Parameter \"threshold\" incorrect, should be numeric", "\n" )
        }

        ## filter breakpoints based on type of filter
        if( filter == "deltaSeg" ) {
            if( !is.null( threshold ) ){
                breakpoints_filtered <- ifelse( abs( object@segmDiff ) > threshold, 1, 0 )
            }else{
                stop( "parameter \"threshold needed\" when deltaSeg chosen as filter" )
            }
        }
        else if( filter == "deltaCall" ) {
            breakpoints_filtered <- ifelse( object@callDiff != 0, 1, 0 )
        }
        else if( filter == "CNA-ass" ) {
            breakpoints_filtered <- ifelse( 
                object@segmDiff != 0, 
                ifelse( object@callDiff != 0 | object@calls != 0, 1, 0), 
                0
            )
        }
        else{
            stop( "parameter \"filter\" has incorrect value" )
        }
        
        ## replace breakpoint data
        object@breakpoints <- breakpoints_filtered
        
        ## recalculate the number of samples that have a breakpoint after filtering
        object@featureData$sampleCount <- apply( object@breakpoints, 1, sum )
        
        ## same for samplenames
        samplesListed <- apply( object@breakpoints, 1, function(x){ (names(x)[ x>0 ]) })
        object@featureData$sampleNamesWithBreakpoints <- sapply( samplesListed, function( x ) { 
          ifelse( any( is.na(x) ), NA, paste( x, collapse = .SEP_CHAR ) )
        } )
        
        ## return updated object
        return( object )
    }
)

#' addGeneAnnotation
#' 
#' @description
#' Maps features to gene locations.
#' @param object An object of class \linkS4class{CopyNumberBreakPoints}
#' @param geneAnnotation An object of class \code{GRanges} or dataframe with at least four columns ("Gene", "Chromosome", "Start", "End")
#' @return Returns an object of class \linkS4class{CopyNumberBreakPointGenes} with gene annotation added.
#' @details The end of the first feature after gene start location up to and including the first feature after gene end location will be defined as gene-associated feaures. For hg18, hg19 and hg38 built-in gene annotation files obtained from ensembl can be used. Please take care to use a matching reference genome for your breakpoint data. In stead of using the built-in gene annotion files, feature-to-gene mapping can be based on an user-defined annotion file. The dataframe should contain at least these four columns: "Gene", "Chromosome", "Start" and "End".
#' @details Note that the chromosome names are labeled similar in the annotation file and in the object of class \linkS4class{CopyNumberBreakPoints} (e.g. "1", "2", ..., "24")
#' @examples
#' data( copynumber.data.chr20 )
#' data( ens.gene.ann.hg18 )
#'
#' ## other buil-in gene annotations are:
#' # data( ens.gene.ann.hg19 )
#' # data( ens.gene.ann.hg38 )
#'
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp ) 
#' # input copynumber.data.chr20 is hg18 based
#' bp <- addGeneAnnotation( bp, ens.gene.ann.hg18 )
#' 
#' ## options to inspect the data
#' bp
#' accessOptions( bp )
#' @aliases addGeneAnnotation
setMethod( "addGeneAnnotation", "CopyNumberBreakPoints",
    function( object, geneAnnotation ) {

    ## input checks
    if ( (class( object ) == 'GRanges') ){

    	geneAnnotationExt <- as.data.frame( geneAnnotation, row.names=NULL)[,grep("start|end|seqnames", colnames(as.data.frame(geneAnnotation)), invert=TRUE)] 

	    geneNames <- names(geneAnnotation)
	    bpChrs    <- gsub( "chr", "", seqnames(geneAnnotation) )
	    bpStart   <- start(geneAnnotation)
	    bpEnd     <- end(geneAnnotation)
	    
	    geneAnnotationCore <- data.frame( 
	    	Gene = geneNames,
	    	Chromosome = bpChrs,
	    	Start = bpStart,
	    	End = bpEnd
	    )
	    geneAnnotation <- cbind( geneAnnotationCore, geneAnnotationExt)
    }

    ## sanity checks
    requiredColnames <- c("Gene", "Chromosome", "Start", "End")
    for (colName in requiredColnames ){
        if( ! colName %in% colnames(geneAnnotation) ){
            stop( "Missing column ", colName, " in geneAnnotation dataframe", "\n" )
        }
    }

    ## check data chromosomes and remove others from gene annotation
    dataChrs <- names(table( featureChromosomes(object) ))
    geneAnnotation <- geneAnnotation[ which( geneAnnotation$Chromosome %in% dataChrs ), ]

    ## setup gene information
    geneData <- data.frame( geneLength = apply( geneAnnotation[, c("Start","End") ], 1, diff ) )
    geneData$remarks <- NA
    geneData$genelength_features <- NA
    geneData$featureTotal <- NA
    geneData$featureNames <- NA

    geneCount <- nrow( geneData )
    sampCount <- length( sampleNames(object) )
    ## setup feature (probe/bin) information
    features <- cbind( object@featureAnnotation, idx=1:nrow( object@featureAnnotation ) )
    ## make list containing geneAnnotation with related featureAnnotation
    gene_features <- list()

    ## define index of first and last chromosomal probe
    start_endChr <- data.frame(
        chr         = features$Chromosome[ c( 1, which( diff( features$Chromosome ) == 1 ) +1 ) ],
        start_index = c( 1, which( diff( features$Chromosome ) == 1 ) +1 ),
        end_index   = c( which( diff( features$Chromosome ) == 1 ), length( features$Chromosome ) ) 
    )

    ## Start analysis
    cat( paste( "Adding of gene annotation started on ", geneCount, " genes by ", sampCount, " samples\n", sep=""))

    progress <- rep( NA, geneCount )
    progress[ c( round( seq( 1, geneCount, by = ( geneCount/4))))] <- c("0%","25%","50%","75%")
    warning_Z <- c() # will be filled if unable to determine situation for one or more genes

    ## Loop for all geneAnnotation
    for( gene_idx in 1:geneCount ) {

        if( !is.na(progress[gene_idx])) { 
            cat( paste( progress[ gene_idx ], "... ") ) 
        }

        feature_idx_chr <- which( features$Chromosome == geneAnnotation$Chromosome[ gene_idx ] )
        tmp_features <- features[ feature_idx_chr, ]
        features_S <- NULL  # first probe after start position gene
        features_E <- NULL  # first probe after end position gene

        features_S <- head( tmp_features$idx[ which(tmp_features$End >= geneAnnotation$Start[ gene_idx ])], 1) # for QDNAseq: take features$Start too; CAVE !!!
        features_E <- head( tmp_features$idx[ which(tmp_features$Start >= geneAnnotation$End[ gene_idx ])], 1)

        geneRelatedFeatures <- NULL

        if( any(features_S) ) {

            if( !any(features_E) ) {                

                features_E <- start_endChr[ which( start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ] ), "end_index"] # make features_E last probe on chr
                geneRelatedFeatures <- rownames(features)[ features_S:features_E ]

                geneData$remarks[gene_idx] <- "E"
                geneData$genelength_features[ gene_idx ] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
                gene_features[[ gene_idx ]] <- geneRelatedFeatures

            }
            else {
                geneRelatedFeatures <- rownames( features )[ features_S:features_E ]

                if( length( geneRelatedFeatures ) == 1 & features_S == start_endChr[ which( start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ]),"start_index"]) {
                    geneData$remarks[gene_idx] <- "A"
                    gene_features[[ gene_idx ]] <- geneRelatedFeatures
                }           

                else if( is.na(geneData$remarks[ gene_idx ]) & length(geneRelatedFeatures) > 1 ) {
                    geneData$remarks[gene_idx] <- "D"
                    start_pos_feature <- features$Start[ features_E ]
                    end_pos_feature <- NA

                    same_chr_idx <- which( start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ] )

                    if( features_S == start_endChr[ same_chr_idx, "start_index" ] ){
                        end_pos_feature <- features$Start[ features_S  ]
                    }
                    else{
                        end_pos_feature <- features$Start[ features_S -1 ]
                    }
                    geneData$genelength_features[gene_idx] <- start_pos_feature - end_pos_feature
                    gene_features[[ gene_idx ]] <- geneRelatedFeatures
                }

                else if(is.na( geneData$remarks[ gene_idx ]) & features_S == features_E) {
                    geneData$remarks[gene_idx] <- "C"
                    geneData$genelength_features[ gene_idx ] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
                    gene_features[[ gene_idx ]] <- geneRelatedFeatures
                }

                else {
                    geneData$remarks[ gene_idx ] <- "Z"
                    warning_Z <- c( "Uncategorized situations (Z) were observed" )
                }
            }
        } 
        else{ # is.na( features_S) or chr not in primary data
            ifelse( any( feature_idx_chr), geneData$remarks[ gene_idx ] <- "B", geneData$remarks[ gene_idx ] <- "X")
            if( geneData$remarks[ gene_idx ] == "B" ){
                geneRelatedFeatures <- rownames( features )[ start_endChr[ which(start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ]),"end_index"] ] # make features_E last probe on chr
            } 
            else{ 
                geneRelatedFeatures <- NA 
            }
            gene_features[[ gene_idx ]] <- geneRelatedFeatures
        }
    } # end for-loop with geneAnnotation

    ## collapse the featurenames per gene to add to geneData dataframe
    ## any() is used to return only one paste result
    geneData$featureNames <- sapply( gene_features, function( x ) { 
        ifelse( any( is.na(x) ), NA, paste( x, collapse = .SEP_CHAR ) ) 
    } )
    geneData$featureTotal <- sapply( gene_features, function(x){ length( x[!is.na(x)] )})

    ## print any warnings (once per warning)
    if( length(warning_Z) ) warning( paste( unique(warning_Z), "\n" ) )

    ## create output object
    output <- new( 'CopyNumberBreakPointGenes', 
        segmDiff    = object@segmDiff,
        callDiff    = object@callDiff,
        segments    = object@segments,
        calls       = object@calls,
        breakpoints = object@breakpoints,
        geneData    = geneData,
        featureData = object@featureData,
        featureAnnotation = object@featureAnnotation,
        geneAnnotation    = geneAnnotation,
        featuresPerGene   = gene_features
    )
    cat( "Adding gene annotation DONE\n" )
    return(output)
})

#' bpGenes
#' 
#' @description
#' Indentifies genes affected by breakpoint locations.
#' @param object An object of class \linkS4class{CopyNumberBreakPointGenes}
#' @return Returns an object of class \linkS4class{CopyNumberBreakPointGenes} with gene-breakpoint information
#' @details This step requires feature-to-gene annotations added to the input object (see \code{?addGeneAnnotation}).
#' @examples
#' data( copynumber.data.chr20 )
#' data( ens.gene.ann.hg18 )
#' 
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp )
#' bp <- addGeneAnnotation( bp, ens.gene.ann.hg18 )
#' bp <- bpGenes( bp )
#' 
#' ## options to inspect the data
#' bp
#' accessOptions( bp )
#' @aliases bpGenes
setMethod( "bpGenes", "CopyNumberBreakPointGenes",
    
    function( object ){

        ## variable setup
        fpg <- object@featuresPerGene
        bps <- object@breakpoints
        geneData <- object@geneData

        sampleNames <- colnames( bps )
        sampleCount <- ncol( bps )
        geneCount <- nrow( geneData )
        geneData$nrOfBreakLocations <- NA
        geneData$sampleCount <- NA

        geneBreakpoints <- matrix( data = 0, nrow = geneCount, ncol = sampleCount )
        colnames(geneBreakpoints) <- sampleNames

        ## ---------------
        ## start analysis
        ## ---------------
        cat( paste("Running bpGenes:", geneCount, "genes and", sampleCount, "samples\n"))
        progress <- rep( NA, geneCount )
        progress[ c( round( seq( 1, geneCount, by = (geneCount/4) ))) ] <- c("0%","25%","50%","75%")

        gene_idx_loop <- which( sapply( fpg, function(x){ length(x[!is.na(x)])}) > 0 )
        for( gene_idx in gene_idx_loop ) {

            if( !is.na( progress[gene_idx] ) ) { 
                cat( paste( progress[gene_idx], "... " ) ) 
            }

            features <- fpg[[ gene_idx ]]
            tmp_bps <- as.matrix( bps[ features, ] )

            # needs transpose at lenght 1 because of as.matrix behaviour
            if( length( features ) == 1 ) tmp_bps <- t( tmp_bps ) 

            geneBreakpoints[ gene_idx, ] <- as.vector( colSums(tmp_bps) )
            geneData$nrOfBreakLocations[ gene_idx ] <- length( which( rowSums(tmp_bps) > 0 ) )
        }
        geneData$sampleCount <- apply( geneBreakpoints, 1, function(x){ length( which( x > 0 ) )} )

        ## add a string with affected samples for each gene
        samplesListed <- apply( geneBreakpoints, 1, function(x){ (names(x)[x>0]) })
        geneData$sampleNamesWithBreakpoints <- sapply( samplesListed, function( x ) { 
            ifelse( any( is.na(x) ), NA, paste( x, collapse = .SEP_CHAR ) )
        } )

        ## add the extra slots in same object class
        object@breakpointsPerGene <- geneBreakpoints
        object@geneData  <- geneData
        geneBreaksTotal  <- sum( object@breakpointsPerGene )
        genesBrokenTotal <- length( which( rowSums( object@breakpointsPerGene ) > 0 ) )

        cat( "bpGenes DONE\n" )
        cat( "A total of ", geneBreaksTotal, " gene breaks in ", genesBrokenTotal, " genes detected\n", sep = "" )

        return(object)
    }
)

## internal function
## correction for co-variates (regression, qual filter)
## lenst = gestandaardiseerde lengte van genen
## NOTE: could add more co-variates
.glmbreak <- function( breakg, lenst = lenst, nprobes = nprobes ) {
    nullmodel <- glm( breakg ~ lenst + nprobes, family = "binomial" )
    pn <- predict( nullmodel, type = "response" )
    return( pn )
}
.cumgf <- function( probs ) {
    p = probs #p: vector of probs P(Y_i=1)
    k = length( probs )
    BS = rep(NA,k+1)
    BS[1] = 1

    for( i in 1:k ) {
        BS2 <- BS[1:i] * (1-p[i])
        BS3 <- BS[1:i] * p[i]
        BS[1:i] <- BS2
        if(i>1) BS[2:(i+1)] <- c( BS[2:i], 0 ) + BS3 else BS[ 2:(i+1) ] <- BS3
    }
    rev( cumsum( rev(BS) ) )
}
.gilbertTest <- function( pvalues=NULL, cumgfs=NULL, fdr.threshold=1 ){
    
    pvsrt <- sort( pvalues, index.return=TRUE )
    pvalssort <- pvsrt$x

    fdrgilbert <- c()
    fdrnot1 <- TRUE
    i <- 1
    fdrprev <- 0
    cumgfscut <- cumgfs
    nc <- ncol( cumgfs )
    while( fdrnot1 ) {
        pvi <- pvalssort[i]
        ## snowfall on apply?
        largsmall <- apply( cumgfscut, 1, function(ceegf) {
            ceegf <- c( ceegf, 0)
            el <- length( ceegf[ ceegf > pvi ]); 
            return( c( ceegf[ el+1 ], el + 1))
        })
        npvsmall <- i
        fdr <- min( 1, sum(largsmall[ 1,]) / npvsmall )
        if( fdr < fdrprev ) fdr <- fdrprev #enforce monotonicity
        elmax <- min( nc, max( largsmall[ 2,] ) )
        cumgfscut <- cumgfscut[ ,1:elmax ]
        i <- i + 1
        if( fdr >= fdr.threshold ) fdrnot1 <- FALSE
        fdrgilbert <- c( fdrgilbert, fdr )
        fdrprev <- fdr
    }

    ## all values above FDR threshold (fdr_threshold) are set to 1
    fdrgilbertfull <- c( fdrgilbert, rep( 1, length(pvalssort) - length(fdrgilbert)))
    fdrs <- fdrgilbertfull[ order(pvsrt$ix) ]
    return(fdrs)
}


#' bpStats
#' 
#' @description
#' Applies cohort-based statistics to identify genes and/or chromosomal locations that are recurrently affected by breakpoints.
#' @param object An object of class \linkS4class{CopyNumberBreakPointGenes}
#' @param level The level at which to operate, this can be either "gene" (correcting for gene length) or "feature" (per probe/bin)
#' @param method The FDR correction method to apply. This can be "BH" (applies Benjamini-Hochberg-type FDR correction) or "Gilbert" (for dedicated Benjamini-Hochberg-type FDR correction)
#' @param fdr.threshold The threshold for FDR correction'
#' @return Returns an object of class \linkS4class{CopyNumberBreakPointGenes} with cohort based statistics added.
#' @details The statistical method on gene-level corrects for covariates that may influence the probability to be a breakpoint gene including number of breakpoints in a profile, number of gene-associated features and gene length by gene-associated feature coverage. The statistical analysis includes multiple testing where standard Benjamini-Hochberg-type FDR correction will be performed by default. This less computational intensive method assumes a similar null-distribution for all candidate breakpoint events and satisfies for analysis on breakpoint location-level. For statistics on gene-level however, we recommend to apply the more comprehensive and powerful dedicated Benjamini-Hochberg-type FDR correction that accounts for discreteness in null-distribution (Gilbert, 2005) following correction for covariates that may influence the probability to be a breakpoint gene including number of breakpoints in a profile, number of gene-associated features and gene length by gene-associated feature coverage.
#' @references Gilbert,P.B. (2005) A modified false discovery rate multiple-comparisons procedure for discrete data, applied to human immunodeficiency virus genetics. Journal of the Royal Statistical Society Series C-Applied Statistics, 54, 143-158.
#' @examples
#' data( copynumber.data.chr20 )
#' data( ens.gene.ann.hg18 )
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp )
#' bp <- addGeneAnnotation( bp, ens.gene.ann.hg18 )
#' bp <- bpGenes( bp )
#' bp <- bpStats( bp )
#' 
#' ## options to inspect the data
#' bp
#' accessOptions( bp )
#' @aliases bpStats
setMethod( "bpStats", "CopyNumberBreakPoints",

    function( object, level="gene", method="BH", fdr.threshold=1 ) {

        allowed.methods <- c( "Gilbert", "BH" )
        allowed.levels <- c( "feature", "gene" )

        ## --------------------
        ## Check input params
        ## --------------------
        if ( ! method %in% allowed.methods ){
            stop( "Parameter \"method\" incorrect, options are: ", paste( allowed.methods, collapse=', ' ), "\n" )
        }
        if ( ! level %in% allowed.levels ){
            stop( "Parameter \"level\" incorrect, options are: ", paste( allowed.levels, collapse=', ' ), "\n" )
        }

        ## --------------------
        ## gene level analysis
        ## --------------------
        if( level == "gene" ) {

            if ( ncol( object@breakpointsPerGene ) < 5 ){
                stop( "A minimum of 5 samples is required for statistical testing...", "\n" )
            }
            if ( level == "gene" & class( object ) != "CopyNumberBreakPointGenes" ){
                stop( "Object of type \"CopyNumberBreakPointGenes\" required for analysis per gene ...", "\n" )
            }

            # setup data
            bps <- object@breakpointsPerGene
            # because of later sum() transform to 1
            bps[ bps > 0 ] <- 1
            # exclude situations A,B,X
            include <- which( !object@geneData$remark %in% c("A","B","X") )
            if ( length(include) < 1 ) stop( "No breakpoint associated genes present, please check...")

            bpt <- bps[ include, ]

            ## start analysis
            cat( paste(
                "Applying statistical test over", 
                ncol(bpt),"samples for:", 
                level, "breakpoints:", method, "test...\n"
            )) # place just before the loop #

            breaktotal <- apply( bpt, 2, sum )
            breakgene <- apply( bpt, 1, sum )

            geneAnn <- object@geneData
            len <- geneAnn$genelength_features[ include ]
            nprobes <- geneAnn$featureTotal[ include ]

            mnlen <- mean( len )
            sdlen <- sd( len )

            lenst <- sapply( len, function(x) ( x - mnlen ) / sdlen )

            nms <- apply( bpt, 2, function( breakg ){ 
                .glmbreak( breakg = breakg, lenst = lenst, nprobes = nprobes)
            }) 
            ### WARNINGS !!!
            # Warning messages:
            # 1: glm.fit: fitted probabilities numerically 0 or 1 occurred
            # 2: glm.fit: algorithm did not converge
            # 3: glm.fit: fitted probabilities numerically 0 or 1 occurred

            cumgfs <- t( apply( nms, 1, function(probs){ .cumgf( probs=probs ) }))

            pvalues <- sapply( 1:length(breakgene), function(i){ 
                cumgfs[ i, breakgene[i] + 1 ] 
            })

            if( method == "BH" ) {
                fdrs <- p.adjust( pvalues, method = "BH")
            }
            else if( method == "Gilbert" ) {
                fdrs <- .gilbertTest( pvalues, cumgfs, fdr.threshold=fdr.threshold )
            }
            else{
                stop( "Chosen method [", method, "] not supported...\n")
            }

            geneAnn$glength_stand <- NA
            geneAnn$glength_stand[ include ] <- lenst
            geneAnn$pvalue <- NA
            geneAnn$pvalue[ include ] <- pvalues
            geneAnn$FDR <- NA
            geneAnn$FDR[ include ] <- fdrs

            object@geneData <- geneAnn
            return( object )
        }

        ## --------------------
        ## feature level analysis
        ## --------------------
        else if( level == "feature" ) {
            bpt <- object@breakpoints

            cat( paste( 
                "Applying statistical test over", ncol(bpt), 
                "samples for", level, 
                "breakpoints:", method, "test...\n"
            ))

            breaktotal <- apply( bpt, 2, sum )
            breakgene <- apply( bpt, 1, sum )
            nprobe <- dim(bpt)[1]

            if( method == "BH" ) {  # dus ook andere pval calculation... 

                probuni <- rep( 1 / nprobe, nprobe )                
                probmat <- probuni %*% t( breaktotal )
                cumgfs <- .cumgf( probmat[1,] ) # bij Gilbert zou dit ook anders moeten ws ???

                pvalues <- sapply( 1:length(breakgene), function(i){ cumgfs[ breakgene[i] + 1 ] } )
                # got error; incorrect dimensions; not array but vector...

                padj <- p.adjust( pvalues, method = "BH") #equivalent to fdr-gilbert, because each probe has the same null-distribution!!!!
            }

            else if( method == "Gilbert" ) {

                len <- object@featureData$featureInterval

                mnlen <- mean( len )
                sdlen <- sd( len )

                lenst <- sapply( len, function(x) ( x - mnlen ) / sdlen )

                nprobes <- rep( 1, nprobe )
                nms <- apply( bpt, 2, function( breakg ){ 
                    .glmbreak( breakg = breakg, lenst = lenst, nprobes = nprobes) # nprobes changes into 1 (= always one feature)
                })

                cumgfs <- t( apply( nms, 1, function(probs){ .cumgf( probs=probs ) }))

                pvalues <- sapply( 1:length(breakgene), function(i){ 
                    cumgfs[i , breakgene[i] + 1 ] 
                })

                padj <- .gilbertTest( pvalues, cumgfs )
            }
            else{
                stop( "Chosen method [", method, "] not supported...\n")
            }

            currentData <- object@featureData
            currentData$nrOfBreakLocations <- breakgene
            currentData$pvalue <- pvalues
            currentData$FDR <- padj
            object@featureData <- currentData

            return( object )
        }
        else{
            stop( "Chosen level [", level, "] not supported...\n")
        }
    }
)

#' bpPlot
#' 
#' @description
#' Plots breakpoint frequencies per chromosome
#' @param object An object of class \linkS4class{CopyNumberBreakPoints} or \linkS4class{CopyNumberBreakPointGenes}
#' @param plot.chr A vector with chromosome(s) to plot. All chromosomes will be plotted when NULL is used.
#' @param plot.ylim An integer giving the max y coordinate.
#' @param fdr.threshold The FDR threshold to label recurrent breakpoint genes with their gene name
#' @param add.jitter Logical. If TRUE, function jitter will be used for the y position of gene labels
#' @return calls plot function
#' @details The plot includes breakpoint locations and breakpoint gene frequencies. Genes that are recurrently affected are labeled with their gene name.
#' @examples
#' data( copynumber.data.chr20 )
#' data( ens.gene.ann.hg18 )
#' bp <- getBreakpoints( copynumber.data.chr20 )
#' bp <- bpFilter( bp )
#' bp <- addGeneAnnotation( bp, ens.gene.ann.hg18 )
#' bp <- bpGenes( bp )
#' bp <- bpStats( bp )
#' 
#' bpPlot( bp, c(20) )
#' @aliases bpPlot
setMethod( "bpPlot", "CopyNumberBreakPoints",
    
    function( object, plot.chr=NULL, plot.ylim=15, fdr.threshold=0.1, add.jitter=FALSE ) {
        cat( paste("Plotting breakpoint frequencies ...\n") )
        
        breakpointGene <- ifelse( class(object) == "CopyNumberBreakPointGenes", 1, 0 )
        stats.gene <- ifelse( class(object) == "CopyNumberBreakPointGenes", 1, 0 )
        stats.feature <- 1

        ### input checks
        if ( class(object) == "CopyNumberBreakPointGenes" ){
        	
        	if( !exists( "sampleNamesWithBreakpoints", where=object@geneData ) ) { 
            	breakpointGene <- 0
            	cat("Breakpoint gene data is not available yet.\n")
        	}
        	if( !exists( "FDR", where=object@geneData ) ) { 
            	stats.gene <- 0
            	cat("Breakpoint gene statistics are not available yet.\n")
        	}	
        }

        if( !exists( "FDR", where=object@featureData ) ) { 
            stats.feature <- 0
        }
        
        ## set general parameters
        nstudy = length( sampleNames(object) )
        ylim = ifelse( plot.ylim < 10, 10, plot.ylim )
        color.gene <- "darkred"
        color.feature <- "black"
        color.feature.chr <- color.feature
        
        ## chromosomes to plot
        chromosomesPresent <- unique( object@featureAnnotation$Chromosome )
        chromosomesToPlot <- chromosomesPresent
        if( !is.null( plot.chr) ){
            for (chr in plot.chr){
                if( ! chr %in% chromosomesPresent ){
                    stop( "Chosen chromosome [", chr, "] not present in object-data\n", sep='' )
                }
            }
            chromosomesToPlot <- plot.chr
        }

        ### check whether chr is in data...
        for( chr in chromosomesToPlot ) {
            cat( "Plotting Chromosome: ", chr, "\n", sep="")
                    
            # subset data:
            chr.feature = which( featureChromosomes(object) == chr )
            featureBreakPerc <- apply( object@breakpoints, 1, function(x) { sum(x)/length(x)*100 } ) [chr.feature]
            recurrent.gene <- NULL
            start.feature <- object@featureAnnotation$Start[chr.feature]

            if( breakpointGene == 1 ) {
                chr.gene = which(geneChromosomes(object)==chr)
                geneBreakPerc <- apply( object@breakpointsPerGene, 1, function(x) { x[x>1]<-1 ; sum(x)/length(x)*100 } ) [chr.gene]
                start.gene <- object@geneAnnotation$Start[chr.gene]
                end.gene <- object@geneAnnotation$End[chr.gene]
                name.gene <- object@geneAnnotation$Gene[chr.gene]
            } 
            
            if( stats.gene == 1 ) {
                recurrent.gene <- which( object@geneData$FDR[chr.gene] < fdr.threshold )
            }
            
            if( stats.feature == 1 ) {
                # color recurrent breakpoints (feature level); this may be a bit overdone
                recurrent.feature <- which( object@featureData$FDR[chr.feature] < fdr.threshold )
                color.feature.chr <- rep( color.feature, length(featureBreakPerc) )
                color.feature.chr[recurrent.feature] <- "darkblue"
            }
            
            # xlim = c(startPlot,endPlot) extra zoomin function ?
            
            plot( start.feature, featureBreakPerc, ylim=c(0, ylim+1), axes=TRUE, type="h", xlab="chromosomal position (Mb)", ylab="breakpoint frequencies (%)", main=paste("BreakPoint frequencyPlot\nchromosome", chr), xaxt="n", yaxt="n",cex.lab=1.3, col= color.feature.chr)
            mtext( paste( "recurrent breakpoint genes are labeled with gene name (FDR<", fdr.threshold, ")", sep=""), side=3, col=color.gene )
            axis( 1, at=as.vector(axis(1, labels=FALSE)), labels=(as.vector(axis(1, labels=FALSE)))/10^6 ) # x-axis in MBs
            axis( 2, at=0:(ylim+1), labels=c(0:(ylim), paste(">", ylim)), las=2 )
            abline( h=c(0,ylim), lty=c(1,5) )
            
            if( breakpointGene == 1 ) {
                above.ylim <- which(geneBreakPerc > ylim)
                
                if( any(above.ylim) ) {
                    frame.plot = segments( start.gene[ -above.ylim ], geneBreakPerc[ -above.ylim ], start.gene[ -above.ylim ], geneBreakPerc[ -above.ylim ], lwd = 4, col = color.gene) # gene(levels)
                    if( stats.gene == 1 & any(recurrent.gene) ) {
                        recurrent.gene.beneath.ylim <- recurrent.gene[ which( !(recurrent.gene %in% above.ylim) ) ]
                        yvals <- geneBreakPerc[ recurrent.gene.beneath.ylim ]
                        if ( add.jitter == TRUE) yvals <- jitter( yvals, factor=1.2,amount=0.2 )
                        text( x=end.gene[ recurrent.gene.beneath.ylim ], y=yvals, name.gene[ recurrent.gene.beneath.ylim ], cex=0.7, pos=2, col=color.gene, font=4 )
                    }
                } 
                else {
                    frame.plot = segments( start.gene, geneBreakPerc, start.gene, geneBreakPerc,lwd=4,col= color.gene) # gene(levels)
                    if( stats.gene == 1 & any(recurrent.gene) ) {
                        yvals <- geneBreakPerc[ recurrent.gene ]
                        if ( add.jitter == TRUE) yvals <- jitter( yvals, factor=1.2,amount=0.2 )
                        text( x=end.gene[ recurrent.gene ], y=yvals, name.gene[ recurrent.gene ], cex=0.7, pos=2, col=color.gene, font=4 )
                    }
                }
                
#                if( stats.gene == 1 & any(recurrent.gene) ) {
#                    yvals <- geneBreakPerc[ recurrent.gene ]
#                    if ( add.jitter == TRUE) yvals <- jitter( yvals, factor=1.2, amount=0.2 )
#                    text( x=end.gene[ recurrent.gene ], y=yvals, name.gene[ recurrent.gene], cex=0.7, pos=2, col=color.gene, font=4 )
#                }

                if( any (above.ylim) ){
                    non.recurrent.gene.above.ylim <- above.ylim[ which( !(above.ylim %in% recurrent.gene) ) ]
                    sign <- rep("", length(geneBreakPerc) )
                    sign[ non.recurrent.gene.above.ylim ] <- "*"
                    frame.plot = segments( start.gene[ above.ylim ], rep( ylim+1 ), end.gene[ above.ylim ], rep( ylim+1 ), lwd = 3, col = color.gene)
                    text( x=end.gene[ above.ylim ], y=rep(ylim+1), paste( name.gene[ above.ylim ]," (", round(geneBreakPerc[ above.ylim ], 1), ")", sign[ above.ylim ], sep=""), cex=0.7, pos=4, col=color.gene, font=4)
                  
                    if ( any(non.recurrent.gene.above.ylim) ) { 
                        mtext("* not significant", side=3, adj=1, cex=0.8) 
                    }
                }
            }
        }
    }
)
