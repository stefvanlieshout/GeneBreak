#' CopyNumber to BreakPoints
#' @description
#' Runs default settings for all workflow steps
#' @param data A "CGHcall" object
#' @return Output of bpStats() an object of class \code{CopyNumberBreakPointGenes}.
#' @examples
#' runWorkflow( cghCallObj )
runWorkflow <- function( object, geneAnnotation ) {
	startTime <- Sys.time()
	cat( "BreakPointGenes workflow started at: ", format(startTime),"\n", sep="" )
	
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

#' CopyNumber to BreakPoints
#' @description
#' Builds the CopyNumberBreakPoints object with CGHcall data 
#' @param data A "CGHcall" object
#' @param first.rm Remove the first breakpoint for each chromosome
#' @return Object of class \code{CopyNumberBreakPoints}
#' @examples
#' getBreakpoints( cghCallObj )
getBreakpoints <- function( data, first.rm=TRUE ) {
	cat( "Breakpoint detection started...", ncol(data), " samples\n", sep="" )
	
	## input checks
	if ( class( data ) != 'cghCall' )
        stop( '[ERR] input data not a cghCall object...' )

    ## setup variables
	## slots are the same for CGHcall and QDNAseq objects
	segmData <- data@assayData$segmented
	callData <- data@assayData$calls
	bpChrs   <- data@featureData@data$Chromosome
	bpStart  <- data@featureData@data$Start
	bpEnd    <- data@featureData@data$End
	featureNames <- rownames( segmData )

	breakpoints <- NULL; segmDiff <- NULL; callDiff <- NULL; startChr <- c()
	featureAnnotation <- NULL;	probeDistance <- c()
	
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
	probeDistance <- c( 0, diff( bpStart ) )
	probeDistance[ startChr ] <- 0
	featureData <- data.frame( 
		probeDistance = probeDistance, 
		row.names = featureNames 
	)
	
	## create object with all required slots
	output <- new( 'CopyNumberBreakPoints', 
		segmDiff = segmDiff,
		callDiff = callDiff,
		featureAnnotation = featureAnnotation,
		featureData = featureData,
		calls = callData,
		segments = segmData,
		breakpoints = ifelse( segmDiff != 0, 1, 0 )
	)
	return( output )
}

#' Build Gene Annotations
#' @description
#' Builds the geneAnnotation data.frame from a text file
#' @param file Path to a TAB delim file with gene annotations. Requires a header to be present with at least these fields: Gene, EnsID, Chromosome, Start, End
#' @param geneCount Only use first geneCount genes (for testing purposes)
#' @return Object of class \code{CopyNumberBreakPoints})
#' @examples
#' geneAnn <- createGeneAnnotations( "/path/to/annotation/file.tsv" )
createGeneAnnotations <- function( file, geneCount=NULL ){
	
	df <- read.delim( file )

	## file speficic altering
	df[,3] <- as.vector( df[,3] )
	df$chromosome_name[ df$chromosome_name == "X"] <- "23"
	df <- df[ which( df[,1] != "" & df[,3] %in% c(1:23)), ]
	colnames( df ) <- c( "Gene", "EnsID", "Chromosome", "Start", "End" )
	df$Chromosome <- as.integer( df$Chromosome )

	## make selection if geneCount is set
	if ( !is.null( geneCount ) ){
		df <- df[ 1:geneCount, ]
	}
	return( df )
}

#' Filter Breakpoints
#' @description
#' Removes breakpoints by a minimum <ADD_TERM!> threshold
#' @param filter Type of filter
#' @return Object of class \code{CopyNumberBreakPoints})
#' @examples
#' bp <- bpFilter( bp, filter = "deltaSeg", threshold = 0.2 )
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

		## return updated object
		object@breakpoints <- breakpoints_filtered
        return( object )
    }
)

## maybe remove this method??
setMethod( "bpSummary", "CopyNumberBreakPoints",
	function( object ) {
		if( nrow(object@breakpoints) > 0 ) {
			cat( "Detected BPs per sample:\n" )
			colSums( object@breakpoints )
		}
	}
)


#' addGeneAnnotation
#' 
#' @description
#' Add description...
#' @param object A "CopyNumberBreakPoints" object [output of getBreakpoints()]
#' @return Object of same class of \code{object}.
#' @examples
#' addGeneAnnotation( breakPointsFiltered, geneAnnotation )
#' @aliases addGeneAnnotation
setMethod( "addGeneAnnotation", "CopyNumberBreakPoints",
	function( object, geneAnnotation ) {
	
	## setup gene information
	geneData <- data.frame( geneLength = apply( geneAnnotation[, c("Start","End") ], 1, diff ) )
	geneData$situation <- NA
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
		features_S <- NULL	# first probe after start position gene
		features_E <- NULL	# first probe after end position gene
		
		features_S <- head( tmp_features$idx[ which(tmp_features$End >= geneAnnotation$Start[ gene_idx ])], 1) # for QDNAseq: take features$Start too; CAVE !!!
		features_E <- head( tmp_features$idx[ which(tmp_features$Start >= geneAnnotation$End[ gene_idx ])], 1)
				
		geneRelatedFeatures <- NULL
		
		if( any(features_S) ) {
			
			if( !any(features_E) ) {				
				
				features_E <- start_endChr[ which( start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ] ), "end_index"] # make features_E last probe on chr
				geneRelatedFeatures <- rownames(features)[ features_S:features_E ]
				
				geneData$situation[gene_idx] <- "E"
				geneData$genelength_features[ gene_idx ] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
				gene_features[[ gene_idx ]] <- geneRelatedFeatures

			}
			else {
				geneRelatedFeatures <- rownames( features )[ features_S:features_E ]
				
				if( length( geneRelatedFeatures ) == 1 & features_S == start_endChr[ which( start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ]),"start_index"]) {
					geneData$situation[gene_idx] <- "A"
					gene_features[[ gene_idx ]] <- geneRelatedFeatures
				}			
					
				else if( is.na(geneData$situation[ gene_idx ]) & length(geneRelatedFeatures) > 1 ) {
					geneData$situation[gene_idx] <- "D"
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
					
				else if(is.na( geneData$situation[ gene_idx ]) & features_S == features_E) {
					geneData$situation[gene_idx] <- "C"
					geneData$genelength_features[ gene_idx ] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
					gene_features[[ gene_idx ]] <- geneRelatedFeatures
				}
				
				else {
					geneData$situation[ gene_idx ] <- "Z"
					warning_Z <- c( "Uncategorized situations (Z) were observed" )
				}
			}
		} else { # is.na( features_S) or chr not in primary data
			ifelse( any( feature_idx_chr), geneData$situation[ gene_idx ] <- "B", geneData$situation[ gene_idx ] <- "X")
			if( geneData$situation[ gene_idx ] == "B" ) {
				geneRelatedFeatures <- rownames( features )[ start_endChr[ which(start_endChr$chr == geneAnnotation$Chromosome[ gene_idx ]),"end_index"] ] # make features_E last probe on chr
			} 
			else { 
				geneRelatedFeatures <- NA 
			}
			gene_features[[ gene_idx ]] <- geneRelatedFeatures
		}
	} # end for-loop with geneAnnotation
	
	## collapse the featurenames per gene to add to geneData dataframe
	## any() is used to return only one paste result
	geneData$featureNames <- sapply( gene_features, function( x ) { 
		ifelse( any( is.na(x) ), NA, paste( x, collapse = "|" ) ) 
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
#' Add description...
#' @param object A "CopyNumberBreakPointGenes" object
#' @param geneAnnotations data.frame with gene annotations
#' @return Object of same class of \code{object}.
#' @examples
#' bpGenes( breakPointsGenes, gene.annotations.hg19 )
#' @aliases bpGenes
setMethod( "bpGenes", "CopyNumberBreakPointGenes",
	
	function( object ){
		
		## variable setup
		fpg <- object@featuresPerGene
		bps <- object@breakpoints
		geneData <- object@geneData

		sampleCount <- ncol( bps )
		geneCount <- nrow( geneData )
		geneData$geneBreaks <- NA
		geneData$samplesWithGeneBreaks <- NA
		
		geneBreakpoints <- matrix( data = 0, nrow = geneCount, ncol = sampleCount )
		
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
			geneData$geneBreaks[ gene_idx ] <- length( which( rowSums(tmp_bps) > 0 ) )
		}
		geneData$samplesWithGeneBreaks <- apply( geneBreakpoints, 1, function(x){ length( which( x > 0 ) )} )
		
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
	
	pvsrt <- sort( pvalues, index.return=T )
	pvalssort <- pvsrt$x
	
	fdrgilbert <- c()
	fdrnot1 <- T
	i <- 1
	fdrprev <- 0
	cumgfscut <- cumgfs
	nc <- ncol( cumgfs )
	while( fdrnot1 ) {
		pvi <- pvalssort[i]
		## snowfall on apply?
		largsmall <- apply( cumgfscut, 1, function(ceegf) 
		    {
				ceegf <- c( ceegf, 0)
		    	el <- length( ceegf[ ceegf > pvi ]); 
		    	return( c( ceegf[ el+1 ], el + 1))
		    }
		)
		npvsmall <- i
		fdr <- min( 1, sum(largsmall[ 1,]) / npvsmall )
		if( fdr < fdrprev ) fdr <- fdrprev #enforce monotonicity
		elmax <- min( nc, max( largsmall[ 2,] ) )
		cumgfscut <- cumgfscut[ ,1:elmax ]
		i <- i + 1
		if( fdr >= fdr.threshold ) fdrnot1 <- F
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
#' Add description...
#' @param object A "CopyNumberBreakPointGenes" object
#' @param level The level ["all" or "gene"]
#' @param method The method ["BH" or "Gilbert"]
#' @return Object of same class of \code{object}.
#' @examples
#' bpStats( breakPointsGenes )
#' @aliases bpStats
setMethod( "bpStats", "CopyNumberBreakPoints",

	function( object, level="gene", method="BH", fdr.threshold=1 ) {
		
		allowed.methods <- c( "Gilbert", "BH" )
		allowed.levels <- c( "all", "gene" )

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

			if ( ncol( object@breakpointsPerGene ) < 2 ){
				stop( "A minimum of 2 samples is required for statistical testing...", "\n" )
			}
			if ( level == "gene" & class( object ) != "CopyNumberBreakPointGenes" ){
				stop( "Object of type \"CopyNumberBreakPointGenes\" required for analysis per gene ...", "\n" )
			}

			# setup data
			bps <- object@breakpointsPerGene
			# because of later sum() transform to 1
			bps[ bps > 0 ] <- 1
			# exclude situations A,B,X
			exclude <- which( object@geneData$situation %in% c("A","B","X") )
			bpt <- bps[ -exclude, ]

			## start analysis
			cat( paste(
				"Applying statistical test over", 
				ncol(bpt),"samples for:", 
				level, "breakpoints:", method, "test...\n"
			)) # place just before the loop #

			breaktotal <- apply( bpt, 2, sum )
			breakgene <- apply( bpt, 1, sum )
			
			geneAnn <- object@geneData
			len <- geneAnn$genelength_features[ -exclude ]
			nprobes <- geneAnn$featureTotal[ -exclude ]
			
			mnlen <- mean( len )
			sdlen <- sd( len )

			lenst <- sapply( len, function(x) ( x - mnlen ) / sdlen )

			nms <- apply( bpt, 2, function( breakg )
				{ 
					.glmbreak( breakg = breakg, lenst = lenst, nprobes = nprobes)
				} 
			) 
			### WARNINGS !!!
			# Warning messages:
			# 1: glm.fit: fitted probabilities numerically 0 or 1 occurred
			# 2: glm.fit: algorithm did not converge
			# 3: glm.fit: fitted probabilities numerically 0 or 1 occurred
			
			cumgfs <- t( apply( nms, 1, function(probs){ .cumgf( probs=probs ) }))
			
			pvalues <- sapply( 1:length(breakgene), function(i)
				{ 
					cumgfs[ i, breakgene[i] + 1 ] 
				}
			)
			
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
			geneAnn$glength_stand[ -exclude ] <- lenst
			geneAnn$pvalue <- NA
			geneAnn$pvalue[ -exclude ] <- pvalues
			geneAnn$FDR <- NA
			geneAnn$FDR[ -exclude ] <- fdrs
			
			object@geneData <- geneAnn
			return( object )
		}

		## --------------------
		## feature level analysis
		## --------------------
		else if( level == "all" ) {
			bpt <- object@breakpoints
			
			cat( paste( 
				"Applying statistical test over", ncol(bpt), 
				"samples for:", level, 
				"breakpoints:", method, "test...\n"
			))

			breaktotal <- apply( bpt, 2, sum )
			breakgene <- apply( bpt, 1, sum )
			nprobe <- dim(bpt)[1]
			
			if( method == "BH" ) {	# dus ook andere pval calculation... 
			
				probuni <- rep( 1 / nprobe, nprobe )				
				probmat <- probuni %*% t( breaktotal )
				cumgfs <- .cumgf( probmat[1,] ) # bij Gilbert zou dit ook anders moeten ws ???

				pvalues <- sapply( 1:length(breakgene), function(i){ cumgfs[ breakgene[i] + 1 ] } )
				# got error; incorrect dimensions; not array but vector...
				
				padj <- p.adjust( pvalues, method="BH") #equivalent to fdr-gilbert, because each probe has the same null-distribution!!!!
			}

			else if( method == "Gilbert" ) {
			
				len <- object@featureData$probeDistance
				
				mnlen <- mean( len )
				sdlen <- sd( len )

				lenst <- sapply( len, function(x) ( x - mnlen ) / sdlen )

				nprobes <- rep(1,nprobe)
				nms <- apply( bpt, 2, function( breakg )
					{ 
						.glmbreak( breakg = breakg, lenst = lenst, nprobes = nprobes) # nprobes changes into 1 (= always one feature)
					} 
				)
				
				cumgfs <- t( apply( nms, 1, function(probs){ .cumgf( probs=probs ) }))
				
				pvalues <- sapply( 1:length(breakgene), function(i)
					{ 
						cumgfs[i , breakgene[i] + 1 ] 
					}
				)
				
				padj <- .gilbertTest( pvalues, cumgfs )
			}
			else{
				stop( "Chosen method [", method, "] not supported...\n")
			}
			
			currentData <- object@featureData
			currentData$nbreak <- breakgene
			currentData$pval <- pvalues
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
#' Add description...
#' @param object A "CopyNumberBreakPoint" object
#' @param plot.chr The chromosome(s) to plot
#' @param plot.ylim The max y
#' @param fdr.threshold The FDR threshold
#' @return Nothing
#' @examples
#' bpPlot( bp_stats, c(1,22) )
#' @aliases bpPlot
setMethod( "bpPlot", "CopyNumberBreakPoints",
	
	function( object, plot.chr=NULL, plot.ylim=15, fdr.threshold=0.1 ) {
		cat( paste("Plotting breakpoint frequencies ...\n") )
		
		### input checks
		breakpointGene <- 1
		stats.gene <- 1
		stats.feature <- 1 

		if( !exists( "samplesWithGeneBreaks", where=object@geneData ) ) { 
			breakpointGene <- 0
			cat("Breakpoint gene data is not available yet.\n")
		}

		if( !exists( "FDR", where=object@geneData ) ) { 
			stats.gene <- 0
			cat("Breakpoint gene statistics are not available yet.\n")
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
			if( !all( plot.chr %in% chromosomesPresent ) ){
				stop( "Chosen chromosome [", chr, "] not present in object-data\n", sep='' )
			}
			chromosomesToPlot <- plot.chr
		}
		#chromosomesToPlot <- ifelse( plot.chr == "all", chromosomesPresent, plot.chr ) # possible to give vector with chromosomes
		


		### check whether chr is in data...
		for( chr in chromosomesToPlot ) {
			cat( "Plotting Chromosome: ", chr, "\n", sep="")
					
			# subset data:
			chr.feature = which( featureChromosomes(object) == chr )
			featureBreakPerc <- apply( object@breakpoints, 1, function(x) { sum(x)/length(x)*100 } ) [chr.feature]
			# featureBreakPerc <- object@featureData$samplesBreaks/nstudy*100 [chr.feature]  # deze kolom bestaat nog niet... 
			recurrent.gene <- NULL
			
			if( breakpointGene == 1 ) {
				chr.gene = which(geneChromosomes(object)==chr)
				start.feature <- object@featureAnnotation$Start[chr.feature]
				#geneBreakPerc <- apply( object@breakpointsPerGene, 1, function(x) { sum(x)/length(x)*100 } ) [chr.gene]
				geneBreakPerc <- apply( object@breakpointsPerGene, 1, function(x) { x[x>1]<-1 ; sum(x)/length(x)*100 } ) [chr.gene] # CHANGED
				# geneBreakPerc <- object@geneData$samplesWithGeneBreaks/nstudy*100 [chr.gene]
				start.gene <- object@geneAnnotation$Start[chr.gene]
				end.gene <- object@geneAnnotation$End[chr.gene]
				name.gene <- object@geneAnnotation$Gene[chr.gene]
			} 
			
			if( stats.gene == 1 ) {
				recurrent.gene <- which( object@geneData$FDR[chr.gene] < fdr.threshold )
			}
			
			if(stats.feature == 1 ) {
				# color recurrent breakpoints (feature level); this may be a bit overdone
				recurrent.feature <- which( object@featureData$FDR[chr.gene] < fdr.threshold )
				color.feature.chr <- rep( color.feature, length(featureBreakPerc) )
				color.feature.chr[recurrent.feature] <- "darkblue"
			}
			
			# xlim = c(startPlot,endPlot) extra zoomin function ?
			
			plot(start.feature, featureBreakPerc, ylim=c(0, ylim+1), axes=T, type="h", xlab="chromosomal position (Mb)", ylab="breakpoint frequencies (%)", main=paste("BP frequencyPlot\nchromosome", chr), xaxt="n", yaxt="n",cex.lab=1.3, col= color.feature.chr)
			axis( 1, at=as.vector(axis(1, labels=F)), labels=(as.vector(axis(1, labels=F)))/10^6 ) # x-axis in MBs
			axis( 2, at=0:(ylim+1), labels=c(0:(ylim), paste(">", ylim)), las=2 )
			abline( h=c(0,ylim), lty=c(1,5) )
			
			if( breakpointGene == 1 ) {
				above.ylim <- which(geneBreakPerc > ylim)

				#frame.plot=segments( start.gene[ -above.ylim ], geneBreakPerc[ -above.ylim ], start.gene[ -above.ylim ], geneBreakPerc[ -above.ylim ],lwd=4,col= color.gene) # gene(levels)
				if( any(above.ylim) ) { # CHANGED
					frame.plot = segments( start.gene[ -above.ylim ], geneBreakPerc[ -above.ylim ], start.gene[ -above.ylim ], geneBreakPerc[ -above.ylim ], lwd = 4, col = color.gene) # gene(levels)
				} else {
					frame.plot = segments( start.gene, geneBreakPerc, start.gene, geneBreakPerc,lwd=4,col= color.gene) # gene(levels)
				}

				
				
				if( stats.gene == 1 & any(recurrent.gene) ) {
					text(x= end.gene[ recurrent.gene ], y=jitter( geneBreakPerc[ recurrent.gene ], factor=1.2,amount=0.2 ), name.gene[ recurrent.gene], cex=0.7, pos=2, col= color.gene, font=4)
				}


					
				else if( any (above.ylim) ){
					frame.plot = segments( start.gene[ above.ylim ], rep( ylim+1 ), end.gene[ above.ylim ], rep( ylim+1 ), lwd = 3, col = color.gene) # gene(levels above ylim)
					text( x=end.gene[ above.ylim ], y=rep(ylim+1), paste( name.gene[ above.ylim ]," (", round(geneBreakPerc[ above.ylim ], 1), ")", sep=""), cex=0.4, pos=4, col=color.gene, font=4)
				}
			}
		}
	}
)