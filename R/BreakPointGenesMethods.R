#' CopyNumber to BreakPoints
#' @description
#' Builds the CopyNumberBreakPoints object with CGHcall data 
#' @param data A "CGHcall" object
#' @param first.rm Remove the first
#' @return Object of class \code{CopyNumberBreakPoints}.
#' @examples
#' getBreakpoints( cghCallObj )
getBreakpoints <- function( data, first.rm=TRUE ) {
	cat( "Breakpoint detection started...", ncol(data), " samples\n", sep="" )
	
	## input checks
	if ( class( data ) != 'cghCall' )
        stop( '[ERR] input data not a cghCall object...' )


	breakpoints <- NULL; segDiff <- NULL; callDiff <- NULL; startChr <- c()
	featureAnnotation <- NULL;	probeDistance <- c()
	
	segmData <- CGHbase::segmented(data)
	callData <- CGHbase::calls(data)
	featureNames <- rownames( segmData )
	
	segDiff  <- rbind( segmData[ 1, ], apply( segmData, 2, diff ) )	
	callDiff <- rbind( callData[ 1, ], apply( callData, 2, diff ) )
	rownames(segDiff) <- featureNames
	rownames(callDiff) <- featureNames
	
	# remove first chromosomal BP: make Bp value 0
	if ( first.rm == T ) {
		startChr <- c( 1, which( diff( chromosomes(data) ) == 1 ) + 1 )
		segDiff[ startChr, ] <- 0
		callDiff[ startChr, ] <- 0
	}
	
	featureAnnotation <- data.frame( 
		Chromosome = chromosomes(data), 
		Start = bpstart(data), 
		End = bpend(data), 
		row.names = featureNames 
	)

	# probe distance is needed for statistics
	probeDistance <- c( 0, diff( bpstart(data) ) )
	probeDistance[startChr] <- 0
	featureData <- data.frame( 
		probeDistance = probeDistance, 
		row.names = featureNames 
	)
	
	# put everything in one object:
	output <- new( 'CopyNumberBreakPoints', 
		segDiff = segDiff,
		callDiff = callDiff,
		featureAnnotation = featureAnnotation,
		featureData = featureData,
		calls = callData,
		segments = segmData,
		breakpoints = ifelse( segDiff != 0,1,0 )
	)
	return(output)
}

createGeneAnnotations <- function( file, geneCount=0 ){
	
	#if ( )
	tmp <- read.delim( file )

	## file speficic altering
	tmp[,3] <- as.vector( tmp[,3] )
	tmp$chromosome_name[ tmp$chromosome_name == "X"] <- "23"
	tmp <- tmp[ which( tmp[,1] != "" & tmp[,3] %in% c(1:23)),]
	colnames( tmp ) <- c( "Gene", "EnsID", "Chromosome", "Start", "End" )
	tmp$Chromosome <- as.integer( tmp$Chromosome )

	## make selection for speedup
	if ( geneCount > 0 ){
		tmp <- tmp[ 1:geneCount, ]
	}
	#AnnotatedDataFrame(data=tmp, varMetadata=metaData)
	#tmp <- AnnotatedDataFrame(data=tmp) # rownames not unique error...
	return( tmp )
}

setMethod( "bpFilter", "CopyNumberBreakPoints",
    function( object, filter, threshold=NULL) {
    	data <- object
        cat("Applying BP selection...\n")
	
		breakpoints_filtered <- NULL
		
		# (delta) seg value threshold:
		if( filter == "deltaSeg" ) {
			# check whether threshold exsist... (error)
			if( !is.null( threshold ) ){
				breakpoints_filtered <- ifelse( abs(data@segDiff) > threshold, 1, 0 )
			}else{
				stop( "parameter \"threshold needed\" when deltaSeg chosen as filter" )
			}
		}
		# delta call value:
		else if( filter == "deltaCall" ) {
			breakpoints_filtered <- ifelse( data@callDiff != 0, 1, 0 )
		}
		# CNA-associated breakpoints:
		else if( filter == "CNA-ass" ) {
			breakpoints_filtered <- ifelse(data@segDiff != 0, 
				ifelse( data@callDiff != 0 | data@calls != 0, 1, 0), 0)
		}
		else{
			stop( "parameter \"filter\" has incorrect value" )
		}

		# we can extent this with numerous filter steps...
		data@breakpoints <- breakpoints_filtered
        return(data)
    }
)

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
	
	data <- object
	
	## setup gene information
	geneData <- data.frame( geneLength = apply( geneAnnotation[, c("Start","End") ], 1, diff ) )
	geneData$situation <- NA
	geneData$genelength_features <- NA
	geneData$featureTotal <- NA
	geneData$featureNames <- NA

	geneCount <- nrow( geneData )
	sampCount <- length( sampleNames(data) )

	## setup probe / bin information
	features <- cbind( data@featureAnnotation, idx=1:nrow( data@featureAnnotation ) )
	#colnames(features)[5] <- "idx"

	## define index of first and last chromosomal probe
	start_endChr <- data.frame(
		chr         = features$Chromosome[ c( 1, which( diff( features$Chromosome ) == 1 ) +1 ) ],
		start_index = c( 1, which( diff( features$Chromosome ) == 1 ) +1 ),
		end_index   = c( which( diff( features$Chromosome ) == 1 ), length( features$Chromosome ) ) 
	)
	
	## make list containing geneAnnotation with related featureAnnotation:
	gene_features <- list()
		
	## loop through single geneAnnotation:
	cat( paste( "Adding of gene annotation started on ", geneCount, " genes by ", sampCount, " samples\n", sep=""))
	
	progress <- rep( NA, geneCount )
	progress[ c( round( seq( 1, geneCount, by = ( geneCount/4))))] <- c("0%","25%","50%","75%")
	
	warning_Z <- c() # will be filled if unable to determine situation for one or more genes
	
	## Loop for all geneAnnotation
	#print(head(genes))

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
		segDiff     = data@segDiff,
		callDiff    = data@callDiff,
		segments    = data@segments,
		calls       = data@calls,
		breakpoints = data@breakpoints,
		featureAnnotation = data@featureAnnotation,
		featureData = data@featureData,
		geneAnnotation    = geneAnnotation, # same as input slot
		geneData = geneData,
		featuresPerGene   = gene_features # list of length nrow genes
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
		breakpointdata <- object
		gene_probes <- breakpointdata@featuresPerGene
		breakpoints <- breakpointdata@breakpoints
		geneData <- object@geneData

		geneData$geneBreaks <- NA
		geneData$samplesWithGeneBreaks <- NA
		geneCount <- nrow( geneData )
		
		## setup genes x samples matrix
		gene_breakpoints <- matrix( data = 0, nrow = geneCount, ncol = ncol(breakpoints) )

		progress <- rep( NA, geneCount )
		progress[ c( round( seq( 1, geneCount, by = ( geneCount/4) ))) ] <- c("0%","25%","50%","75%")

		cat(paste("Breakpoints per gene:", geneCount, "genes and", ncol(breakpoints), "samples\n"))
		
		gene_idx_loop <- which( sapply( gene_probes, function(x){ length(x[!is.na(x)])}) > 0 )
		for( gene_idx in gene_idx_loop ) {
			
			if( !is.na(progress[gene_idx]) ) { 
				cat( paste( progress[gene_idx],"... " ) ) 
			}

			tmp_bps <- NULL
			tmp_bps <- as.matrix( breakpoints[ gene_probes[[ gene_idx ]], ] ) # as.matrix() in stead of rbind () ?
			if( length( gene_probes[[ gene_idx ]] ) == 1 ){
				tmp_bps <- t( tmp_bps )
			}
			
			gene_breakpoints[ gene_idx, ] <- as.vector( colSums(tmp_bps) )
			geneData$geneBreaks[ gene_idx ] <- length( which( rowSums(tmp_bps) > 0 ) )
		}
		geneData$samplesWithGeneBreaks <- apply( gene_breakpoints,1,function(x){ length( which( x > 0 ) )})
		
		## just add the extra slots in same object class
		object@breakpointsPerGene <- gene_breakpoints
		object@geneData <- geneData
		geneBreaksTotal <- sum( object@breakpointsPerGene )
        genesBrokenTotal <- length( which( rowSums( object@breakpointsPerGene ) > 0 ) )
		
		cat( "bpGenes DONE\n" )
		cat( "A total of ", geneBreaksTotal, " gene breaks in ", genesBrokenTotal, " genes detected\n", sep = "" )
		
		return(object)
	}
)



## internal function
## correction for co-variates (regression)
## NOTE: could add more co-variates
.glmbreak <- function( breakg, lenst = lenst, nprobes = nprobes ) {
	nullmodel <- glm( breakg ~ lenst + nprobes, family = "binomial" )
	pn <- predict( nullmodel, type = "response" )
	return( pn )
}
# lenst = gestandaardiseerde lengte genen
# regressie stap (qual filter)


## internal function
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

# ## internal function
# .pval_gene <- function( i, breakgene=breakgene, cumgfs=cumgfs ) {
# 	#i<-1
# 	#nr_idx <- length(i)
# 	#cat( "NR idx", nr_idx, "\n" )
# 	obs <- breakgene[i]
# 	pv <- cumgfs[ i, obs + 1 ]
# 	return( pv )
# }

# ## internal function
# .pval_probe <- function( i, breakgene=breakgene, cumgfs=cumgfs ) {
# 	obs <- breakgene[i]
# 	pv <- cumgfs[ , obs + 1 ]
# 	return( pv )
# }

.gilbertTest <- function( pvalues=NULL, cumgfs=NULL ){
	
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
				ceegf <- c(ceegf,0)
		    	el <- length(ceegf[ceegf>pvi]); 
		    	return(c(ceegf[el+1],el+1))
		    }
		)
		npvsmall <- i
		fdr <- min( 1, sum(largsmall[ 1,]) / npvsmall )
		if( fdr < fdrprev ) fdr <- fdrprev #enforce monotonicity
		elmax <- min( nc, max( largsmall[ 2,] ) )
		cumgfscut <- cumgfscut[ ,1:elmax ]
		i <- i + 1
		#print(c(i,elmax,fdr))
		if( fdr >= fdr_threshold ) fdrnot1 <- F
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
#' @param level The level
#' @param method The method
#' @return Object of same class of \code{object}.
#' @examples
#' bpStats( breakPointsGenes )
#' @aliases bpStats
setMethod( "bpStats", "CopyNumberBreakPoints",

	function( object, level="gene", method="BH", fdr_threshold=1 ) {
		
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
		

		## copy main object
		breakpointData <- object
		
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
			bps <- breakpointData@breakpointsPerGene
			# because of later sum() transform to 1
			bps[ bps > 0 ] <- 1
			# exclude situations A,B,X
			exclude <- which( breakpointData@geneData$situation %in% c("A","B","X") )
			bpt <- bps[ -exclude, ]

			## start analysis
			cat( paste(
				"Applying statistical test over", 
				ncol(bpt),"samples for:", 
				level, "breakpoints:", method, "test...\n"
			)) # place just before the loop #

			breaktotal <- apply( bpt, 2, sum )
			breakgene <- apply( bpt, 1, sum )
			
			geneAnn <- breakpointData@geneData
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
				fdrs <- .gilbertTest( pvalues, cumgfs )
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
			
			breakpointData@geneData <- geneAnn
			return(breakpointData)
		}

		## --------------------
		## feature level analysis
		## --------------------
		else if( level == "all" ) {
			bpt <- breakpointData@breakpoints
			
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

				pvalues <- sapply( 1:length(breakgene), function(i)
					{ 
						cumgfs[ breakgene[i] + 1 ] # got error; incorrect dimensions; not array but vector...
					}
				)
				
				padj <- p.adjust( pvalues, method="BH") #equivalent to fdr-gilbert, because each probe has the same null-distribution!!!!
			}

			else if( method == "Gilbert" ) {
			
				len <- breakpointData@featureData$probeDistance
				
				mnlen <- mean( len )
				sdlen <- sd( len )

				lenst <- sapply( len, function(x) ( x - mnlen ) / sdlen )

				nprobes <- rep(1,nprobe)
				nms <- apply( bpt, 2, function( breakg )
					{ 
						.glmbreak( breakg = breakg, lenst = lenst, nprobes = nprobes) # nprobes changes into 1 (= always one feature)
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
						cumgfs[i , breakgene[i] + 1 ] 
					}
				)
				
				padj <- .gilbertTest( pvalues, cumgfs )
			}
			else{
				stop( "Chosen method [", method, "] not supported...\n")
			}
			
			currentData <- breakpointData@featureData
			currentData$nbreak <- breakgene
			currentData$pval <- pvalues
			currentData$FDR <- padj
			breakpointData@featureData <- currentData

			return( breakpointData )
		}
		else{
			stop( "Chosen level not supported...\n")
		}
	}
)
