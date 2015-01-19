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
	
	segmData <- segmented(data)
	callData <- calls(data)
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
		row.names = featureNames(data) 
	)
	# probe distance is needed for statistics
	probeDistance <- c( 0, diff( bpstart(data) ) )
	probeDistance[startChr] <- 0
	featureAnnotation$probeDistance <- probeDistance
	
	# put everything in one object:
	output <- new( 'CopyNumberBreakPoints', 
		segDiff = segDiff,
		callDiff = callDiff,
		featureAnnotation = featureAnnotation,
		calls = calls(data),
		segments = segmented(data),
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
				stop( "parameter \"treshold needed\" when deltaSeg chosen as filter" )
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
	genes <- geneAnnotation
	genes$geneLength <- apply( genes[, c("Start","End") ], 1, diff )
	genes$situation <- NA
	genes$genelength_features <- NA
	genes$featureTotal <- NA
	genes$featureNames <- NA

	geneCount <- nrow( genes )


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
	cat( paste( "Breakpoint annotation started...", geneCount, "genes\n"))
	
	progress <- rep( NA, geneCount )
	progress[ c( round( seq( 1, length(genes[,"Gene"]), by=(length(genes[,"Gene"])/4))))] <- c("0%","25%","50%","75%")
	
	warning_Z <- NULL # will be set if unable to determine situation for one or more genes
	
	## Loop for all geneAnnotation
	#print(head(genes))

	for( gene_idx in 1:geneCount ) {

		if( !is.na(progress[gene_idx])) { 
			cat( paste( progress[ gene_idx ], "... ") ) 
		}
		
		feature_idx_chr <- which( features$Chromosome == genes$Chromosome[ gene_idx ] )
		tmp_features <- features[ feature_idx_chr, ]
		
		features_S <- NULL	# first probe after start position gene
		features_E <- NULL	# first probe after end position gene
		
		features_S <- head( tmp_features$idx[ which(tmp_features$End >= genes$Start[ gene_idx ])], 1) # for QDNAseq: take features$Start too; CAVE !!!
		features_E <- head( tmp_features$idx[ which(tmp_features$Start >= genes$End[ gene_idx ])], 1)
				
		geneRelatedFeatures <- NULL
		
		if( any(features_S) ) {
			
			if( !any(features_E) ) {				
				
				features_E <- start_endChr[ which( start_endChr$chr == genes$Chromosome[ gene_idx ] ), "end_index"] # make features_E last probe on chr
				geneRelatedFeatures <- rownames(features)[ features_S:features_E ]
				
				genes$situation[gene_idx] <- "E"
				genes$genelength_features[ gene_idx ] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
				gene_features[[ gene_idx ]] <- geneRelatedFeatures

			} else {
				geneRelatedFeatures <- rownames( features )[ features_S:features_E ]
				
				if( length( geneRelatedFeatures ) == 1 & features_S == start_endChr[ which( start_endChr$chr == genes$Chromosome[ gene_idx ]),"start_index"]) {
					genes$situation[gene_idx] <- "A"
					gene_features[[ gene_idx ]] <- geneRelatedFeatures
				}			
					
				else if( is.na(genes$situation[ gene_idx ]) & length(geneRelatedFeatures) > 1 ) {
					genes$situation[gene_idx] <- "D"
					genes$genelength_features[gene_idx] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
					gene_features[[ gene_idx ]] <- geneRelatedFeatures
				}
					
				else if(is.na( genes$situation[ gene_idx ]) & features_S == features_E) {
					genes$situation[gene_idx] <- "C"
					genes$genelength_features[ gene_idx ] <- features$Start[ features_E ] - features$Start[ features_S - 1 ]
					gene_features[[ gene_idx ]] <- geneRelatedFeatures
				}
				
				else {
					genes$situation[ gene_idx ] <- "Z"
					warning_Z <- "NB: uncategorized situations (Z) were observed"
				}
			}
		} else { # is.na( features_S) or chr not in primary data
			ifelse( any( feature_idx_chr), genes$situation[ gene_idx ] <- "B", genes$situation[ gene_idx ] <- "X")
			if( genes$situation[ gene_idx ] == "B" ) {
				geneRelatedFeatures <- rownames( features )[ start_endChr[ which(start_endChr$chr == genes$Chromosome[ gene_idx ]),"end_index"] ] # make features_E last probe on chr
			} 
			else { 
				geneRelatedFeatures <- NA 
			}
			gene_features[[gene_idx]] <- geneRelatedFeatures
		}
	} # end for-loop with geneAnnotation
	
	genes$featureTotal <- sapply( gene_features, function(x){ length( x[!is.na(x)] )})
	## for the names as.vector used, because sapply returns lists...?
	genes$featureNames <- as.vector( sapply( gene_features, function(x){ ifelse( is.na(x), NA, paste( x, collapse = "|" ))}) )
		
	cat("Breakpoint annotation DONE\n")
	
	if(any(warning_Z)) { 
		cat(paste(warning_Z,"\n")) 
	}
	
	## create output object
	output <- new( 'CopyNumberBreakPointGenes', 
		segDiff     = segDiff(data),
		callDiff    = callDiff(data),
		segments    = segments(data),
		calls       = calls(data),
		breakpoints = breakpoints(data),
		featureAnnotation = featureAnnotation(data),
		geneAnnotation    = genes, # same as input slot with extra colums
		featuresPerGene   = gene_features # list of length nrow genes
	)
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
		
		breakpointdata <- object
		Genes <- breakpointdata@geneAnnotation
		gene_probes <- breakpointdata@featuresPerGene
		#gene_probes <- annotations_gene[["gene_probes"]]
		
		breakpoints <- breakpointdata@breakpoints
		
		Genes$geneBreaks <- NA
		Genes$samplesWithGeneBreaks <- NA
		
		# genes x samples:
		gene_breakpoints <- matrix( data=0, nrow=nrow(Genes), ncol=ncol(breakpoints) )

		progress <- rep( NA,nrow(Genes) )
		progress[ c( round( seq( 1, nrow(Genes), by=(nrow(Genes)/4) ))) ] <- c("0%","25%","50%","75%")
		
		cat(paste("Breakpoints per gene:", nrow(Genes), "genes and", ncol(breakpoints), "samples\n"))
		
		gene_idx_loop <- which( sapply( gene_probes, function(x){ length(x[!is.na(x)])}) > 0 )
		for(gene_idx in gene_idx_loop) {
			# print(gene_idx)
			if( !is.na(progress[gene_idx]) ) { 
				cat(paste(progress[gene_idx],"... ")) 
			}

			tmp_bps <- NULL
			tmp_bps <- as.matrix( breakpoints[ unlist(gene_probes[gene_idx]), ] ) # as.matrix() in stead of rbind () ?
			
			gene_breakpoints[ gene_idx, ] <- as.vector( colSums(tmp_bps) )
			Genes$geneBreaks[ gene_idx ] <- length( which( rowSums(tmp_bps) > 0 ) )
		}
		Genes$samplesWithGeneBreaks <- apply( gene_breakpoints,1,function(x){ length( which( x > 0 ) )})
		
		# BPs_perGene <- list(
		# 	annotation=breakpointdata@featureAnnotation,
		# 	geneAnnotation=Genes,
		# 	breakpoints=breakpoints,
		# 	breakpointsPerGene=gene_breakpoints
		# )

		## create output object
		# output <- new( 'CopyNumberBreakPointGenes', 
		# 	segDiff     = segDiff(data),
		# 	callDiff    = callDiff(data),
		# 	segments    = segments(data),
		# 	calls       = calls(data),
		# 	breakpoints = breakpoints(data),
		# 	featureAnnotation = featureAnnotation(data),
		# 	geneAnnotation    = genes, # same as input slot with extra colums
		# 	featuresPerGene   = gene_features # list of length nrow genes
		# )
		
		## just add the extra slots in same class object
		object@breakpointsPerGene <- gene_breakpoints
		object@geneAnnotation <- Genes
		
		cat( "bpGenes DONE\n" )
		
		return(object)
		#return(output)
		#return(BPs_perGene)
	}
)



## internal function
.glmbreak <- function( breakg, lenst=lenst, nprobes=nprobes ) {
	nullmodel <- glm( breakg ~ lenst + nprobes, family="binomial" )
	pn <- predict( nullmodel, type="response" )
	return( pn )
}

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

## internal function
.pval_gene <- function( i, breakgene=breakgene, cumgfs=cumgfs ) {
	#i<-1
	obs <- breakgene[i]
	pv <- cumgfs[ i, obs + 1 ]
	return( pv )
}

## internal function
.pval_probe <- function( i, breakgene=breakgene, cumgfs=cumgfs ) {
	obs <- breakgene[i]
	pv <- cumgfs[ obs + 1 ] # slidely different
	return( pv )
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
setMethod( "bpStats", "CopyNumberBreakPointGenes",
	
	function( object, level, method ) {
		
		allowed.methods <- c( "Gilbert", "BH" )
		allowed.levels <- c( "all", "gene" )

		## Check input params
		if ( ! method %in% allowed.methods ){
			stop( "Parameter \"method\" incorrect, options are: ", paste( allowed.methods, collapse=', ' ), "\n" )
		}
		if ( ! level %in% allowed.levels ){
			stop( "Parameter \"level\" incorrect, options are: ", paste( allowed.levels, collapse=', ' ), "\n" )
		}

		## copy main object
		breakpointData <- object
				
		## gene level analysis
		if( level == "gene" ) {
			#bps <- breakpointData@breakpointsPerGene
			bps <- breakpointData@breakpointsPerGene
			
			bps[ bps > 0 ] <- 1
		
			# don't take situation A,B,X into account
			exclude <- which( breakpointData@geneAnnotation$situation %in% c("A","B","X") )
			
			bpt <- bps[ -exclude, ]

			cat( paste(
				"Applying statistical test over", 
				ncol(bpt),"samples for:", 
				level, "breakpoints:", method, "test...\n"
			)) # place just before the loop #

			breaktotal <- apply( bpt, 2, sum )
			breakgene <- apply( bpt, 1, sum )
			
			geneAnn <- breakpointData@geneAnnotation
			len <- geneAnn$genelength_features[ -exclude ]
			nprobes <- geneAnn$featureTotal[ -exclude ]
			
			mnlen <- mean( len )
			sdlen <- sd( len )
			lenst <- sapply( len, function(x) ( x - mnlen ) / sdlen )
			
			cat( "LENST made by now\n")

			nms <- apply( bpt, 2, function( breakg ){ .glmbreak( breakg=breakg, lenst=lenst, nprobes=nprobes) } ) ### WARNINGS !!!
				# Warning messages:
				# 1: glm.fit: fitted probabilities numerically 0 or 1 occurred
				# 2: glm.fit: algorithm did not converge
				# 3: glm.fit: fitted probabilities numerically 0 or 1 occurred
			# save(nms,file="nms.Rdata")
			cumgfs <- t( apply( nms, 1, function(probs){ .cumgf(probs=probs) }))
			# save(cumgfs,file="cumgfs.Rdata")       
			
			pvalues <- sapply( 1:length(breakgene), function(i){ .pval_gene( i=i, breakgene=breakgene, cumgfs=cumgfs) })
			
			if( method == "BH" ) {
				fdrs <- p.adjust( pvalues, method="BH")
			}
			else if( method == "Gilbert" ) {
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
					largsmall <- apply( cumgfscut, 1, function(ceegf) 
					    {
							ceegf <- c(ceegf,0)
					    	el <- length(ceegf[ceegf>pvi]); 
					    	return(c(ceegf[el+1],el+1))
					    }
					)
					npvsmall <- i
					fdr <- min(1,sum(largsmall[1,])/npvsmall)
					if( fdr < fdrprev ) fdr <- fdrprev #enforce monotonicity
					elmax <- min( nc, max(largsmall[2,] ))
					cumgfscut <- cumgfscut[ ,1:elmax]
					i <- i + 1
					#print(c(i,elmax,fdr))
					if( fdr >= 1 ) fdrnot1 <- F
					fdrgilbert <- c( fdrgilbert, fdr )
					fdrprev <- fdr
				}
				
				fdrgilbertfull <- c( fdrgilbert, rep( 1, length(pvalssort) - length(fdrgilbert) ))
				
				fdrs <- fdrgilbertfull[ order(pvsrt$ix) ]
				# allres <- data.frame(name=name,nbreak=breakgene,glength=len,glength_stand=lenst,nprobe=nprobes,pvalue=pvalues,fdr=fdrs)
			}
			else{
				stop( "WRONG METHOD\n")
			}
			
			geneAnn$glength_stand <- NA
			geneAnn$glength_stand[ -exclude ] <- lenst
			geneAnn$pvalue <- NA
			geneAnn$pvalue[ -exclude ] <- pvalues
			geneAnn$FDR <- NA
			geneAnn$FDR[ -exclude ] <- fdrs
			
			FDR_gene <- list( setdiff( ls(breakpointData), "geneAnnotation" ), geneAnnotation=geneAnn)
			return(FDR_gene)
		}
		
		# gene #
		
		##################
		
		# # # # # probe # # # # # Gilbert equivalent ...
		
		else if( level == "all" ) {
			bpt <- breakpointData@breakpoints
			
			cat( paste( 
				"Applying statistical test over", ncol(bpt), 
				"samples for:", level, 
				"breakpoints:", method, "test...\n"
			)) # place just before the loop #

			nprobe <- dim(bpt)[1]
			probuni <- rep( 1 / nprobe, nprobe )
			
			breaktotal <- apply( bpt, 2, sum )
			breakgene <- apply( bpt, 1, sum )
			
			probmat <- probuni %*% t(breaktotal)
			
			cumgfs <- cumgf( probmat[1,] )
			
			pvalues <- sapply( 1:length(breakgene), function(i){ .pval_probe( i=i, breakgene=breakgene, cumgfs=cumgfs) } )
			# same till here
			
			if(method=="BH") {
				padj <- p.adjust( pvalues, method="BH") #equivalent to fdr-gilbert, because each probe has the same null-distribution!!!!
			}
			
			FDR_probe <- data.frame( bptann=breakpointData@featureAnnotation, nbreak=breakgene, pval = pvalues, fdr = padj)
			return(FDR_probe)
		}
		else{
			stop( "WRONG LEVEL\n")
		}
	}
)
