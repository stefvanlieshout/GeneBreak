#' CopyNumber to BreakPoints
#' @description
#' Builds the CopyNumberBreakPoints object with CGHcall data 
#' @param data A "CGHcall" object
#' @param first.rm Remove the first
#' @return Object of class \code{CopyNumberBreakPoints}.
#' @examples
#' getBP( cghCallObj )
getBP <- function( data, first.rm=TRUE ) {
	cat( "Breakpoint detection started...", ncol(data), " samples\n", sep="" )
	
	## input checks
	if ( class( data ) != 'cghCall' )
        stop( '[ERR] input data not a cghCall object...' )

	breakpoints <- NULL; segDiff <- NULL; callDiff <- NULL; startChr <- c()
	featureAnnotation <- NULL;	probeDistance <- c()
	
	segDiff <- rbind(segmented(data)[1,], apply( segmented(data),2,diff) )	
	callDiff <- rbind(calls(data)[1,], apply( calls(data),2,diff) )
	rownames(segDiff) <- rownames( segmented(data) )
	rownames(callDiff) <- rownames( segmented(data) )
	
	# remove first chromosomal BP: make Bp value 0
	if (first.rm==T) {
		startChr <- c(1,which(diff(chromosomes(data))==1)+1)
		segDiff[startChr,] <- 0
		callDiff[startChr,] <- 0
	}
	
	featureAnnotation <- data.frame( 
		Chromosome = chromosomes(data), 
		Start = bpstart(data), 
		End = bpend(data), 
		row.names = featureNames(data) 
	)
	# probe distance is needed for statistics
	probeDistance <- c(0,diff(bpstart(data)))
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
	
	tmp <- read.delim( file )

	## file speficic altering
	tmp[,3] <- as.vector( tmp[,3] )
	tmp$chromosome_name[ tmp$chromosome_name=="X"] <- "23"
	tmp <- tmp[ which( tmp[,1]!="" & tmp[,3]%in% c(1:23)),]
	colnames( tmp ) <- c("Gene","EnsID","Chromosome","Start","End")
	tmp$Chromosome <- as.integer( tmp$Chromosome )

	## make selection for speedup
	if ( geneCount > 0 )
	tmp <- tmp[1:geneCount,]
	return( tmp )
}

setMethod( "bpFilter", "CopyNumberBreakPoints",
    function( object, filter="none", threshold=NULL) {
    	data <- object
        cat("Applying BP selection...\n")
	
		breakpoints_filtered <- NULL
		
		# (delta) seg value threshold:
		if( filter=="deltaSeg" ) {
			# check whether threshold exsist... (error)
			if( !is.null( threshold ) ){
				breakpoints_filtered <- ifelse( abs(data@segDiff) > threshold, 1, 0 )
			}else{
				stop( "parameter \"treshold needed\" when deltaSeg chosen as filter" )
			}
		}
		# delta call value:
		else if( filter=="deltaCall" ) {
			breakpoints_filtered <- ifelse( data@callDiff != 0,1,0 )
		}
		# CNA-associated breakpoints:
		else if( filter=="CNA-ass" ) {
			breakpoints_filtered <- ifelse(data@segDiff != 0, 
				ifelse( data@callDiff !=0 | data@calls != 0, 1,0), 0)
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


#' geneToProbeAnn
#' 
#' @description
#' \code{geneToProbeAnn} returns a list with.....
#' @param object A "CopyNumberBreakPoints" object [output of getBP()]
#' @return Object of same class of \code{object}.
#' @examples
#' getBP()
#' BPsample( bp, filter="none" )
#' @aliases geneToProbeAnn
setMethod( "addGeneAnnotation", "CopyNumberBreakPoints",
	function( object, geneAnnotations ) {
	
	data <- object
	Genes <- geneAnnotations
	Genes$geneLength <- apply( Genes[, c("Start","End") ], 1, diff )
	Genes$situation <- NA
	Genes$genelength_probes <- NA
	Genes$probesTotal <- NA
	Genes$probeNames <- NA

	geneCount <- nrow(Genes)
	
	Probes <- cbind( data@featureAnnotation, idx=1:nrow(data@featureAnnotation) )

	## index of first and last chromosomal probe
	start_endChr <- data.frame(
		chr         = Probes$Chromosome[ c(1,which(diff(Probes$Chromosome)==1)+1) ],
		start_index = c( 1,which(diff(Probes$Chromosome)==1)+1 ),
		end_index   = c( which(diff(Probes$Chromosome)==1),length(Probes$Chromosome)) 
	)
	
	# make list containing genes with related probes:
	gene_probes <- list()
		
	# loop through single genes:
	cat( paste( "Breakpoint annotation started...", geneCount, "genes\n"))
	
	progress <- rep(NA,length(Genes[,"Gene"]))
	progress[c(round(seq(1,length(Genes[,"Gene"]),by=(length(Genes[,"Gene"])/4))))] <- c("0%","25%","50%","75%")
	
	warning_Z <- NULL
	
	# Loop for all genes
	for( gene_idx in 1:length(Genes[,"Gene"])) {
		# print(gene_idx)
		if(!is.na(progress[gene_idx])) { cat(paste(progress[gene_idx],"... ")) }
		
		probe_idx_chr <- NULL
		tmp_probes <- NULL
		probe_idx_chr <- which(Probes$Chromosome==Genes$Chromosome[gene_idx])
		tmp_probes <- Probes[probe_idx_chr,]
		
		probes_S <- NULL	# first probe after start position gene
		probes_E <- NULL	# first probe after end position gene
		
		# probes_S <- head(which(Probes$Chromosome==Genes$Chromosome[gene_idx] & Probes$End>=Genes$Start[gene_idx]),1) # for QDNAseq: take Probes$Start too; CAVE !!!
		# probes_E <- head(which(Probes$Chromosome==Genes$Chromosome[gene_idx] & Probes$Start>=Genes$End[gene_idx]),1)
		
		probes_S <- head(tmp_probes$idx[which(tmp_probes$End>=Genes$Start[gene_idx])],1) # for QDNAseq: take Probes$Start too; CAVE !!!
		probes_E <- head(tmp_probes$idx[which(tmp_probes$Start>=Genes$End[gene_idx])],1)
				
		geneRelatedProbes <- NULL
		
		if(any(probes_S)) {
			
			if(!any(probes_E)) {				
				probes_E<-start_endChr[which(start_endChr$chr==Genes$Chromosome[gene_idx]),"end_index"] # make probes_E last probe on chr
				geneRelatedProbes <- rownames(Probes)[probes_S:probes_E]
				
				Genes$situation[gene_idx] <- "E"
				Genes$genelength_probes[gene_idx] <- Probes$Start[probes_E]-Probes$Start[probes_S-1]
				# Genes$probeNames[gene_idx] <- paste(geneRelatedProbes,collapse="&")
				
				gene_probes[[gene_idx]] <- geneRelatedProbes
				# Genes$probesTotal[gene_idx] <- length(geneRelatedProbes)
			} else {
				# in case sit != E|B :
				geneRelatedProbes <- rownames(Probes)[probes_S:probes_E]
				
				if(length(geneRelatedProbes)==1 & probes_S==start_endChr[which(start_endChr$chr==Genes$Chromosome[gene_idx]),"start_index"]) {
					Genes$situation[gene_idx] <- "A"
					# Genes$probeNames[gene_idx] <- paste(geneRelatedProbes,collapse="&") 
					
					gene_probes[[gene_idx]] <- geneRelatedProbes
					# Genes$probesTotal[gene_idx] <- length(geneRelatedProbes)
					}			
					
				else if(is.na(Genes$situation[gene_idx]) & length(geneRelatedProbes)>1) {
					Genes$situation[gene_idx] <- "D"
					Genes$genelength_probes[gene_idx] <- Probes$Start[probes_E]-Probes$Start[probes_S-1]
						# CAVE 1st probe of a chr-arm (centromere) ; suggestion: calc FDR based on genelength_probes in stead of geneLength
					# Genes$probeNames[gene_idx] <- paste(geneRelatedProbes,collapse="&")
					
					gene_probes[[gene_idx]] <- geneRelatedProbes
					# Genes$probesTotal[gene_idx] <- length(geneRelatedProbes)
					}
					
				else if(is.na(Genes$situation[gene_idx]) & probes_S==probes_E) { # & probes_S!=start_endChr[which(start_endChr$chr==Genes$Chromosome[gene_idx]),"start_index"]) {
					Genes$situation[gene_idx] <- "C"
					Genes$genelength_probes[gene_idx] <- Probes$Start[probes_E]-Probes$Start[probes_S-1]
						# CAVE 1st probe of a chr-arm
					# Genes$probeNames[gene_idx] <- paste(geneRelatedProbes,collapse="&")
					
					gene_probes[[gene_idx]] <- geneRelatedProbes
					# Genes$probesTotal[gene_idx] <- length(geneRelatedProbes)
					}
				else {
					Genes$situation[gene_idx] <- "Z"
					warning_Z <- "NB: uncategorized situations (Z) were observed"
				}
			}
		} else {	# is.na(probes_S) or chr not in primary data
			ifelse(any(probe_idx_chr),Genes$situation[gene_idx]<-"B",Genes$situation[gene_idx]<-"X")
			if(Genes$situation[gene_idx]=="B") {
				geneRelatedProbes <- rownames(Probes)[ start_endChr[which(start_endChr$chr==Genes$Chromosome[gene_idx]),"end_index"] ] # make probes_E last probe on chr
				} 
				else { geneRelatedProbes <- NA }
			gene_probes[[gene_idx]] <- geneRelatedProbes
			# Genes$probesTotal[gene_idx] <- NA
		}
	} # end for-loop with genes
	
	Genes$probesTotal <- sapply(gene_probes,function(x){length(x[!is.na(x)])})
	Genes$probeNames <- sapply(gene_probes,function(x){paste(x,collapse="&")})
		
	# put everything in one annotation-object:
	annotations <- list( geneAnnotation=Genes, gene_probes=gene_probes ) 
	cat("Breakpoint annotation done\n")
	
	if(any(warning_Z)) { 
		cat(paste(warning_Z,"\n")) 
	}
	
	# put everything in one object:
	output <- new( 'CopyNumberBreakPointGenes', 
		segDiff = segDiff(data),
		callDiff = callDiff(data),
		featureAnnotation = featureAnnotation(data),
		calls = calls(data),
		segments = segments(data),
		breakpoints = breakpoints(data),
		geneAnnotation = Genes,
		geneProbes = gene_probes
	)
	return(output)
	#return(annotations)
})