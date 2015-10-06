### R code from vignette source 'GeneBreak.Rnw'

###################################################
### code chunk number 1: loadingPackage
###################################################
library(GeneBreak)


###################################################
### code chunk number 2: settingOptions
###################################################
options("GeneBreak::verbose"=NA)
options(width=75)


###################################################
### code chunk number 3: loadingCopynumberData
###################################################
library(CGHcall)
data( "copynumber.data.chr20" )


###################################################
### code chunk number 4: displayCopynumberData
###################################################
copynumber.data.chr20


###################################################
### code chunk number 5: getBreakpoints
###################################################
breakpoints <- getBreakpoints( data = copynumber.data.chr20 )
breakpoints


###################################################
### code chunk number 6: getBreakpointsAlternative (eval = FALSE)
###################################################
## library(CGHcall)
## cgh <- copynumber.data.chr20
## 
## segmented <- data.frame( Chromosome=chromosomes(cgh), Start=bpstart(cgh),
##  End=bpend(cgh), FeatureName=rownames(cgh), segmented(cgh))
## called <- data.frame( Chromosome=chromosomes(cgh), Start=bpstart(cgh),
##  End=bpend(cgh), FeatureName=rownames(cgh), calls(cgh))
## 
## breakpoints <- getBreakpoints( data = segmented, data2 = called )


###################################################
### code chunk number 7: bpFilter
###################################################
breakpointsFiltered <- bpFilter( breakpoints, filter = "CNA-ass" )
breakpointsFiltered


###################################################
### code chunk number 8: loadingAnnotation
###################################################
data( "ens.gene.ann.hg18" )


###################################################
### code chunk number 9: showHeadGeneAnnotation
###################################################
head( ens.gene.ann.hg18 )


###################################################
### code chunk number 10: addGeneAnnotation
###################################################
breakpointsAnnotated <- addGeneAnnotation( breakpointsFiltered, ens.gene.ann.hg18 )


###################################################
### code chunk number 11: showFeatures_PCMTD2
###################################################
featuresPerGene ( breakpointsAnnotated , geneName = "PCMTD2" )


###################################################
### code chunk number 12: showGeneAssociatedFeatureInfo
###################################################
geneFeatures <- geneInfo( breakpointsAnnotated )
head( geneFeatures[ , 
  c("Gene", "Chromosome", "Start", "End", "featureTotal", 
    "featureNames", "remarks") ] )


###################################################
### code chunk number 13: bpGenes
###################################################
breakpointGenes <- bpGenes( breakpointsAnnotated )


###################################################
### code chunk number 14: headBreakpointGenes
###################################################
result_BreakpointGenes <- geneInfo ( breakpointGenes )
head( result_BreakpointGenes[ which ( result_BreakpointGenes$sampleCount > 0 ) , 
  c( "Gene", "Chromosome", "Start", "End", "featureTotal", "nrOfBreakLocations",
     "sampleCount", "sampleNamesWithBreakpoints") ] )


###################################################
### code chunk number 15: bpStats
###################################################
breakpointStatistics <- bpStats( breakpointGenes, 
  level = "gene", method = "Gilbert" )


###################################################
### code chunk number 16: recurrentGenes
###################################################
head( recurrentGenes( breakpointStatistics ) )


###################################################
### code chunk number 17: bpStats
###################################################
breakpointStatistics <- bpStats( 
  breakpointStatistics, level = "feature", method = "BH" )


###################################################
### code chunk number 18: showStatsObject
###################################################
breakpointStatistics


###################################################
### code chunk number 19: recurrentGenes
###################################################
head( featureInfo( breakpointStatistics ) )


###################################################
### code chunk number 20: GeneBreak.Rnw:251-252
###################################################
pdf("bpPlot.png", width=10)


###################################################
### code chunk number 21: bpPlot
###################################################
bpPlot( breakpointStatistics, fdr.threshold = 0.1 )


###################################################
### code chunk number 22: GeneBreak.Rnw:257-258
###################################################
dev.off()


###################################################
### code chunk number 23: createAnnotationExample (eval = FALSE)
###################################################
## # gene annotations obtained via Biomart. 
## # HUGO gene names (HGNC symbol), Ensembl_ID and chromosomal location
## 
## # Used (and most) recent releases:
## # HG18: release54
## # HG19: release75
## # HG38: release80 (date: 150629)
## 
## library(biomaRt)
## 
## ensembl54 = useMart( 
##   host = 'may2009.archive.ensembl.org', 
##   biomart = 'ENSEMBL_MART_ENSEMBL', 
##   dataset = "hsapiens_gene_ensembl" 
## )
## ensembl75 = useMart( 
##   host = 'feb2014.archive.ensembl.org', 
##   biomart = 'ENSEMBL_MART_ENSEMBL', 
##   dataset = "hsapiens_gene_ensembl" 
## )
## ensembl80 = useMart( 
##   "ensembl", 
##   dataset = "hsapiens_gene_ensembl" 
## )
## 
## createAnnotationFile <- function( biomartVersion ) {
##   biomart_result <- getBM( 
##     attributes =  c( 
##       "hgnc_symbol", "ensembl_gene_id", "chromosome_name",  
##       "start_position", "end_position", "band", "strand" 
##     ), 
##     mart = biomartVersion
##   )
## 
##   biomart_result[ ,3] <- as.vector( biomart_result[ ,3] )
##   idx_x <- biomart_result$chromosome_name == "X"
##   idx_y <- biomart_result$chromosome_name == "Y"
##   biomart_result$chromosome_name[ idx_x ] <- "23"
##   biomart_result$chromosome_name[ idx_y ] <- "24"
##   
##   biomart_genes <- biomart_result[ which(biomart_result[ ,1] != "" & 
##     biomart_result[ ,3] %in% c(1:24)) , ]
##   colnames(biomart_genes)[1:5] <- c("Gene","EnsID","Chromosome","Start","End")
##   
##   cat( 
##     c( "Biomart version:", biomartVersion@host, 
##           "including:", dim(biomart_genes)[1], "genes\n"
##     )
##   ) 
##   
##   return( biomart_genes )
## }
## 
## ens.gene.ann.hg18 <- createAnnotationFile( ensembl54 )
## ens.gene.ann.hg19 <- createAnnotationFile( ensembl75 )
## ens.gene.ann.hg38 <- createAnnotationFile( ensembl80 )
## 


###################################################
### code chunk number 24: sessionInfo
###################################################
sessionInfo()


