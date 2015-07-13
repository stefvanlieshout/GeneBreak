### R code from vignette source 'GeneBreak.Rnw'

###################################################
### code chunk number 1: loadingPackage
###################################################
library(GeneBreak)


###################################################
### code chunk number 2: settingOptions
###################################################
options("GeneBreak::verbose"=NA)
options(width=40)


###################################################
### code chunk number 3: loadingCopynumberData
###################################################
data( "copynumber.data.chr20" )


###################################################
### code chunk number 4: displayCopynumberData
###################################################
copynumber.data.chr20


###################################################
### code chunk number 5: getBreakpoints
###################################################
breakpoints <- getBreakpoints( data = copynumber.data.chr20 )


###################################################
### code chunk number 6: loadingAnnotation
###################################################
data( "ens.gene.ann.hg18" )


###################################################
### code chunk number 7: bpFilter
###################################################
breakpointsFiltered <- bpFilter( breakpoints )


###################################################
### code chunk number 8: addGeneAnnotation
###################################################
breakpointsAnnotated <- addGeneAnnotation( breakpointsFiltered, ens.gene.ann.hg18 )


###################################################
### code chunk number 9: bpGenes
###################################################
breakpointGenes <- bpGenes( breakpointsAnnotated )


###################################################
### code chunk number 10: bpStats
###################################################
breakpointStatistics <- bpStats( breakpointGenes )


###################################################
### code chunk number 11: showStatsObject (eval = FALSE)
###################################################
## breakpointStatistics


###################################################
### code chunk number 12: recurrentGenes
###################################################
head( recurrentGenes( breakpointStatistics ) )


###################################################
### code chunk number 13: createAnnotationExample (eval = FALSE)
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
## ensembl54 = useMart( host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl" )
## ensembl75 = useMart( host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl" )
## ensembl80 = useMart( "ensembl", dataset="hsapiens_gene_ensembl" )
## 
## createAnnotationFile <- function( biomartVersion ) {
##   biomart_result <- getBM(attributes =  c("hgnc_symbol", "ensembl_gene_id", "chromosome_name",  "start_position", "end_position", "band", "strand"), mart = biomartVersion)
## 
##   biomart_result[ ,3] <- as.vector( biomart_result[ ,3] )
##   biomart_result$chromosome_name[ biomart_result$chromosome_name=="X" ] <- "23"
##   biomart_result$chromosome_name[ biomart_result$chromosome_name=="Y" ] <- "24"
##   
##   biomart_genes <-biomart_result[ which(biomart_result[ ,1]!="" & biomart_result[ ,3] %in% c(1:24)) , ]
##   colnames(biomart_genes)[1:5]<-c("Gene","EnsID","Chromosome","Start","End")
##   
##   cat( c("Biomart version:", biomartVersion@host, "including:", dim(biomart_genes)[1], "genes\n") ) 
##   return( biomart_genes )
## }
## 
## ens.gene.ann.hg18 <- createAnnotationFile( ensembl54 )
## ens.gene.ann.hg19 <- createAnnotationFile( ensembl75 )
## ens.gene.ann.hg38 <- createAnnotationFile( ensembl80 )
## 


###################################################
### code chunk number 14: sessionInfo
###################################################
sessionInfo()


