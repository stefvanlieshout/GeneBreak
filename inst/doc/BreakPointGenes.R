### R code from vignette source 'BreakPointGenes.Rnw'

###################################################
### code chunk number 1: BreakPointGenes.Rnw:21-22
###################################################
library(BreakPointGenes)


###################################################
### code chunk number 2: BreakPointGenes.Rnw:25-27
###################################################
options("BreakPointGenes::verbose"=NA)
options(width=40)


###################################################
### code chunk number 3: BreakPointGenes.Rnw:34-35
###################################################
data( "copynumber.data.chr20" )


###################################################
### code chunk number 4: BreakPointGenes.Rnw:74-75
###################################################
data( "ens.gene.ann.hg18" )


###################################################
### code chunk number 5: BreakPointGenes.Rnw:97-98
###################################################
breakpointsFiltered <- bpFilter( breakpoints )


###################################################
### code chunk number 6: BreakPointGenes.Rnw:103-104
###################################################
breakpointsAnnotated <- addGeneAnnotation( breakpointsFiltered, ens.gene.ann.hg18 )


###################################################
### code chunk number 7: BreakPointGenes.Rnw:109-110
###################################################
breakpointGenes <- bpGenes( breakpointsAnnotated )


###################################################
### code chunk number 8: BreakPointGenes.Rnw:115-116
###################################################
breakpointStats <- bpStats( breakpointGenes )


###################################################
### code chunk number 9: BreakPointGenes.Rnw:120-121 (eval = FALSE)
###################################################
## breakpointStats


###################################################
### code chunk number 10: BreakPointGenes.Rnw:128-129
###################################################
head( recurrentGenes( breakpointStatistics ) )


###################################################
### code chunk number 11: BreakPointGenes.Rnw:138-170 (eval = FALSE)
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
### code chunk number 12: BreakPointGenes.Rnw:179-180
###################################################
sessionInfo()


