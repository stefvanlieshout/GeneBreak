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
data( gene.annotation.hg19.chr20 )


###################################################
### code chunk number 5: BreakPointGenes.Rnw:98-99
###################################################
bp <- getBreakpoints( data = copynumber.data.chr20 )


###################################################
### code chunk number 6: BreakPointGenes.Rnw:104-105
###################################################
bp <- bpFilter( bp )


###################################################
### code chunk number 7: BreakPointGenes.Rnw:110-111
###################################################
bp <- addGeneAnnotation( bp, gene.annotation.hg19.chr20 )


###################################################
### code chunk number 8: BreakPointGenes.Rnw:116-117
###################################################
bp <- bpGenes( bp )


###################################################
### code chunk number 9: BreakPointGenes.Rnw:122-123
###################################################
bp <- bpStats( bp )


###################################################
### code chunk number 10: BreakPointGenes.Rnw:127-128 (eval = FALSE)
###################################################
## bp


###################################################
### code chunk number 11: BreakPointGenes.Rnw:141-142
###################################################
sessionInfo()


