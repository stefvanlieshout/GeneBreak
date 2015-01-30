### R code from vignette source 'BreakPointGenes.Rnw'

###################################################
### code chunk number 1: BreakPointGenes.Rnw:22-23
###################################################
library(BreakPointGenes)


###################################################
### code chunk number 2: BreakPointGenes.Rnw:26-28
###################################################
options("BreakPointGenes::verbose"=NA)
options(width=40)


###################################################
### code chunk number 3: BreakPointGenes.Rnw:63-71 (eval = FALSE)
###################################################
## bp <- getBreakpoints( data = cghCallObject )
## # using the output of CGHcall or QDNAseq
## 
## # or
## 
## bp <- getBreakpoints( data = matrix(), data2 = matrix() )
## # providing segments and call values directly (NOT YET POSSIBLE)
## 


###################################################
### code chunk number 4: BreakPointGenes.Rnw:76-90
###################################################
library( "BreakPointGenes" )

# load built-in dataset (CGHcall)
data( "LGG150.data" )
# load built-in gene annotation dataset
data( gene.annotation.hg19 )
# setup the breakpoint data
bp <- getBreakpoints( data = LGG150.data )
# optionally filter the data
bp <- bpFilter( bp, filter = "deltaSeg", threshold = 0.2 )
# setup the gene data 
bp <- addGeneAnnotation( bp, gene.annotation.hg19 )
# perform gene analysis
bp <- bpGenes( bp )


###################################################
### code chunk number 5: BreakPointGenes.Rnw:99-100
###################################################
sessionInfo()


