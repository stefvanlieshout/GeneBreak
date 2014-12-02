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
### code chunk number 3: BreakPointGenes.Rnw:57-65 (eval = FALSE)
###################################################
## breakPoints <- getBP( cghCallObject )
## # all files ending in .bam from the current working directory
## 
## # or
## 
## breakPoints <- getBP( "file.txt", type='txt' )
## # file 'file.txt' from the current working directory
## 


###################################################
### code chunk number 4: BreakPointGenes.Rnw:73-75
###################################################
breakPoints <- "example_code"
breakPoints


###################################################
### code chunk number 5: BreakPointGenes.Rnw:81-82
###################################################
cgh <- matrix()


###################################################
### code chunk number 6: BreakPointGenes.Rnw:91-92
###################################################
sessionInfo()


