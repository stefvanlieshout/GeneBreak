BreakPointGenes R package
====================

Breakpoints of copy number aberrations (CNA) indicate underlying DNA breaks and thereby regions involved in structural variants. The availability of large copy number data cohorts enables to identify genes recurrently affected by structural variants. This package allows you to systematically identify recurrent CNA-associated breakpoint genes. Breakpoint detection and a tailored annotation approach for breakpoint-to-gene mapping are implemented, which takes the gene position in relation to breakpoint interval into account. Furthermore, statistics with multiple testing correction is incorporated to obtain recurrent breakpoint events. The method can be applied both on copy number data obtained from microarrays, such as array Comparative Genomic Hybridization (CGH), and on whole genome sequencing.

Installation
---------------------

Requires package "devtools":

```R
devtools::install_github( "stefvanlieshout/BreakPointGenes" )
```

Sample workflow
---------------------

This package builds on to the Copy Number analysis workflows of CGHcall (for cgh data) and QDNAseq (for NGS data). The objects created in those packages can serve as the input of BreakPointGenes (importing other data-sources will be added later).

The LGG150 data used in the example contains only one sample, which makes no sense in the context of searching for recurrently affected genes. We will add a example dataset with more samples later.

```R
library( "BreakPointGenes" )

# using build in data (CGHcall)
data( "LGG150.data" )
# using build in gene annotation data
data( gene.annotation.hg19 )
# setup the breakpoint data
bp <- getBreakpoints( data = LGG150.data )
bp <- bpFilter( object = bp, filter = "deltaSeg", threshold = 0.2 )
bp <- addGeneAnnotation( object = bp, geneAnn_max )
# perform gene analysis
bpGenes <- bpGenes( bp )
# get recurrent breakpoints
bpStats <- bpStats( bpGenes )
```

More information or help
---------------------

Contact us at...

Background
---------------------

More background...

