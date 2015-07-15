GeneBreak R package
====================

Breakpoints of copy number aberrations (CNA) indicate underlying DNA breaks and thereby regions involved in structural variants. The availability of large copy number data cohorts enables to identify genes recurrently affected by structural variants. This package allows you to systematically identify recurrent CNA-associated breakpoint genes. Breakpoint detection and a tailored annotation approach for breakpoint-to-gene mapping are implemented, which takes the gene position in relation to breakpoint interval into account. Furthermore, statistics with multiple testing correction is incorporated to obtain recurrent breakpoint events. The method can be applied both on copy number data obtained from microarrays, such as array Comparative Genomic Hybridization (CGH), and on whole genome sequencing.

Installation
---------------------

```R
devtools::install_github( "stefvanlieshout/GeneBreak" )
```

This requires package "devtools":

```R
install.packages("devtools")
```

Sample workflow
---------------------

Output from the Copy Number analysis workflows of [CGHcall] for cgh data and [QDNAseq] for NGS data can serve as input for the analysis of GeneBreak.

The test-data used in this example contains data from one chromosome for a total of 200 samples.

[CGHcall]: http://www.bioconductor.org/packages/release/bioc/html/CGHcall.html
[QDNAseq]: http://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html

```R
library( "GeneBreak" )

# read vignette for more explanation about the workflow
vignette( "GeneBreak")

# explore built-in data
data( package="GeneBreak" )

# get more information about built-in data
help( "copynumber.data.chr20" )

# load built-in dataset (object from CGHcall)
data( "copynumber.data.chr20" )

# load built-in gene annotation dataset (hg19 and hg38 are also availabe)
data( "ens.gene.ann.hg18" )

# setup the breakpoint data
breakpoints <- getBreakpoints( data = copynumber.data.chr20 )

# print some information about the object
breakpoints

# take a peek at the data access options
accessOptions( breakpoints )

# optionally filter the data
breakpointsFiltered <- bpFilter( breakpoints )

# add/setup the gene data 
breakpointsFiltered <- addGeneAnnotation( breakpointsFiltered, ens.gene.ann.hg18 )

# perform gene analysis
breakpointGenes <- bpGenes( breakpointsFiltered )

# get recurrent breakpoints
breakpointStatistics <- bpStats( breakpointGenes )

# print object information to screen
breakpointStatistics

# print some information of top 5 recurrently affected genes
head( recurrentGenes( breakpointStatistics ) )
#  A total of 14 recurrent breakpoint genes (at FDR < 0.1)
#           Gene sampleCount featureTotal        pvalue           FDR
# 25468   PCMTD2          64            4 1.350385e-103 8.899035e-101
# 25481 C20orf69          33            3  5.522293e-44  1.819595e-41
# 14751    BFSP1           8            5  3.941447e-07  8.658045e-05
# 16066   ABHD12          10            9  5.756361e-05  9.483605e-03
# 15305 C20orf26           7           18  2.748743e-04  3.622843e-02
# 3493      HAO1           5            5  6.528175e-04  3.961727e-02

# plot breakpoint frequencies of one chromosome
bpPlot( breakpointStatistics, plot.chr=c(20) )
```

More information or help
---------------------

Contact us at...

Background
---------------------

More background...

