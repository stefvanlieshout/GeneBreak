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

This package builds on to the Copy Number analysis workflows of [CGHcall] for cgh data and [QDNAseq] for NGS data. The objects created in those packages can serve as the input of BreakPointGenes (importing other data-sources will be added later).

The test-data used in the example contains only one chromosome, but a total of 200 samples.

[CGHcall]: http://www.bioconductor.org/packages/release/bioc/html/CGHcall.html
[QDNAseq]: http://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html

```R
library( "BreakPointGenes" )

# get better understanding of the package workflow
vignette( "BreakPointGenes")

# explore built-in data
data( package="BreakPointGenes" )

# get more information about built-in data
help( "copynumber.data.chr20" )

# load built-in dataset (CGHcall)
data( "copynumber.data.chr20" )

# load built-in gene annotation dataset
data( gene.annotation.hg19.chr20 )

# setup the breakpoint data
bp <- getBreakpoints( data = copynumber.data.chr20 )

# optionally filter the data
bp <- bpFilter( bp )

# setup the gene data 
bp <- addGeneAnnotation( bp, gene.annotation.hg19.chr20 )

# perform gene analysis
bp <- bpGenes( bp )

# get recurrent breakpoints
bp <- bpStats( bp )

# print object information
bp
# --- Object Info ---
# This is an object of class "CopyNumberBreakPointGenes"
# 3653 features by 200 samples
# A total of 985 breakpoints
# A total of 1029 gene breaks in 241 genes
# A total of 14 recurrent breakpoint genes (FDR < 0.1)
# See accessOptions(object) for how to access data in this object

# print some information of top 5 recurrently affected genes
recurrentGenes( bp )[ 1:5, ]
#  A total of 14 recurrent breakpoint genes (at FDR < 0.1)
#           Gene sampleCount featureTotal        pvalue           FDR
# 25468   PCMTD2          64            4 1.350385e-103 8.899035e-101
# 25481 C20orf69          33            3  5.522293e-44  1.819595e-41
# 14751    BFSP1           8            5  3.941447e-07  8.658045e-05
# 16066   ABHD12          10            9  5.756361e-05  9.483605e-03
# 15305 C20orf26           7           18  2.748743e-04  3.622843e-02

# plot results of one chromosome
bpPlot( bp, plot.chr=c(20) )
```

More information or help
---------------------

Contact us at...

Background
---------------------

More background...

