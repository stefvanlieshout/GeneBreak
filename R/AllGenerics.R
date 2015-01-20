##setGeneric("applyFilters", function(object, residual=TRUE, blacklist=TRUE,
##    mappability=NA, bases=NA, chromosomes=c("X", "Y"))
##    standardGeneric("applyFilters"))

#setGeneric("BreakPointGenes", function(object, ...) standardGeneric("BreakPointGenes"))

## NOTES:
## any "new" method defined need a generic here too

#setGeneric( "counts", function(object) standardGeneric("counts") )
#setGeneric( "mCounts", function(object) standardGeneric("mCounts") )
#setGeneric( "multiplyCounts", function(object, mFactor=2) standardGeneric("multiplyCounts") )
#setGeneric( "addMultipliedCounts", function(object, mFactor=2) standardGeneric("addMultipliedCounts") )

setGeneric( "sampleNames", function(object) standardGeneric("sampleNames") )
setGeneric( "geneNames", function(object) standardGeneric("geneNames") )


## -----------------
## CNBPgene generics
## -----------------

## set
setGeneric( "bpFilter", function(object, filter, ...) standardGeneric("bpFilter") )
setGeneric( "bpSummary", function(object) standardGeneric("bpSummary") )
setGeneric( "addGeneAnnotation", function(object, geneAnnotation ) standardGeneric("addGeneAnnotation") )
setGeneric( "bpGenes", function(object) standardGeneric("bpGenes") )
setGeneric( "bpStats", function(object, level="gene", method="BH", fdr_treshold=1 ) standardGeneric("bpStats") )

## get
setGeneric( "callDiff", function(object) standardGeneric("callDiff") )
setGeneric( "segDiff", function(object) standardGeneric("segDiff") )
setGeneric( "callData", function(object) standardGeneric("callData") )
setGeneric( "segmentData", function(object) standardGeneric("segmentData") )
setGeneric( "breakpointData", function(object) standardGeneric("breakpointData") )
setGeneric( "featureAnnotation", function(object) standardGeneric("featureAnnotation") )
setGeneric( "geneAnnotation", function(object) standardGeneric("geneAnnotation") )
setGeneric( "featuresPerGene", function(object, geneName) standardGeneric("featuresPerGene") )
setGeneric( "breakpointsPerGene", function(object) standardGeneric("breakpointsPerGene") )
setGeneric( "featureChromosomes", function(object) standardGeneric("featureChromosomes") )
setGeneric( "geneChromosomes", function(object) standardGeneric("geneChromosomes") )

# EOF
