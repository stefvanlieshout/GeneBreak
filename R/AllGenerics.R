##setGeneric("applyFilters", function(object, residual=TRUE, blacklist=TRUE,
##    mappability=NA, bases=NA, chromosomes=c("X", "Y"))
##    standardGeneric("applyFilters"))

#setGeneric("BreakPointGenes", function(object, ...) standardGeneric("BreakPointGenes"))

## NOTES:
## any "new" method defined need a generic here too

setGeneric( "counts", function(object) standardGeneric("counts") )
setGeneric( "mCounts", function(object) standardGeneric("mCounts") )
setGeneric( "multiplyCounts", function(object, mFactor=2) standardGeneric("multiplyCounts") )
setGeneric( "addMultipliedCounts", function(object, mFactor=2) standardGeneric("addMultipliedCounts") )

setGeneric( "sampleNames", function(object) standardGeneric("sampleNames") )
setGeneric( "geneNames", function(object) standardGeneric("geneNames") )


## CNBPgene generics
setGeneric( "bpFilter", function(object, ...) standardGeneric("bpFilter") )
setGeneric( "bpSummary", function(object) standardGeneric("bpSummary") )
setGeneric( "addGeneAnnotation", function(object, geneAnnotations ) standardGeneric("addGeneAnnotation") )

setGeneric( "callDiff", function(object) standardGeneric("callDiff") )
setGeneric( "segDiff", function(object) standardGeneric("segDiff") )
setGeneric( "calls", function(object) standardGeneric("calls") )
setGeneric( "segments", function(object) standardGeneric("segments") )
setGeneric( "breakpoints", function(object) standardGeneric("breakpoints") )
setGeneric( "featureAnnotation", function(object) standardGeneric("featureAnnotation") )
setGeneric( "geneAnnotation", function(object) standardGeneric("geneAnnotation") )

# EOF
