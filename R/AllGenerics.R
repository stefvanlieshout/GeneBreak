## --------------------
## NOTES:
## any "new" method defined need a generic here too
## --------------------

## set
setGeneric( "bpFilter", function(object, filter="CNA-ass", threshold=NULL ) standardGeneric("bpFilter") )
setGeneric( "addGeneAnnotation", function(object, geneAnnotation ) standardGeneric("addGeneAnnotation") )
setGeneric( "bpGenes", function(object) standardGeneric("bpGenes") )
setGeneric( "bpStats", function(object, level="gene", method="BH", fdr.threshold=1 ) standardGeneric("bpStats") )

## get
setGeneric( "accessOptions", function(object) standardGeneric("accessOptions") )
setGeneric( "sampleNames", function(object) standardGeneric("sampleNames") )
setGeneric( "geneNames", function(object) standardGeneric("geneNames") )
setGeneric( "callDiff", function(object) standardGeneric("callDiff") )
setGeneric( "segDiff", function(object) standardGeneric("segDiff") )
setGeneric( "callData", function(object) standardGeneric("callData") )
setGeneric( "segmentData", function(object) standardGeneric("segmentData") )
setGeneric( "breakpointData", function(object) standardGeneric("breakpointData") )
setGeneric( "featuresPerGene", function(object, geneName) standardGeneric("featuresPerGene") )
setGeneric( "breakpointsPerGene", function(object) standardGeneric("breakpointsPerGene") )
setGeneric( "featureChromosomes", function(object) standardGeneric("featureChromosomes") )
setGeneric( "namesFeatures", function(object) standardGeneric("namesFeatures") )
setGeneric( "geneChromosomes", function(object) standardGeneric("geneChromosomes") )
setGeneric( "geneChromosomes<-", function(object,value) standardGeneric("geneChromosomes<-") )

setGeneric( "recurrentGenes", function(object, fdr.threshold=0.1, summarize=TRUE, order.column="FDR") standardGeneric("recurrentGenes") )
setGeneric( "geneInfo", function(object) standardGeneric("geneInfo") )
setGeneric( "featureInfo", function(object) standardGeneric("featureInfo") )

## output
setGeneric( "bpPlot", function(object, plot.chr=NULL, plot.ylim=15, fdr.threshold=0.1, add.jitter=FALSE) standardGeneric( "bpPlot" ))

# EOF
