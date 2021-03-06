% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/GeneBreakAccessors.R
\docType{methods}
\name{featureChromosomes,CopyNumberBreakPoints-method}
\alias{featureChromosomes}
\alias{featureChromosomes,CopyNumberBreakPoints-method}
\title{Access Object featureChromosomes. This method returns a vector with feature chromosomes.}
\usage{
\S4method{featureChromosomes}{CopyNumberBreakPoints}(object)
}
\arguments{
\item{object}{An object of class \code{CopyNumberBreakPoints}}
}
\value{
a vector with feature chromosomes
}
\description{
Access Object featureChromosomes. This method returns a vector with feature chromosomes.
}
\examples{
data( copynumber.data.chr20 )
bp <- getBreakpoints( copynumber.data.chr20 )
featureChromosomes( bp )
}

