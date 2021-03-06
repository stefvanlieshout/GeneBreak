% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/GeneBreakAccessors.R
\docType{methods}
\name{accessOptions,CopyNumberBreakPoints-method}
\alias{accessOptions}
\alias{accessOptions,CopyNumberBreakPoints-method}
\title{Access Object Data. This method lists possible functions to access the data of the object.}
\usage{
\S4method{accessOptions}{CopyNumberBreakPoints}(object)
}
\arguments{
\item{object}{An object of class \code{CopyNumberBreakPoints} or \code{CopyNumberBreakPointGenes}}
}
\value{
prints text to screen
}
\description{
Access Object Data. This method lists possible functions to access the data of the object.
}
\examples{
data( copynumber.data.chr20 )
bp <- getBreakpoints( copynumber.data.chr20 )

accessOptions( bp )
}

