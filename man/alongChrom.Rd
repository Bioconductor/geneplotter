\name{cumChrom}
\alias{cumChrom}
\title{A function for plotting cumulative expression data from an
  exprset for a given chromosome. }
\description {
  Given a particular exprset object, a chromLocation object, and a
  chromosome name - will plot the cumulative sum of the data for all
  the genes from the desired chromosome for each patient in the exprset.
}
\usage{
cumChrom(eSet, chrom, specChrom, lTypes="1", colors="red")
}
\arguments{
  \item{eSet}{The exprset object to be used.}
  \item{chrom}{The desired chromosome.}
  \item{specChrom}{An object of type chromLocation for the species being
    represented.}
  \item{lTypes}{A vector of line types, which will be cycled. }
  \item{colors}{A vector of colors for the plots, which will be cycled.}
}
\details{
  The genes on the chromosome of interst are extracted from the
  chromLocation object passed in, which are then intersected with the
  genes listed in the exprset.  These remaining genes will then have
  their cumulative sum taken for each patient in the exprset - these
  values are then plotted using the line types and colors provided.
}
\author{Jeff Gentry}
\examples{
   data(eset)
   data(hgByChroms)
   data(cLengths)
   data(hgChromClass)
   lty <- c(1, 2, 3, 4, 5)
   cols <- c("red", "green", "blue", "orange", "magenta", "black")
   cumChrom(eset, "1", hgChromClass, lty, cols)
}
\keyword{utilities}

