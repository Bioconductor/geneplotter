\name{alongChrom}
\alias{alongChrom}
\title{A function for plotting cumulative expression data from an
  exprset for a given chromosome. }
\description {
  Given a particular exprset object, a chromLocation object, and a
  chromosome name - will plot the cumulative sum of the data for all
  the genes from the desired chromosome for each patient in the exprset.
}
\usage{
alongChrom(eSet, chrom, specChrom, index=TRUE, cumul=TRUE, lTypes="1", colors="red",...)
}
\arguments{
  \item{eSet}{The exprset object to be used.}
  \item{chrom}{The desired chromosome.}
  \item{specChrom}{An object of type chromLocation for the species being
    represented.}
  \item{index}{Determines whether the X axis is displayed by gene index,
    or with the genes marked relative to their position on the
    chromosome.  Default is index mode.}
  \item{cumul}{Determines whether the Y axis is displaying cumulative
    sums or raw expression values.  Default is for cumulative sums.}
  \item{lTypes}{A vector of line types, which will be cycled. }
  \item{colors}{A vector of colors for the plots, which will be cycled.}
  \item{...}{Any remaining graphics commands may be passed along as per plot()}
}
\details{
  The genes on the chromosome of interst are extracted from the
  chromLocation object passed in, which are then intersected with the
  genes listed in the exprset.  These remaining genes will then have
  either their raw expression values or the cumulative sums of these
  values plotted.
}
\author{Jeff Gentry}
\examples{
   data(eset)
   data(hgu95AffyProbeData)
   lty <- c(1, 2, 3, 4, 5)
   cols <- c("red", "green", "blue", "orange", "magenta", "black")
   alongChrom(eset, "1", hgu95AffyProbeData, TRUE, TRUE, lty, cols)
   alongChrom(eset, "1", hgu95AffyProbeData, TRUE, FALSE, lty, cols)
   alongChrom(eset, "1", hgu95AffyProbeData, FALSE, TRUE, lty, cols)
   alongChrom(eset, "1", hgu95AffyProbeData, FALSE, FALSE, lty, cols)
}
\keyword{utilities}

