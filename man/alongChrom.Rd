\name{alongChrom}
\alias{alongChrom}
\title{A function for plotting cumulative expression data from an
  exprset for a given chromosome. }
\description {
  Given a particular exprset object, a chromLocation object, and a
  chromosome name - will plot the cumulative sum of the data for all
  the genes from the desired chromosome for each patient in the
  exprset.  This function returns an environment which contains
  information that can be used by the function identifyLines(), which
  will apply an identify() on the plot.
}
\usage{
alongChrom(eSet, chrom, specChrom,  xloc=c("equispaced", "physical")[1],
plotFormat=c("cumulative", "local")[1],
scale=c("none","zscale","rankscale","rangescale","zrobustscale")[1], lty=1, colors="red",...)
}
\arguments{
  \item{eSet}{The exprset object to be used.}
  \item{chrom}{The desired chromosome.}
  \item{specChrom}{An object of type chromLocation for the species being
    represented.}
  \item{xloc}{Determines whether the X axis points (gene names) will be
    displayed according to their relative position on the chromosome
    (physical), or spaced evenly (equispaced).  Default is equispaced.}
  \item{plotFormat}{Determines whether the cumulative sum is to be
    generated from the expression data, or just the raw data.  Default
    is cumulative.}
  \item{scale}{Determines what method of scaling will be applied to the
    data.  Default is none.}
  \item{lty}{A vector of line types, which will be cycled. }
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
   data(hgu95ProbLocs)
   lty <- c(1, 2, 3, 4, 5)
   cols <- c("red", "green", "blue", "orange", "magenta", "black")
   cols <- cols[eset$cov3]
   if (interactive()) {
      par(ask=TRUE)
   }
   for (xl in c("equispaced", "physical"))
     for (sc in c("none","rangescale"))
     {
        alongChrom(eset, "1", hgu95ProbLocs, xloc=xl, plotFormat="cumulative",
                   scale=sc,lty, cols)
     }
}
\keyword{utilities}

