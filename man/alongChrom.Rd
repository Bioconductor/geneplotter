\name{alongChrom}
\alias{alongChrom}
\title{A function for plotting expression data from an exprset for a
  given chromosome. }
\description {
  Given a particular exprset object, a chromLocation object, and a
  chromosome name - will plot the cumulative sum of the data for all
  the genes from the desired chromosome for each patient in the
  exprset.  This function returns an environment which contains
  information that can be used by the function identifyLines(), which
  will apply an identify() on the plot.
}
\usage{
   alongChrom(eSet, chrom, specChrom,  xlim, whichGenes, 
   xloc=c("equispaced", "physical"), 
   plotFormat=c("cumulative", "local","image"), 
   scale=c("none","zscale","rankscale","rangescale","zrobustscale"),
   geneSymbols=FALSE, byStrand=FALSE, lty=1, type="S",colors="red",...)}
\arguments{
  \item{eSet}{The exprset object to be used.}
  \item{chrom}{The desired chromosome.}
  \item{specChrom}{An object of type chromLocation for the species being
    represented.}
  \item{xlim}{A pair of values - either character or integer, which will
    denote the range of genes to display (based on base pair: either directly
    in the case of integers, or using the locations of the named genes if
    character).  If not supplied, the entire chromosome is used.}
  \item{whichGenes}{If supplied, will limit the displayed genes to the ones
    provided in this vector.}
  \item{xloc}{Determines whether the X axis points (gene names) will be
    displayed according to their relative position on the chromosome
    (physical), or spaced evenly (equispaced).  Default is equispaced.}
  \item{plotFormat}{Determines the method which to plot the data.}
  \item{scale}{Determines what method of scaling will be applied to the
    data.  Default is none.}
  \item{geneSymbols}{Notes whether to use Affy IDs or Gene Symbols,
    default is Affy IDs}
  \item{byStrand}{Determines whether to show the entire plot at once, or
    a split plot by strands.  Default is a singular plot}
  \item{lty}{A vector of line types, which will be cycled. }
  \item{type}{Plot type, from par.  Defaults to "S".}
  \item{colors}{A vector of colors for the plots, which will be cycled.}
  \item{...}{Any remaining graphics commands may be passed along as per plot()}
}
\details{
  The genes on the chromosome of interest are extracted from the
  chromLocation object passed in, which are then intersected with the
  genes listed in the exprset.  These remaining genes will then be
  plotted according to the 'plotFormat' argument.  If "image" is
  specified, an image plot is created showing the expression levels of
  the samples by gene, using a colour map to denote the levels.  If
  "cumulative" is chosen, the cumulative expression level is plotted
  against the genes for each sample.  Likewise, if "local" is used, the
  raw data is plotted for each sample against the genes.
}
\author{Jeff Gentry}
\examples{
   data(eset)
   data(hgu95AProbLocs)
   lty <- c(1, 2, 3, 4, 5)
   cols <- c("red", "green", "blue", "orange", "magenta", "black")
   cols <- cols[eset$cov3]
   if (interactive()) {
      par(ask=TRUE)
   }

   ## Here we're using xlim to denote a physical region to display
   xlim <- c(87511280,127717880)
   for (xl in c("equispaced", "physical"))
     for (sc in c("none","rangescale"))
     {
        alongChrom(eset, "1", hgu95AProbLocs, xlim=xlim, xloc=xl,
           plotFormat="cumulative", scale=sc,lty=lty, col=cols)
     }

  ## Here we're looking for specific genes
  which <- c("31540_at","31583_at", "31508_at", "31529_at", "31439_f_at",
            "31729_at")
  ## Gene "31529_at" does not exist in the current set of genes,
  ## here it demonstrates ohw genes not available are dropped.
   for (xl in c("equispaced", "physical"))
     for (sc in c("none","rangescale"))
     {
        alongChrom(eset, "1", hgu95AProbLocs, which=which, xloc=xl,
           plotFormat="cumulative", scale=sc,lty=lty, col=cols)
     }

  ## Do an image plot
  for (bs in c(TRUE,FALSE))
     alongChrom(eset, "1", hgu95AProbLocs, xlim=xlim, plotFormat="image",
               scale="zscale", byStrand=bs)
}
\keyword{utilities}

