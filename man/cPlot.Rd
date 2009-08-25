\name{cPlot}
\alias{cPlot}
\title{A plotting function for chromosomes.}
\description{
  Given a chromLocation object, will plot all the gene locations from
  that object.
}
\usage{
cPlot(plotChroms, useChroms=chromNames(plotChroms),
      scale=c("relative","max"), fg="white", bg="lightgrey", 
      glen=0.4, xlab="", ylab="Chromosome",  main = organism(plotChroms))
}
\arguments{
  \item{plotChroms}{An object of type chromLocation which contains all
    the gene information to be plotted.}
  \item{useChroms}{A vector of chromosome names to be used in the plot.
  Default is to use all the chromosomes from the plotChroms object.}
  \item{scale}{Passed on to cScale as it's scale argument.  Determines
    whether the graph is scaled on a relative or absolute basis.}
  \item{fg}{The colour to be used for the genes.  Default is white.}
  \item{bg}{The colour to be used for the background of the plot.
    Defaults to lightgrey.}
  \item{glen}{A scaling factor applied to the plotted length of each
    gene.  Defaults to 0.4 - it is recommended that this not be set
    larger then 0.5 as it will cause overlap between chromosomes.}
  \item{xlab}{A label for the x axis.}
  \item{ylab}{A label for the y axis.}
  \item{main}{A main label for the plot.}
}
\details{
  This function will first use the lengths of the chromosomes, stored in
  the object to create scaling factors for the X axis.  Once the
  scaling factors are determined, the \code{chromLocation} object which is
  passed in is used to determine all the gene locations/strand
  information/etc, which is then plotted for the user.
}
\author{Jeff Gentry}
\seealso{\code{\link{cScale}}, \code{\link{cColor}},
  \code{\link[annotate]{chromLocation-class}}}
\examples{
   ## A bit of a hack to not have a package dependency on hgu95av2
   ## but need to fiddle w/ the warn level to not fail the example anyways.

   curWarn <- options(warn=0)
   on.exit(options(curWarn), add=TRUE)
   if (require("hgu95av2.db")) {
       z <- buildChromLocation("hgu95av2")

       if (interactive()) {
           curPar <- par(ask=TRUE)
           on.exit(par(curPar), add=TRUE)
       }

       for (sc in c("max","relative"))
           cPlot(z,c("1","5","10","X","Y"),sc)
   } else print("This example can not be run without hgu95av2 data package")
}
\keyword{utilities}
