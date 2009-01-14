\name{cScale}
\alias{cScale}
\title{A function for mapping chromosome length to a number of points.}
\description{
  Given a number of points (generally representing the number of points
  on a plot's axis), and a vector of chromosome lengths - will generate
  a vector of the same length as the one passed in containing scaling
  factors for each chromosome.
}
\usage{
cScale(points, cLengths, method=c("max", "relative"), chrom)
}
\arguments{
  \item{points}{The number of points to scale the chromosome length to.}
  \item{cLengths}{A vector of chromosome lengths.}
  \item{method}{Determines whether to use relative or absolute scaling.
    Default is "max" (absolute).}
  \item{chrom}{Which chrom to determine the scale for}
}
\details{
  The scale factor is calculated in a manner based on the \code{method}
  argument.  If method is \code{max}, the factor is derived by dividing the
  points argument by each chromosome's length (in base pairs).  If the
  method chosen is \code{relative}, then the scale is determined by dividing
  the points argument by the maximum chromsome length, and applying that
  value to each chromosome.
}
\author{Jeff Gentry}
\seealso{\code{\link{cPlot}}}
\examples{
   ## A bit of a hack to not have a package dependency on hgu95av2
   ## but need to fiddle w/ the warn level to not fail the example anyways.
   curWarn <- options(warn=0)
   on.exit(options(warn=curWarn), add=TRUE)
   if (require("hgu95av2.db")) {
      z <- buildChromLocation("hgu95av2")

     for (sc in c("max","relative"))
       scale <- cScale(1000, chromLengths(z),sc,"Y")
   } else print("This example needs the hgu95av2 data package")
}
\keyword{utilities}
