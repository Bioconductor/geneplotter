\name{cPlot}
\alias{cPlot}
\title{A plotting function for chromosomes. }
\description {
  Given a chromLocation object, will plot all the gene locations from
  that object.
}
\usage{
cPlot(plotChroms, scale=c("max","relative")[1], ccol="lightgrey")
}
\arguments{
  \item{plotChroms}{An object of type chromLocation which contains all
    the gene information to be plotted. }
  \item{scale}{Passed on to cScale as it's scale argument.  Determines
    whether the graph is scaled on a relative or absolute basis.}
  \item{ccol}{The colour to be used for the genes.  Default is lightgrey.}
}
\details{
  This function will first use the lengths of the chromosomes, stored in
  the object to create scaling factors for the X axis.  Once the
  scaling factors are determined, the chromLocation object which is
  passed in is used to determine all the gene locations/strand
  information/etc, which is then plotted for the user.
}
\author{Jeff Gentry}
\seealso{\code{\link{cScale}}}
\examples{
     data(hgu95ProbLocs)

     if (interactive()) {
       par(ask=TRUE)
     }

     for (sc in c("max","relative"))
       cPlot(hgu95ProbLocs,sc,"lightgrey")
}
\keyword{utilities}

