\name{cPlot}
\alias{cPlot}
\title{A plotting function for chromosomes. }
\description {
  Given a chromLocation object, will plot all the gene locations from
  that object.
}
\usage{
cPlot(plotChroms)
}
\arguments{
  \item{plotChroms}{An object of type chromLocation which contains all
    the gene information to be plotted. }
}
\details{
  This function will first use the lengths of the chromosomes, stored in
  the object to create scaling factors for the X axis.  Once the
  scaling factors are determined, the chromLocation object which is
  passed in is used to determine all the gene locations/strand
  information/etc, which is then plotted for the user.
}
\author{Jeff Gentry}
\examples{
     data(hgChromClass)
     cPlot(hgChromClass)
}
\keyword{utilities}

