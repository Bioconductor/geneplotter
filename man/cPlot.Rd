\name{cPlot}
\alias{cPlot}
\title{A plotting function for chromosomes. }
\description {
  Given a number of points for an X axis and a chromLocation object,
  will plot all the gene locations from that object.
}
\usage{
cPlot(xPoints, plotChroms)
}
\arguments{
  \item{xPoints}{The number of points to represent on the X axis, determines
    the resolution.}
  \item{plotChroms}{An object of type chromLocation which contains all
    the gene information to be plotted. }
}
\details{
  The X argument is used to scale the X axis on the plot.  Once the
  scaling factors are determined, the chromLocation object which is
  passed in is used to determine all the gene locations/strand
  information/etc, which is then plotted for the user.
}
\author{Jeff Gentry}
\examples {
     data(hgByChroms)
     data(cLengths)
     plotChroms <- buildChromClass("Human", "Chen Li's HG data",
byChroms, cLengths)
     cPlot(1000,plotChroms)
}
\keyword{utilities}

