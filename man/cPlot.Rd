\name{cPlot}
\alias{cPlot}
\title{A plotting function for chromosomes. }
\description {
  Given a data file containing gene information, including the
  gene name, chromosome number, start location, and strand type,
  will create a graphic depicting this information.
}
\usage{
cPlot(xPoints, fileName)
}
\arguments{
  \item{xPoints}{The number of points to represent on the X axis, determines
    the resolution.}
  \item{fileName}{The filename of the data file.  Stored in
    geneplotter/data.}
}
\details{
  The X argument is used to scale the X axis on the plot.  The other
  determining factor is the length of the chromosome - the scale is
  simply the number of points divided by the length.  This scale is then
  applied to all genes on the particular chromosome.
}
\author{Jeff Gentry}
\examples {
     cPlot(1000,"data/hg_byChroms")
}
\keyword{utilities}

