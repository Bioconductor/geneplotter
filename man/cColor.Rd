\name{cColor}
\alias{cColor}
\title{A function for marking specific genes on a cPlot. }
\description {
  Given a set of genes, will highlight them in the color desired on
  a plot which has already been created via the function cPlot().
}
\usage{
cColor(genes, color, xPoints, plotChroms)
}
\arguments{
  \item{genes}{The genes that are being highlighted.}
  \item{color}{The color to highlight the genes.}
  \item{xPoints}{The number of points to represent on the X axis, determines
    the resolution.}
  \item{plotChroms}{An object of type chromLocation which contains all
    the gene information to be plotted. }
}
\details{
  It is important to call the function cPlot() first, and to use the
  same xPoints parameter in each - otherwise the scaling will be off
  for the X axis.  This function will search for the specific locations
  of the genes desired, which are contained within the plotChroms
  instance of a chromLocation class.  It will then pass these on to the
  plotting routine to highlight the desired locations.
}
\author{Jeff Gentry}
\examples{
     data(hgByChroms)
     data(cLengths)
     plotChroms <- buildChromClass("Human", "Cheng Li's HG data",
byChroms, cLengths)
     cPlot(1000,plotChroms)
     genes <- c("266_s_at", "31411_at", "610_at", "failExample")
     cColor(genes, "red", 1000, plotChroms)
     genes2 <- c("960_g_at", "41807_at", "931_at", "39032_at")
     cColor(genes2, "blue", 1000, plotChroms)
}
\keyword{utilities}

