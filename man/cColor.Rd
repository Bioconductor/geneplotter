\name{cColor}
\alias{cColor}
\title{A function for marking specific genes on a cPlot. }
\description {
  Given a set of genes, will highlight them in the color desired on
  a plot which has already been created via the function cPlot().
}
\usage{
cColor(genes, color, plotChroms)
}
\arguments{
  \item{genes}{The genes that are being highlighted.}
  \item{color}{The color to highlight the genes.}
  \item{plotChroms}{An object of type chromLocation which contains all
    the gene information to be plotted. }
}
\details{
  It is important to call the function cPlot() first.  This function
  will then search for the specific locations of the genes desired,
  which are contained within the plotChroms instance of a chromLocation
  class.  It will then pass these on to the plotting routine to
  highlight the desired locations.
}
\author{Jeff Gentry}
\examples{
     data(hgByChroms)
     data(cLengths)
     data(hgChromClass)
     cPlot(hgChromClass)
     genes <- c("266_s_at", "31411_at", "610_at", "failExample")
     cColor(genes, "red", hgChromClass)
     genes2 <- c("960_g_at", "41807_at", "931_at", "39032_at")
     cColor(genes2, "blue", hgChromClass)
}
\keyword{utilities}

