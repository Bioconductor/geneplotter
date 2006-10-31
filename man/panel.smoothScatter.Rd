\name{panel.smoothScatter}
\alias{panel.smoothScatter}
\title{ Use lattice graphics for smoothScatter }
\description{
 This function allows the user to place \code{smoothScatter} plots into
lattice graphics.
}
\usage{
panel.smoothScatter(x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white", brewer.pal(9, "Blues"))), nrpoints = 100, transformation = function(x) x^0.25, pch = ".", cex = 1, ...)
}
\arguments{
\item{x}{Numeric vector containing x-values or n by 2 matrix containing x
  and y values.}
\item{y}{Numeric vector containing y-values (optional). The length of \code{x}
  must be the same as that of \code{y}.}
\item{nbin}{Numeric vector of length 1 (for both directions) or
  2 (for x and y separately) containing the number of equally spaced
  grid points for the density estimation.}
\item{bandwidth}{Numeric vector: the smoothing bandwidth. If missing,
  these functions come up with a more or less useful guess. This
  parameter then gets passed on to the function
  \code{\link[KernSmooth::bkde2D]{bkde2D}}.}
\item{colramp}{Function accepting an integer \code{n} as an argument and
  returning \code{n} colors.}
\item{nrpoints}{Numeric vector of length 1 giving number of points to
  be superimposed on the density image. The first \code{nrpoints} points
  from those areas of lowest regional densities will be plotted. Adding
  points to the plot allows for the identification of outliers.
  If all points are to be plotted, choose \code{nrpoints = Inf}.}
\item{transformation}{Function that maps the density scale to the color scale.}
\item{pch}{\code{\link{points}} parameter setting}
\item{cex}{character expansion parameter setting
   (see \code{\link{par}}) passed to the \code{\link{points}} phase
   of the function}
\item{\dots}{Further arguments that are passed on to \code{\link{image}}.}
}

\details{

This is essentially a set of calls to the \code{smoothScatter} function, 
once the lattice device has been set up.
}
\value{
 The function is called for its side effects, namely the production
of the appropriate plots on a graphics device.
}
\author{ D. Sarkar }
\seealso{ \code{\link{smoothScatter}}}
\examples{
##something here
}
\keyword{ hplot }
