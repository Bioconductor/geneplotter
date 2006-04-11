\name{smoothScatter}
\alias{smoothScatter}
\alias{densCols}
\title{Scatterplots with smoothed densities color representation}

\description{\code{smoothScatter} produces a smoothed color density
  representation of the scatterplot, obtained through a kernel
  density estimate. \code{densCols} produces a vector containing
  colors which encode the local densities at each point in a
  scatterplot.}

\usage{
smoothScatter(x, y,
              nbin = 128,
              bandwidth,
              colramp = colorRampPalette(c("white", brewer.pal(9, "Blues"))),
              nrpoints = 100,
              transformation = function(x) x^.25,
              xlab, ylab, ...) 

     densCols(x, y,
              nbin = 128,
              bandwidth,
              colramp = colorRampPalette(brewer.pal(9, "Blues")[-(1:3)]))
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
  parameter then gets passed on to function.
  \code{\link[KernSmooth::bkde2D]{bkde2D}}.}
\item{colramp}{Function accepting an integer \code{n} as an argument and
  returning \code{n} colors.}
\item{nrpoints}{Numeric vector of length 1 giving number of points to
  be superimposed on the density image.
  If all points are to be plotted, choose \code{nrpoints = Inf}.}
\item{transformation}{Function that maps the density scale to the color scale.}
\item{xlab}{Character. Gets passed on to \code{\link{image}}}
\item{ylab}{Character. Gets passed on to \code{\link{image}}}
\item{...}{Further arguments that are passed on to \code{\link{image}}.}
}

\details{
  These functions are convenience wrappers around
  \code{\link[KernSmooth::bkde2D]{bkde2D}}.
}

\value{\code{smoothScatter} is called for its side-effect, producing a
  plot on the current graphics device. \code{densCols} returns
  a vector of length \code{nrow(x)} that contains colors to be used
  in a subsequent scatterplot. Each color represents the local density
  around the corresponding point.  
}

\seealso{\code{\link[KernSmooth::bkde2D]{bkde2D}}}

\author{Florian Hahne <f.hahne@dkfz.de>}

\keyword{hplot}

\examples{
if(interactive()) {
  x1  <- matrix(rnorm(1e4), ncol=2)
  x2  <- matrix(rnorm(1e4, mean=3, sd=1.5), ncol=2)
  x   <- rbind(x1,x2)

  layout(matrix(1:4, ncol=2, byrow=TRUE))
  smoothScatter(x, nrpoints=0)
  smoothScatter(x)
  smoothScatter(x, nrpoints=Inf, colramp=colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")), bandwidth=40)  

  colors  <- densCols(x)
  plot(x, col=colors, pch=20)
}
}



