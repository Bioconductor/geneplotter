\name{GetColor}
\alias{GetColor}
\alias{dChip.colors}
\alias{greenred.colors}
\title{ A function to get the Red-Blue color scheme used by dChip }
\description{
 A simple, vectorized function that computes a Red/Blue color for
 plotting microarray expression data.
}
\usage{
GetColor(value, GreenRed=FALSE, DisplayRange=3)
dChip.colors(n)
greenred.colors(n)
}

\arguments{
  \item{value}{The vector of expression values. }
  \item{GreenRed}{If \code{TRUE} the Green-Red colors are produced,
    otherwise Red-Blue are procduced. }
  \item{DisplayRange}{ A parameter controlling the range of
    \code{value}'s that will be plotted. }
  \item{n}{An integer saying how many colors to be in the palette.}
}
\details{
  \code{GetColor} is a simple mapping into RGB land provided by Cheng
  Li.
  \code{dChip.colors} provides functionality similar to that of
  \code{\link[grDevices]{topo.colors}} for the red--blue colors used for genome
  plots. \code{greenred.colors} does the same for the green-black-red gradient. 
}
\value{
 A vector of RGB colors suitable for plotting in R.
}

\author{R. Gentleman, based on an original by C. Li. }

\examples{
  set.seed(10)
  x <- rnorm(10)
  GetColor(x)
  dChip.colors(10)
}
\keyword{ manip }

