\name{GetColor}
\alias{GetColor}
\title{ A function to get the Red-Blue color scheme used by dChip }
\description{
 A simple, vectorized function that computes a Red/Blue color for
 plotting microarray expression data.
}
\usage{
GetColor(value, GreenRed=FALSE, DisplayRange=3)
}

\arguments{
  \item{value}{The vector of expression values. }
  \item{GreenRed}{If \code{TRUE} the Green-Red colors are produced,
    otherwise Red-Blue are procduced. }
  \item{DisplayRange}{ A parameter controlling the range of
    \code{value}'s that will be plotted. }
}
\details{
  A simple mapping into RGB land provided by Cheng Li.
}
\value{
 A vector of RGB colors suitable for plotting in R.
}

\author{R. Gentleman, based on an original by C. Li. }

\examples{
  set.seed(10)
  x <- rnorm(10)
  GetColor(x)
}
\keyword{ manip }

