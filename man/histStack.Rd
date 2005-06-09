\name{histStack}
\alias{histStack}
\title{Stacked histogram}
\description{Stacked histogram}
\usage{histStack(x, breaks, midsFun=paste, ...) }
\arguments{
  \item{x}{A list of numeric vectors.}
  \item{breaks}{Histogram breaks, as in
    \code{\link[graphics:hist]{hist}}}
  \item{midsFun}{Function, can be used to control the formatting and
    number of the bar labels on the x-axis. See example.}
  \item{...}{Further arguments that get passed to
    \code{\link[graphics:barplot]{barplot}}}
}

\details{The function calls \code{\link[graphics:hist]{hist}}
  for each element of \code{x} and plots the frequencies
  as a stacked barplot using
  \code{\link[graphics:barplot]{barplot}} with \code{beside=FALSE}.}

\value{The function is called for its side effect, producing a barplot
  on the active graphics device. It returns the result of the call to
  \code{\link[graphics:barplot]{barplot}}.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
  x    <- list(rnorm(42), rnorm(42)+2)
  br   <- seq(-3, 5, length=13)
  cols <- c("#1D267B", "#ceffc0")
  histStack(x, breaks=br, col=cols)

  histStack(x, breaks=br, col=cols,
            midsFun=function(z) {
              z <- paste(round(z, 3))
              z[seq(2, length(z), by=2)] <- NA
              return(z)
            } )
}
\keyword{hplot}

