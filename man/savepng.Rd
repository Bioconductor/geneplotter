\name{savepng}
\alias{savepdf}
\alias{saveeps}
\alias{savepng}
\alias{savetiff}

\title{Save the contents of the current graphics device to a file}
\description{Save the contents of the current graphics device to file}

\usage{
  savepdf(fn, dir, width=6, asp=1)
  saveeps(fn, dir, width=6, asp=1) 
  savepng(fn, dir, width=480, asp=1)
  savetiff(fn, dir, density=360, keeppdf=TRUE, ...)
}

\arguments{
  \item{fn}{character: name of the output file (without extension). 
    An extension \code{.pdf}, \code{.eps}, \code{.png}, or \code{.tiff}
    will be added automatically.}
  \item{dir}{character: directory to which the file should be written.}
  \item{width}{numeric: width of the image in pixels (png) or inches
    (pdf, eps).}
  \item{asp}{numeric: aspect ratio; height=width*asp.}
  \item{density}{pixels per inch (see Details).}
  \item{keeppdf}{Should the intermediate PDF file (see Details)
    be kept? If \code{FALSE}, it is deleted before the function
    returns.}
  \item{...}{Further arguments that are passed on to \code{savepdf}
    (see Details).}
}

\details{The functions are called for their side effect, writing a
  graphics file.

  \code{savepdf}, \code{savepng}, and \code{saveeps} use the
  devices \code{\link{pdf}}, \code{\link{png}}, and
  \code{\link{postscript}}, respectively.

  There is currently no TIFF device for R, so \code{savetiff}
  works differently. It relies on the external tool \code{convert} from
  the ImageMagick software package. First, \code{savetiff} produces
  a PDF files with \code{savepdf}, then uses \code{\link{system}} to
  invoke \code{convert} with the parameter \code{density}.
  \code{savetiff} does \bold{not} check for the existence of 
  \code{convert} or the success of the system call, and returns silently
  no matter what.
}

\value{Character: name of the file that was written.}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}

\seealso{\code{\link{dev.copy}}, \code{\link{pdf}}, \code{\link{png}},
  \code{\link{postscript}}}

\examples{
  x = seq(0, 20*pi, len=1000)
  plot(x*sin(x), x*cos(x), type="l")

  try({   ## on some machines, some of the devices may not be available
    c(
      savepdf("spiral", dir=tempdir()),
      savepng("spiral", dir=tempdir()),
      saveeps("spiral", dir=tempdir()),
      savetiff("spiral", dir=tempdir())
     )
  })    
}

\keyword{programming}
\keyword{error}
