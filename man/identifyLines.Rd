\name{identifyLines}
\alias{identifyLines}
\title{A function for identifying sample lines on an alongChrom()
  plot. }
\description{
  Given an environment returned from alongChrom(), will call identify()
  on the alongChrom() plot.
}
\usage{
identifyLines(identEnvir,...)
}
\arguments{
  \item{identEnvir}{An environment created by alongChrom.}
  \item{...}{Any extra arguments are passed on to identify().}
}
\details{
  The environment contains two vectors of points - "X" and "Y" are their
  labels, and contain the X and Y points respectively.  These vectors
  are extracted out of the environment and fed into an identify() call
  (with some extra processing to get the labeling correct.)  As with
  identify(), one need only to right click to exit from the function.
}
\author{Jeff Gentry}
\seealso{
  \code{\link{alongChrom}}, \code{\link{identify}}.
}
\examples{
   data(eset)

   ## A bit of a hack to not have a package dependency on hgu95av2
   ## but need to fiddle w/ the warn level to not fail the example anyways.
   curWarn <- options(warn=0)
   on.exit(options(warn=curWarn), add=TRUE)
   if (require(hgu95av2)) {
     z <- buildChromLocation("hgu95av2")
     lty <- c(1, 2, 3, 4, 5)
     cols <- c("red", "green", "blue", "orange", "magenta", "black")
     identEnv <- alongChrom(eset,"1", z, xloc="equispaced",
                          plotFormat="cumulative", scale="none",
                          lty,cols)
     if (interactive()) {
       identifyLines(identEnv)
     }
  } else print("This example needs the hgu95av2 data package")
}
\keyword{utilities}
