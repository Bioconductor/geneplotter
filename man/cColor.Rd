\name{cColor}
\alias{cColor}
\title{A function for marking specific probes on a cPlot. }
\description{
  Given a set of probes, will highlight them in the color desired on
  a plot which has already been created via the function cPlot().
}
\usage{
cColor(probes, color, plotChroms, scale=c("relative","max"),glen=0.4)
}
\arguments{
  \item{probes}{The probes that are being highlighted.}
  \item{color}{The color to highlight the probes.}
  \item{plotChroms}{An object of type \code{chromLocation} which contains all
    the gene information to be plotted. }
  \item{scale}{Whether to plot the graph scaled absolutely or relative
    by chromosome.  Default is absolute.}
  \item{glen}{The length of the gene line plotted}
}
\details{
  It is important to call the function \code{cPlot()} first.  This function
  will then search for the specific locations of the probes desired,
  which are contained within the \code{plotChroms} instance of a
  \code{chromLocation} class.  It will then pass these on to the
  plotting routine to highlight the desired locations.  NOTE:  It
  is important that \code{plotChroms}, \code{scale} and \code{glen}
  parameters are the same as used for \code{cPlot()}.
}
\author{Jeff Gentry}
\seealso{
  \code{\link{cPlot}}, \code{\link{chromLocation-class}}
}
\examples{
     z <- buildChromLocation("hgu95av2")
     cPlot(z)
     probes <- c("266_s_at", "31411_at", "610_at", "failExample")
     cColor(probes, "red", z)
     probes2 <- c("960_g_at", "41807_at", "931_at", "39032_at")
     cColor(probes2, "blue", z)
}
\keyword{utilities}

