\name{multiecdf}
\alias{multiecdf}
\alias{multiecdf.default}
\alias{multiecdf.formula}

\alias{multidensity}
\alias{multidensity.default}
\alias{multidensity.formula}

\title{Multiple empirical cumulative distribution functions (ecdf)
  and densities}

\description{Plot multiple empirical cumulative distribution functions (ecdf)
  and densities with user interface similar to that of \code{\link{boxplot}}.}

\usage{
multiecdf(x, \dots)
\method{multiecdf}{formula}(formula, data = NULL, \dots, na.action = NULL)
\method{multiecdf}{default}(x, xlim, col, main="", do.points=FALSE, \dots)

multidensity(x, \dots)
\method{multidensity}{formula}(formula, data = NULL, \dots, na.action = NULL)
\method{multidensity}{default}(x, xlim, col, main="", \dots)

}

\arguments{
  \item{formula}{a formula, such as \code{y ~ grp}, where \code{y} is a
    numeric vector of data values to be split into groups according to
    the grouping variable \code{grp} (usually a factor).}
  \item{data}{a data.frame (or list) from which the variables in
    \code{formula} should be taken.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is to ignore missing
    values in either the response or the group.}
  \item{x}{A list of numeric vectors.}
  \item{xlim}{Range of the x axis.}
  \item{col}{Line colors. If missing, the
    \code{\link[RColorBrewer::brewer.pal]{Set1}} palette from RColorBrewer
  is used.}
  \item{main}{Plot title.}
  \item{do.points}{logical; if \code{TRUE}, also draw points at the knot locations.}
  \item{...}{Further arguments that get passed on.}
}

\seealso{
  \code{\link[graphics::boxplot]{boxplot}},
  \code{\link[stats::ecdf]{ecdf}}
  \code{\link[stats::density]{density}}
}
\details{}
\value{}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
  f = 1 + (runif(1000)>0.5)
  x = rnorm(length(f), mean=f, sd=f)
  
  multiecdf(x~f)
  multidensity(x~f)
}
\keyword{hplot}

