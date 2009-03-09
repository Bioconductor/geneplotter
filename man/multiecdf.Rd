\name{multiecdf}
\alias{multiecdf}
\alias{multiecdf.list}
\alias{multiecdf.formula}
\alias{multiecdf.matrix}
\alias{multidensity}
\alias{multidensity.list}
\alias{multidensity.formula}
\alias{multidensity.matrix}

\title{Multiple empirical cumulative distribution functions (ecdf) and densities}

\description{Plot multiple empirical cumulative distribution functions (ecdf)
  and densities with user interface similar to that of \code{\link{boxplot}}.}

\usage{
multiecdf(x, \dots)
\method{multiecdf}{formula}(formula, data = NULL, xlab, na.action = NULL, \dots)
\method{multiecdf}{matrix}(x, xlab, ...) 
\method{multiecdf}{list}(x,
          xlim,
          col = brewer.pal(9, "Set1"),
          main = "ecdf",
          xlab,
          do.points = FALSE,
          subsample = 1000L,
          legend = list(
            x = "right",
            legend = if(is.null(names(x))) paste(seq(along=x)) else names(x),
            fill = col),
          \dots)

multidensity(x, \dots)
\method{multidensity}{formula}(formula, data = NULL, xlab, na.action = NULL, \dots)
\method{multidensity}{matrix}(x, xlab, ...) 
\method{multidensity}{list}(x,
             bw = "nrd0",
             xlim,
             ylim,
             col  = brewer.pal(9, "Set1"),
             main = if(length(x)==1) "density" else "densities",
             xlab,
             lty  = 1L,
             legend = list(
               x = "topright",
               legend = if(is.null(names(x))) paste(seq(along=x)) else names(x),
               fill = col),
             \dots)
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
  \item{x}{a list of numeric vectors.}
  \item{bw}{the smoothing bandwidth, see the manual page for
    \code{\link{density}}. If \code{bw} is a character string specifying a rule to choose
    the bandwidth, this rule is applied to \code{x[[1]]} and then
    the same numerical value of \code{bw} is used throughout.}    
  \item{xlim}{Range of the x axis. If missing, the data range is used.}
  \item{ylim}{Range of the y axis. If missing, the range of the density
    estimates is used.}
  \item{col, lty}{Line colors and line type.}
  \item{main}{Plot title.}
  \item{xlab}{x-axis label.}
  \item{do.points}{logical; if \code{TRUE}, also draw points at the knot
    locations.}
  \item{subsample}{numeric or logical of length 1. If numeric, and
   larger than 0, subsamples of that size are used to compute and plot
   the ecdf for those elements of \code{x} with more than that number of
   observations. If logical and \code{TRUE}, a value of 1000 is used for
   the subsample size.}
  \item{legend}{a list of arguments that is passed to the function \code{\link[graphics::legend]{legend}}.}
  \item{...}{Further arguments that get passed on to the \code{plot} functions.}
}

\seealso{
  \code{\link[graphics::boxplot]{boxplot}},
  \code{\link[stats::ecdf]{ecdf}},
  \code{\link[stats::density]{density}}
}

\details{The usefulness of \code{multidensity} can vary,
  depending on the data and because of smoothing
  artifacts. \code{multiecdf} will in many cases be preferable.}

\value{For the \code{multidensity} functions, a list of \code{\link{density}} objects.}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
  words = strsplit(packageDescription("geneplotter")$Description, " ")[[1]]
  factr = factor(sample(words, 2000, replace = TRUE))
  x = rnorm(length(factr), mean=as.integer(factr))
  
  multiecdf(x ~ factr)
  multidensity(x ~ factr)
}
\keyword{hplot}

