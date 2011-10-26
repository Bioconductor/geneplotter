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
  and densities with a user interface similar to that of \code{\link{boxplot}}.
  The usefulness of \code{multidensity} is variable, depending on the data and because of smoothing
  artifacts. \code{multiecdf} will in many cases be preferable. Please
  see Details.
}

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
  \item{x}{methods exist for: \code{formula}, \code{matrix}, \code{data.frame}, \code{list} of numeric vectors.}
  \item{bw}{the smoothing bandwidth, see the manual page for
    \code{\link{density}}. The length of \code{bw} needs to be either 1
    (in which case the same is used for all groups)
    or the same as the number of groups in \code{x} (in which case the
  corresponding value of \code{bw} is used for each group).}
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
  \item{legend}{a list of arguments that is passed to the function
  \code{\link[graphics]{legend}}.} 
  \item{...}{Further arguments that get passed on to the \code{plot} functions.}
}

\seealso{
  \code{\link[graphics]{boxplot}},
  \code{\link[stats]{ecdf}},
  \code{\link[stats]{density}}
}

\details{The choice of the smoothing bandwidths in \code{multidensity}
  can be problematic, in particular, if the different groups vary with
  respect to range and/or number of data points. If curves look
  excessively wiggly or overly smooth, try varying the arguments
  \code{xlim} and \code{bw}; note that the argument \code{bw} can be a
  vector, in which case it is expect to align with the groups.}

\value{For the \code{multidensity} functions, a list of
  \code{\link{density}} objects.} 

\author{Wolfgang Huber}

\examples{
  words = strsplit(packageDescription("geneplotter")$Description, " ")[[1]]
  factr = factor(sample(words, 2000, replace = TRUE))
  x = rnorm(length(factr), mean=as.integer(factr))
  
  multiecdf(x ~ factr)
  multidensity(x ~ factr)
}
\keyword{hplot}

