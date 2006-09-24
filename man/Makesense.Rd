\name{Makesense}
\alias{Makesense}
\alias{Makesense,matrix,character-method}
\alias{Makesense,exprSet,character-method}
\alias{Makesense,exprSet,missing-method}
\alias{Makesense,eSet,missing-method}
\alias{Makesense,eSet,character-method}
\title{Produce Smoothed Sense/Anti-sense For All Chromosomes}
\description{
  'Makesense' takes either an exprSet object or a matrix of gene
  expressions and will produce a smoothed positive and negative strands
  for all chromosomes.
}
\usage{
Makesense(expr, lib, ...)
}
\arguments{
  \item{expr}{One of an \code{exprSet}, an \code{eSet} or a \code{matrix} of 
    gene expressions with genes as rows and columns as samples.}
  \item{lib}{The name of the Bioconductor annotation data package that
    will be used to provide mappings from probes to chromosomal
    locations, such as \code{hgu95av2} or \code{hgu133a}.  If
    \code{expr} is an \code{exprSet}, the argument defaults to the
    \code{annotation} slot of the \code{exprSet}.}
  \item{...}{Currently, the only optional argument is \code{f}, the
    smoother span to be passed to 'lowess'. Its value should
    be in the interval of (0,1). This gives the proportion of points in
    the plot which influence the smooth at each value. Larger values
    give more smoothness.  The default value for this argument is 1/10.}
}
\details{
  The \code{expr} argument can either be of class \code{exprSet} or
  \code{matrix}, where the latter represents the matrix of gene
  expressions.

  If the \code{expr} argument is an \code{exprSet}, the \code{lib}
  argument will use the \code{annotation} slot.  Users can override this
  behaviour and supply their own \code{lib} argument if they wish.  If
  the \code{exprSet} has no value associated with the \code{annotation}
  slot (which should not happen, but is possible) then the user must
  supply the \code{lib} argument manually or the function will throw an
  error.
  }
\value{
A list of 2 components:
\item{ans2}{a list, whose components correspond to samples 
in the same order as appearing in the columns of 
'expr'. Each component is also a list, named by chromosomes 
"1"-"22", "X" and "Y". Each named component is again a list 
with two elements named "posS" and "negS", corresponding to the
positive and negative strands of a chromosome, each of which is an object returned by 'lowess'.}
\item{lib}{A string giving the name of the annotation data package to
use.  Optional if \code{expr} is an \code{exprSet}.}
}
\author{Robert Gentleman and Xiaochun Li}

\seealso{ \code{\link{plotChr}}}
\examples{
if(require("hgu133a") )
{
  data(expressionSet133a)

  esetobj <- Makesense(exprs(expressionSet133a), "hgu133a")

  esetobj2 <- Makesense(expressionSet133a[1:200, ])
}
}
\keyword{dplot}
