\name{Makesense}
\alias{Makesense}
\title{Produce Smoothed Sense/Anti-sense For All Chromosomes}
\description{
'Makesense' takes a matrix of gene expressions and produces 
 smoothed positive and negative strands for all chromosomes.
}
\usage{
Makesense(expr, lib, f = 1/10)
}
\arguments{
  \item{expr}{a matrix of gene expressions, with genes as rows and columns as samples.}
  \item{lib}{The name of the Bioconductor annotation data package that
    will be used to provide mappings from probes to chromosomal
    locations, such as \code{hgu95av2} or \code{hgu133a}.}
  \item{f}{the smoother span to be passed to 'lowess'. Its value should be in the interval of (0,1). This gives the proportion of points in the plot which influence the smooth at each value. Larger values give more smoothness.}
}
\value{
A list of 2 components:
\item{ans2}{a list, whose components correspond to samples 
in the same order as appearing in the columns of 
'expr'. Each component is also a list, named by chromosomes 
"1"-"22", "X" and "Y". Each named component is again a list 
with two elements named "posS" and "negS", corresponding to the
positive and negative strands of a chromosome, each of which isan object returned by 'lowess'.}
\item{lib}{as passed in.}
}
\author{Robert Gentleman and Xiaochun Li}

\seealso{ \code{\link{plotChr}}}
\examples{
if (!is.installed("hgu133a"))
stop("\nThis example requires the hgu133a package which is not installed.",
     "\nPlease install the package (perhaps using: ",
     "install.packages2(\"hgu133a\")) and try again")

data(eset133a)

esetobj <- Makesense(exprs(eset133a), "hgu133a")
}
\keyword{dplot}
