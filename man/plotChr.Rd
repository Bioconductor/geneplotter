\name{plotChr}
\alias{plotChr}
\title{Plot Smoothed Sense/Anti-sense of Specified Chromosomes}
\description{
For a given chromosome, plot the smooths of the sense and the 
anti-sense from 5' to 3' (left to right on x-axis). 
}
\usage{
plotChr(chrN, senseObj, cols = rep("blue", length(senseObj[[1]])), log = FALSE, xloc = c("equispaced", "physical"), geneSymbols = FALSE, ngenes = 20, lines.at = NULL, lines.col = "red")
}
\arguments{
  \item{chrN}{The desired chromosome, e.g. for humans it would be a character string in the set of c(1:22, "X", "Y").}
  \item{senseObj}{The result of \code{Makesense}.}
  \item{cols}{A vector of colors for the lines in the plot, typically specified according to a certain pheotype of samples.}
  \item{log}{Logical, whether log-transformation should be taken on the smoothed expressions.}
  \item{xloc}{Determines whether the "Representative Genes" 
	will be displayed according to their relative 
	positions on the chromosome (physical), or spaced 
	evenly (equispaced).  Default is equispaced.}
  \item{geneSymbols}{Logical, whether to use Affy IDs or 
	Gene Symbols for "Representative Genes", default is 
	Affy IDs.}
  \item{ngenes}{Desired number of "Representative Genes". The 
	number of actual displayed genes may differ.}
  \item{lines.at}{A vector of Affy IDs. Vertical lines will 
	be drawn at specified genes.}
  \item{lines.col}{A vector of colors associated with 
	\code{lines.at}.}
}
\author{Robert Gentleman and Xiaochun Li}
\seealso{\code{\link{Makesense}}}
\examples{
example(Makesense)

op <- par(ask=TRUE)

cols <- c("red", "green")[unclass(eset133a$cov1)]
plotChr("21", esetobj, cols)

# plot on log-scale:

plotChr("21", esetobj, cols, log=TRUE)

# genesymbol instead of probe names:

plotChr("21", esetobj, cols, log=TRUE, geneSymbols=TRUE)

# add vertical lines at genes of interest:

gs <- c("220372_at", "35776_at", "200943_at")
plotChr("21", esetobj, cols, log=TRUE, geneSymbols=FALSE, lines.at=gs)

# add vertical lines at genes of interest
# with specified colors:

gs <- c("220372_at", "35776_at", "200943_at")
cc <- c("blue", "cyan","magenta")
plotChr("21", esetobj, cols, log=TRUE, geneSymbols=FALSE, lines.at=gs,
lines.col=cc)
par(op)
}
\keyword{hplot} 
