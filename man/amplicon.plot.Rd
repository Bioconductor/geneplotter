\name{amplicon.plot}
\alias{amplicon.plot}
\title{Create an amplicon plot  }
\description{
  Given a two-sample test statistic and an exprSet this function plots
  regions of the genome that are either highly expressed (in red) or
  have low expression (blue) differentially in the two groups.
}
\usage{
amplicon.plot(ESET, FUN, genome)
}
\arguments{
  \item{ESET}{an object of class \code{exprSet}}
  \item{FUN}{A two sample test function suitable for \code{\link{esApply}}. }
  \item{genome}{A character string of the base name for the annotation. }
}
\details{
  In some genetic studies we are interested in finding regions of the
  genome where there are a set of highly expressed genes in some
  subgroup of the population. This set of highly (or lowly) expressed
  genes is often of great interest. For example in breast cancer the
  HER--2 gene is on an amplicon. In some patients approximately 5 genes
  located near HER--2 are all amplified.

  These plot should help in the search for such regions.
}
\value{
  No value is returned. This function is executed purely for side
  effect.
}

\author{Robert Gentleman }

\seealso{\code{\link{esApply}}, \code{\link{make.chromOrd}}}

\examples{
  ##none yet; takes too long
}
\keyword{ hplot }
