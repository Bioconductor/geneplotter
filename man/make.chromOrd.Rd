\name{make.chromOrd}
\alias{make.chromOrd}
\title{Make a chromOrd object  }
\description{
  This function makes a chromOrd object.
}
\usage{
make.chromOrd(genome, gnames)
}
\arguments{
  \item{genome}{A character string. }
  \item{gnames}{A character vector of the genes to be selected. }
}
\details{
  This function reads in a lot of annotation data and creates a list
  with one element for each chromosome. The elements of this list are
  indices indicating the order of the genes that are on that chromosome
  (and in the annotation data set being used).
}
\value{
 A list of chromOrd type. One element for each chromosome. Suitable for
 reordering other values according to the chromosomal location.
}

\author{ Robert Gentleman }

\seealso{ \code{\link{amplicon.plot}} }

\examples{
  data(eset)
  make.chromOrd("hgu95A", geneNames(eset))
}
\keyword{ utilities }
