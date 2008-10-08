\name{groupedHeatmap}
\alias{groupedHeatmap}

\title{Heatmap of a matrix with grouped rows and columns}
\description{The function uses \code{\link[grid]{grid.rect}}
  and \code{\link[grid]{grid.rect}} to draw a heatmap with grouped rows and columns.
}

\usage{
groupedHeatmap(z, frow, fcol,
  fillcolours = c("#2166ac","#4393c3","#92c5de","#d1e5f0","#fefefe","#fddbc7","#f4a582","#d6604d","#b2182b"),
  bordercolour = "#e0e0e0",
  zlim = range(z, na.rm=TRUE))
}

\arguments{
  \item{z}{A matrix with row and column names.}
  \item{frow}{A \code{factor} of length \code{nrow(z)} indicating the
    row grouping.}
  \item{fcol}{A \code{factor} of length \code{ncol(z)} indicating the
    column grouping.}
  \item{fillcolours}{A \code{character} vector of colours from which the
    colour map is obtained through interpolation.} 
  \item{bordercolour}{Either a \code{character} vector of length 1,
    specifying the border colour of the heatmap tiles, or \code{NULL} or
    \code{NA}, which indicates that the border colour should match the
    fill colour.} 
  \item{zlim}{Lower and upper limit of \code{z} values represented in the colour
    map.}
}

\details{The function can be called within other drawing operations
  from the grid package, e.g. within a viewport.
}

\value{The function is called for its side effect, drawing text and
  rectangles on the current viewport.}

\seealso{\code{\link[grid]{grid.text}}, \code{\link[grid]{grid.rect}}}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{

data("mtcars")

groupedHeatmap(
  scale(mtcars),
  frow = factor(sapply(strsplit(rownames(mtcars), " "), "[", 1)),
  fcol = factor(round(seq_len(ncol(mtcars))/3)))
}

