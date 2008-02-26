\name{plotExpressionGraph}
\alias{plotExpressionGraph}
\alias{getPlotExpressionColors}
\alias{IMCAEntrezLink}
\alias{defMapFun}
\title{A function to plot a graph colored by expression data}
\description{
  Given a graph and expression data for one entity, will plot the graph
  with the nodes colored according to the expression levels provided.
}
\usage{
plotExpressionGraph(graph, nodeEGmap, exprs, ENTREZIDenvir, mapFun, log = FALSE, nodeAttrs = list(), ...)
}
\arguments{
  \item{graph}{The graph to plot}
  \item{nodeEGmap}{A \code{list} with element names being node names and the
    elements being EntrezLink IDs corresponding to those node names.}
  \item{exprs}{A \code{vector} of expression data, with names being Affymetrix
    IDs and values being the expression level.}
  \item{ENTREZIDenvir}{An \code{environment} mapping Affymetrix IDs to EntrezLink
  IDs, such as the ones provided in the xxx2ENTREZID environments from
  the Bioconductor data packages (where xxx) is a data package).}
  \item{mapFun}{A function to map expression levels to colors.}
  \item{log}{Whether or not the expression data.}
  \item{nodeAttrs}{A \code{list} of node attributes, as per \code{plot.graph}.}
  \item{\dots}{Any extra arguments to be passed to \code{plot.graph}.}
}
\details{
  This function can be used to plot a graph and have the nodes colored
  according to expression levels provided by the user.  The
  \code{graph} parameter is a \code{graph} object from the \code{graph}
  package.

  The \code{nodeEGmap} parameter is a list that maps the nodes of the
  graphs to EntrezLink IDs.  An example of this is the
  \code{IMCAEntrezLink} object in the
  \code{integrinMediatedCellAdhesion} data set in the
  \code{graph} package.

  The \code{exprs} argument is a vector mapping expression levels to
  Affymetrix IDs.  One way to generate an appropriate vector is to
  extract a single column from an \code{ExpressionSet}.

  The \code{ENTREZIDenvir} environment maps Affymetrix IDs to EntrezLink
  IDs.  The simplest way to provide this argument is to load the
  preferred Bioconductor data package (e.g. \code{hgu95av2.db}) and pass in
  that package's \code{xxx2ENTREZID}, where \code{xxx} is the name of the
  package.

  The \code{mapFun} function defaults to the function \code{defMapFun},
  which maps nodes to be either blue, green or red depending for
  expression ranges of 0-100, 101-500, and 501+.  In the case where
  \code{log} is \code{TRUE} these ranges are modified with
  \code{\link{log2}}.  Custom versions of this function can be supplied
  by the user - it must take two parameters, first the expression vector
  and a boolean value (\code{log}) specifying if the data has had a
  \code{log2} applied to it.  The function must return a vector with the
  same names as the expression vector, but the values of the vector will
  be color strings.

  The \code{nodeAttrs} list can be specified if any other node
  attributes are desired to be set by the user.  Please see the
  \code{\link[Rgraphviz]{plot.graph}} man page for more information on this.  The
  other attribute list (\code{attrs} and \code{edgeAttrs}) can be passed
  in via the \code{...} parameter.
  
  The IMCAEntrezLink data structure was created for the purpose of
   illustrating this program.  On Sept 24 2007, the current version
   of \code{hgu95av2.db} was used to map from the nodes of IMCAGraph
   (in graph package) to Entrez identifiers.
}
\author{Jeff Gentry}
\seealso{\code{\link[Rgraphviz]{plot.graph}}, \code{integrinMediatedCellAdhesion}}
\examples{
   if ((require("Rgraphviz")) && (require("hgu95av2.db")) &&
        (require("fibroEset")) && (interactive())) {     
            data(integrinMediatedCellAdhesion)
	    data(IMCAEntrezLink)
            data(fibroEset)
            attrs=getDefaultAttrs()
            attrs$graph$rankdir <- "LR"
            plotExpressionGraph(IMCAGraph, IMCAEntrezLink,
                                exprs(fibroEset)[,1],
                                hgu95av2ENTREZID, attrs=attrs)
   }
}
\keyword{utilities}
\keyword{graphs}
\keyword{hplot}
