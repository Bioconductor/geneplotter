\name{imageMap-methods}
\docType{methods}

\alias{imageMap}
\alias{imageMap-methods}
\alias{imageMap,matrix-method}
\alias{imageMap,matrix,connection,list,character-method}

\title{Write an HTML IMG tag together with a MAP image map.}
\description{Write an HTML IMG tag together with a MAP image map.}
\usage{
  \S4method{imageMap}{matrix,connection,list,character}(object, con, tags, imgname)
}
\arguments{
  \item{object}{Matrix with 4 columns, specifying the coordinates
    of the mouse-sensitive region . Each row specifies the corners of a 
    rectangle within the image, in the following order: (left x,
    lower y, right x, upper y). Note that the point (x=0, y=0) is 
    at the left upper side of the image.}
  \item{con}{Connection to which the image map is written.}
  \item{tags}{Named list whose elements are named character vectors.
    Names must correspond to node names in \code{object}. See details.}
  \item{imgname}{Character. Name of the image file (for example PNG
    file) that contains the plot.}
}

\details{The most important tags are \code{TITLE}, \code{HREF},
  and \code{TARGET}. If the list \code{tags} contains an element
  with name \code{TITLE}, then this must be a named character vector
  containing the tooltips that are to be displayed when the mouse moves
  over a node. The names of the nodes are specified in the \code{names}
  attribute  of the character vector and must match those of
  \code{object}.
  
  Similarly, \code{HREF} may be used to specify hyperlinks that the
  browser can follow when the mouse clicks on a node, and \code{TARGET}
  to specify the target browser window.

  Currently, only rectangular regions are implemented; the actual
  shape of the nodes as specified in \code{object} is ignored.
  Also, tags for edges of the graph are currently not supported.

  This function is typically used with the following sequence
  of steps:
  \enumerate{
    \item generate your graphic and save it as a bitmap file, e.g.
    using the \code{\link{jpeg}}, \code{\link{png}}, or
    \code{\link{bitmap}} devices. At this stage, you also need to
    figure out the pixel coordinates of the interesting regions
    within your graphic. Since the mapping between device coordinates
    and pixel coordinates is not obvious, this may be a little tricky.
    See the examples below, and for a more complex example, see the
    source code of the function \code{\link{plotPlate}}.
    \item open an HTML page for writing and write HTML header,
    e.g. using the \code{\link{openHtmlPage}} function.
    \item Call the \code{\link{imageMap}} function.
    \item Optionally, write further text into the HTML connection.
    \item Close HTML file, e.g. using the \code{\link{closeHtmlPage}} function.
  }
}

\value{The function is called for its side effect, which is writing text into
the connection \code{con}.}

\seealso{\code{\link[prada::plotPlate]{plotPlate}}, 
   \code{\link[base::writeLines]{writeLines}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{
f1  = paste(tempfile(), ".html", sep="")
f2  = paste(tempfile(), ".html", sep="")
fpng = tempfile()

if(capabilities()["png"] && interactive()) {
  ## create the image
  width  = 512
  height = 80
  png(fpng, width=width, height=height)
  par(mai=rep(0,4))
  plot(0, xlim=c(0,width-1),ylim=c(0,height-1),xaxs="i",yaxs="i",type="n",bty="n")
  coord=cbind(c( 12,200,378),
              c(  4, 28, 14),
              c(126,300,506),
              c( 56, 67, 47))
  rect(coord[,1], height-coord[,2], coord[,3], height-coord[,4],
       col=c("#1f77b4", "#b2df8a","#E41A1C"))
  dev.off()

  ## create the frame set
  cat("<html><head><title>Hello world</title></head>\n",
      "<frameset rows=\"120,*\" border=\"0\">\n",
      "<frame name=\"banner\" src=\"file://", f2, "\">\n",
      "<frame name=\"main\" scrolling=\"auto\">",
      "</frameset>", sep="",file=f1)

  ## create the image map
  con = file(f2, open="w")
  imageMap(coord, con,
    list(HREF=paste("http://", c("www.bioconductor.org", "www.r-project.org","www.ebi.ac.uk"), sep=""),
         TITLE=c("Blue box", "Green box", "Red box"),
         TARGET=c("main", "main","new")), fpng)
  close(con)

  cat("Now have a look at file ", f1, " with your browser.\n", sep="")
}
}


