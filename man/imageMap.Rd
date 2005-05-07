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
    source code of the function \code{\link[prada]{plotPlate}}.
    \item open an HTML page for writing and write HTML header,
    e.g. using the \code{\link{openHtmlPage}} function.
    \item Call the \code{\link{imageMap}} function.
    \item Optionally, write further text into the HTML connection.
    \item Close HTML file, e.g. using the \code{\link{closeHtmlPage}} function.
  }
}

\value{The function is called for its side effect, which is writing text into
the connection \code{con}.}

\seealso{\code{\link[prada]{plotPlate}}, 
   \code{\link[base]{writeLines}}}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{manip}
\examples{
f1  = paste(tempfile(), ".html", sep="")
f2  = paste(tempfile(), ".html", sep="")
fpng = tempfile()

if(capabilities()["png"] && interactive()) {
  ## create the image
  colors = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")
  width  = 512
  height = 256
  png(fpng, width=width, height=height)
  par(mai=rep(0,4))
  plot(0,xlim=c(0,width-1),ylim=c(0,height-1),xaxs="i",yaxs="i",type="n",bty="n")
  cx=floor(runif(100)*(width-11))
  cy=floor(runif(100)*(height-11))
  coord=cbind(cx, cy, cx+10, cy+10)
  rect(coord[,1], height-coord[,2], coord[,3], height-coord[,4],
       col=sample(colors, 100, replace=TRUE))
  text(width/2, height-3, "Klick me!", adj=c(0.5, 1), font=2)
  dev.off()

  ## create the frame set
  cat("<html><head><title>Hello world</title></head>\n",
      "<frameset rows=\"280,*\" border=\"0\">\n",
      "<frame name=\"banner\" src=\"file://", f2, "\">\n",
      "<frame name=\"main\" scrolling=\"auto\">",
      "</frameset>", sep="",file=f1)

  ## create the image map
  href  =sample(c("www.bioconductor.org", "www.r-project.org"),nrow(coord),replace=TRUE)
  title =sample(as.character(packageDescription("geneplotter")),nrow(coord),replace=TRUE)
  con = file(f2, open="w")
  imageMap(coord, con,
    list(HREF=paste("http://", href, sep=""),
         TITLE=title, TARGET=rep("main", nrow(coord))), fpng)
  close(con)

  cat("Now have a look at file ", f1, " with your browser.\n", sep="")
}
}


