\name{colorRampPalette}
\alias{colorRampPalette}
\alias{colorRamp}
\title{Extend a color palette to a color ramp }
\description{
These functions are useful for converting hand-designed `sequential' or `diverging' color schemes into continous color ramps eg for image and filled contour plots.
}
\usage{
colorRampPalette(palette, bias = 1,method=c("linear","spline"))
colorRamp(palette, bias=1,method=c("spline","linear"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{palette}{A vector of colors}
  \item{bias}{Bias >1 puts more colors at high values, <1 puts more at low values}
\item{method}{Interpolation method. Spline interpolation may give smoother results but can end up outside the valid RGB space, especially with highly saturated palettes or those based on a small starting set.}
}

\value{
\code{colorRampPalette} returns a palette generating function like \code{\link{terrain.colors}} or \code{\link{heat.colors}} that takes an integer argument and generates a palette with that many colors.

\code{colorRamp} returns a color generating function, that takes a number between 0 and 1 and generates the corresponding color.

}
\author{Thomas Lumley}

\seealso{\url{www.colorbrewer.org} }

\examples{
   m<-outer(1:20,1:20,function(x,y) sin(sqrt(x*y)/3))

   redblue<-colorRampPalette(c("red","orange","blue"))
   filled.contour(m,color.palette=redblue)

   yb<-colorRampPalette(c("yellow2","goldenrod","darkred"))
   filled.contour(m,color.palette=yb)

   yb<-colorRampPalette(c("yellow2","goldenrod","darkred"),bias=0.5)
   filled.contour(m,color.palette=yb)

   ## Two ColorBrewer sequential palettes
   ylOrBn4<-c("#FFFFD4", "#FED98E", "#FE9929","#CC4C02")
   ylGnBl5<-c("#FFFFCC","#C7E9B4","#7FCDBB","#40B6C4","#2C7FB8" ,"#253494")

   ylOrBn<-colorRampPalette(ylOrBn4)
   ylGnBl<-colorRampPalette(ylGnBl5,method="spline")
   filled.contour(m,color.palette=ylOrBn)
   filled.contour(m,color.palette=ylGnBl)

   if(require("RColorBrewer")){
	RdBu7<-colorRampPalette(brewer.pal(7,"RdBu"),method="spline")
        RWB<-colorRampPalette(c("red","white","blue"))
        redblue<-colorRampPalette(c("red","blue"))
	filled.contour(m,color.palette=redblue,main="naive red-blue")
        filled.contour(m,color.palette=RWB,main="naive red-white-blue")
        filled.contour(m,color.palette=RdBu7,main="ColorBrewer RdBu7")
   }

}
\keyword{color}% at least one, from doc/KEYWORDS
