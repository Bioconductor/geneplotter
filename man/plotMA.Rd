\name{plotMA-methods}
\docType{methods}
\alias{plotMA}
\alias{plotMA-methods}
\alias{plotMA,data.frame-method}
\alias{plotMA,ExonCountSet-method}
\title{Generate an MA plot}

\description{Generate a plot of log fold change versus mean expression (MA plot)
}
\usage{

\S4method{plotMA}{data.frame}( object, ylim = NULL,
  colNonSig = "gray32", colSig = "red3", colLine = "#ff000080",
  log = "x", cex=0.45, xlab="mean expression", ylab="log fold change", ... )
}
\arguments{
   \item{object}{
      A \code{data.frame} with (at least) three columns, the first containing the 
      mean expression values (for the x-axis), the second the logarithmic fold change 
      (for the-y axis) and the third 
      a logical vector indicating significance (for the colouring of the dots).
   }
   \item{ylim}{
      The limits for the y-axis. If missing, an attempt is made to choose a sensible value.
      Dots exceeding the limits will be displayed as triangles at the limits, pointing outwards.
   }
   \item{colNonSig}{
      colour to use for non-significant data points.
   }
   \item{colSig}{
      colour to use for significant data points. 
   }
   \item{colLine}{
      colour to use for the horizontal (y=0) line.
   }
   \item{log}{
      which axis/axes should be logarithmic; will be passed to \code{\link{plot}}.
   }
  \item{cex}{
     The \code{cex} parameter for \code{\link{plot}}.
   }
  \item{xlab}{
     The x-axis label.
   }
  \item{ylab}{
     The y-axis label.
   }
  \item{...}{
     Further parameters to be passed through to \code{\link{plot}}.
   }
}

