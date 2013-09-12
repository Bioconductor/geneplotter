setMethod( "plotMA", signature( object="data.frame" ),
function( object, ylim = NULL,
  colNonSig = "gray32", colSig = "red3", colLine = "#ff000080",
  log = "x", cex=0.45, xlab="mean expression", ylab="log fold change", ... )
{
   if( !( ncol(object) == 3 & inherits( object[[1]], "numeric" ) & inherits( object[[2]], "numeric" )
         & inherits( object[[3]], "logical" ) ) ) {
      stop( "When called with a data.frame, plotMA expects the data frame to have 3 columns, two numeric ones for mean and log fold change, and a logical one for significance.")
   }
   colnames(object) <- c( "mean", "lfc", "sig" )
   object = subset( object, mean != 0 )
   py = object$lfc
   if( is.null(ylim) )
      ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
   plot(object$mean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=ifelse( object$sig, colSig, colNonSig ), xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline( h=0, lwd=4, col=colLine )
}
)
