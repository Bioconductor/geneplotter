histStack <- function(x, breaks, breaksFun=paste, ylab="frequency", ...) {
  if(!is.list(x))
    stop("'x' must be a list.")
  bars <- NULL
  for (i in 1:length(x)) {
    if(!is.numeric(x[[i]]))
      paste("Element", i, "of 'x' is not numeric.")
    h <- hist(x[[i]], breaks=breaks, plot=FALSE)
    bars <- rbind(bars, h$counts)
  }

  barplot(bars, names.arg=NULL, space=0, ylab=ylab, ...)

  at     = seq(along=h$breaks)-1
  modulo = ceiling(length(at)/10)
  sel    = (at %% modulo == 0)
  axis(side=1,at=at[sel],labels=breaksFun(h$breaks)[sel])

}
