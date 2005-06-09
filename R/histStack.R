histStack <- function(x, breaks, breakFun=paste, ylab="frequency", ...) {
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

  at     = seq(along=h$breaks)
  modulo = floor(length(at)/10)
  sel    = ((at-1) %% modulo == 0)
  axis(side=1,at=at[sel],labels=breakFun(h$breaks)[sel])
  browser()
}
