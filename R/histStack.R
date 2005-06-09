histStack <- function(x, breaks, midsFun=paste, ...) {
  if(!is.list(x))
    stop("'x' must be a list.")
  bars <- mids <- NULL
  for (i in 1:length(x)) {
    if(!is.numeric(x[[i]]))
      paste("Element", i, "of 'x' is not numeric.")
    h <- hist(x[[i]], breaks=breaks, plot=FALSE)
    bars <- rbind(bars, h$counts)
    if(is.null(mids))
      mids <- h$mids
    else
      stopifnot(identical(mids, h$mids))
  }

  mids <- midsFun(mids)
  mb <- max(colSums(bars))
  ylim <- c(0,8)
  if(mb>8)
    ylim <- c(0,12)
  if(mb>12)
    ylim <- c(0,mb)
  barplot(bars, names.arg=mids, ylim=ylim, ...)  
}
