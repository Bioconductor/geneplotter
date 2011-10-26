multiecdf = function(x, ...)
  UseMethod("multiecdf")
multidensity = function(x, ...)
  UseMethod("multidensity")


multiecdf.formula = function(formula, data = NULL,
  xlab,
  na.action = NULL,
  ...) {

  if(missing(xlab))
     xlab = deparse(substitute(formula))
  if(missing(formula) || (length(formula) != 3))
    stop("'formula' missing or incorrect")
  m = match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data = as.data.frame(data)
  m$... = m$xlab = NULL
  m$na.action = na.action ## force use of default for this method
  m[[1]] = as.name("model.frame")
  mf = eval(m, parent.frame())
  response = attr(attr(mf, "terms"), "response")
  multiecdf(split(mf[[response]], mf[-response]), xlab=xlab, ...)
}

multidensity.formula = function(formula, data = NULL,
  xlab,
  na.action = NULL,
  ...){
  
  if(missing(xlab))
     xlab = deparse(substitute(formula))
  if(missing(formula) || (length(formula) != 3))
    stop("'formula' missing or incorrect")
  m = match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data = as.data.frame(data)
  m$... = m$xlab = NULL
  m$na.action = na.action ## force use of default for this method
  m[[1]] = as.name("model.frame")
  mf = eval(m, parent.frame())
  response = attr(attr(mf, "terms"), "response")
  multidensity(split(mf[[response]], mf[-response]), xlab=xlab, ...)
}


multiecdf.matrix = function(x, xlab, ...) {
  if(missing(xlab))
     xlab = deparse(substitute(x))
  l = lapply(seq_len(ncol(x)), function(j) x[,j])
  names(l) = colnames(x)
  multiecdf(l, xlab=xlab, ...)
}

multidensity.matrix = function(x, xlab, ...) {
  if(missing(xlab))
     xlab = deparse(substitute(x))
  l = lapply(seq_len(ncol(x)), function(j) x[,j])
  names(l) = colnames(x)
  multidensity(l, xlab=xlab, ...)
}

multiecdf.data.frame = function(x, xlab, ...) {
  if(missing(xlab))
     xlab = deparse(substitute(x))
  multiecdf.list(x, xlab=xlab, ...)
}

multidensity.data.frame = function(x, xlab, ...) {
  if(missing(xlab))
     xlab = deparse(substitute(x))
  multidensity.list(x, xlab=xlab, ...)
}

multiecdf.list = function(x,
  xlim,
  col = brewer.pal(9, "Set1"),
  main = "ecdf",
  xlab,
  do.points = FALSE,
  subsample = 1000L,
  legend = list(
    x = "right",
    legend = if(is.null(names(x))) paste(seq(along=x)) else names(x),
    fill = col),
  ...) {
  
  if(missing(xlab))
     xlab = deparse(substitute(x))

  stopifnot(length(x)>=1, length(subsample)==1)
  
  if(is.logical(subsample))
    subsample = if(subsample) 1000L else 0L
  stopifnot(is.numeric(subsample))
  if( (!is.na(subsample)) && (subsample>0) )
    for(i in seq(along=x))
      if(length(x[[i]])>subsample)
        x[[i]] = x[[i]][sample(1:length(x[[i]]), subsample)]
  
  ef = lapply(x, ecdf)
  if(missing(xlim))
    xlim = range(unlist(x), na.rm=TRUE)
  plot(ef[[1]], xlim=xlim, xlab=xlab, main=main, col=col[1], do.points=do.points, ...)
  m = match.call(expand.dots = FALSE) # avoid warnings for invalid arguments
  m$... = m$...[!names(m$...) %in% c("main", "xlab", "ylab", "ylim")]  

  for(j in seq(along=ef)[-1]) {
    mycol = col[1+((j-1)%%length(col))]
    args = c(list(x=ef[[j]], col=mycol, do.points=do.points), m$...)
    do.call(lines, args)
  }

  if(is.list(legend))
    do.call(graphics::legend, legend)
  
  invisible(ef)
}

multidensity.list = function(x,
  bw   = "nrd0",
  xlim,
  ylim,
  col  = brewer.pal(9, "Set1"),
  main = if(length(x)==1) "density" else "densities",
  xlab,
  lty  = 1L ,
  legend = list(
    x = "topright",    
    legend = if(is.null(names(x))) paste(seq(along=x)) else names(x),
    fill = col),     
  ...) {

  ## process argument 'xlab'
  if(missing(xlab))
     xlab = deparse(substitute(x))
  
  ## process argument 'bw': 
  if(length(bw)==1)
    bw = rep(bw, length(x))
  if(length(bw)!=length(x))
    stop("'length(bw)' needs to be either 1 or the same as 'length(x)'.")
  
  ## process argument 'x'
  stopifnot(length(x)>=1)
  if(missing(xlim))
    xlim = range(unlist(x), na.rm=TRUE)
  x = lapply(x, function(z) z[(z>=xlim[1]) & (z<=xlim[2])])
  
  ef = vector(mode="list", length=length(x))
  for(j in seq(along=x))
    ef[[j]] = density(x[[j]], na.rm=TRUE, bw=bw[j])
  
  if(missing(ylim))
    ylim = range(unlist(lapply(ef, "[[", "y")), na.rm=TRUE)
  
  plot(ef[[1]], xlim=xlim, ylim=ylim, xlab=xlab, main=main, col=col[1], lty=lty[1], ...)
  m = match.call(expand.dots = FALSE) ## avoid warnings for invalid arguments
  m$... = m$...[!names(m$...) %in% c("main", "xlab", "ylab", "ylim")]  
  for(j in seq(along=ef)[-1]) {
    args = c(list(x=ef[[j]]), col=col[1+((j-1)%%length(col))],
                               lty=lty[1+((j-1)%%length(lty))], m$...)
    do.call(lines, args)
  }

  if(is.list(legend))
    do.call(graphics::legend, legend)

  invisible(ef)
}

