multiecdf = function(x, ...)
  UseMethod("multiecdf")

multiecdf.formula <-
    function(formula, data = NULL, ..., na.action = NULL)
{
    if(missing(formula) || (length(formula) != 3))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action # force use of default for this method
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    multiecdf(split(mf[[response]], mf[-response]), ...)
}

multiecdf.default <- function(x, xlim, col, main="", do.points=FALSE,
                              subsample=TRUE, ...) {
  stopifnot(length(x)>=1)
  if(subsample)
    for(i in 1:length(x))
      if(length(x[[i]])>1000)
        x[[i]] <- x[[i]][sample(1:length(x[[i]]), 1000)]
  ef = lapply(x, ecdf)
  if(missing(xlim))
    xlim = range(unlist(x))
  if(missing(col))
    col = brewer.pal(9, "Set1")
  plot(ef[[1]], xlim=xlim, col.hor=col[1], col.vert=col[1], main=main, do.points=do.points, ...)
  for(j in 2:length(ef)) {
    mycol = col[1+((j-1)%%length(col))]
    lines(ef[[j]], col.hor=mycol, col.vert=mycol, do.points=do.points, ...)
  }
}

multidensity = function(x, ...)
  UseMethod("multidensity")

multidensity.formula <-
    function(formula, data = NULL, ..., na.action = NULL)
{
    if(missing(formula) || (length(formula) != 3))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action # force use of default for this method
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    multidensity(split(mf[[response]], mf[-response]), ...)
}

multidensity.default = function(x, xlim, col, main="", ...) {
  stopifnot(length(x)>=1)
  ef = lapply(x, density)
  if(missing(xlim))
    xlim = range(unlist(x))
  if(missing(col))
    col = brewer.pal(9, "Set1")
  plot(ef[[1]], xlim=xlim, col=col[1], main=main, ...)
  for(j in 2:length(ef)) {
    lines(ef[[j]], col=col[1+((j-1)%%length(col))], ...)
  }
}

