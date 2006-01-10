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

multiecdf.default <- function(x, xlim, col, do.points=FALSE,
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
  plot(ef[[1]], xlim=xlim, col.hor=col[1], col.vert=col[1], do.points=do.points, ...)
  m <- match.call(expand.dots = FALSE) # avoid warnings for invalid arguments
  m$... <- m$...[!names(m$...) %in% c("main", "xlab", "ylab", "ylim")]  

  for(j in 2:length(ef)) {
    mycol = col[1+((j-1)%%length(col))]
    args <- c(list(x=ef[[j]], col.hor=mycol, col.vert=mycol, do.points=do.points), m$...)
    do.call("lines", args)

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

multidensity.default = function(x, xlim, col, ...) {
  stopifnot(length(x)>=1)
  ef = lapply(x, density)
  if(missing(xlim))
    xlim = range(unlist(x))
  if(missing(col))
    col = brewer.pal(9, "Set1")
  plot(ef[[1]], xlim=xlim, col=col[1],  ...)
  m <- match.call(expand.dots = FALSE) # avoid warnings for invalid arguments
  m$... <- m$...[!names(m$...) %in% c("main", "xlab", "ylab", "ylim")]  
  for(j in 2:length(ef)) {
    args <- c(list(x=ef[[j]], col=col[1+((j-1)%%length(col))]), m$...)
    do.call("lines", args)
  }
}

