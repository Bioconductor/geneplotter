multiecdf = function(x, ...)
  UseMethod("multiecdf")

multiecdf.formula <- function(formula, data = NULL, ..., na.action = NULL)
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

multiecdf.matrix <- function(x, ...)
  multiecdf(x~col(x), ...)
            
multiecdf.default <- function(x, xlim, col, do.points=FALSE,
                              subsample=TRUE, ...) {
  stopifnot(length(x)>=1)
  if(subsample)
    for(i in seq(along=x))
      if(length(x[[i]])>1000)
        x[[i]] <- x[[i]][sample(1:length(x[[i]]), 1000)]
  
  ef = lapply(x, ecdf)
  if(missing(xlim))
    xlim = range(unlist(x), na.rm=TRUE)
  if(missing(col))
    col = brewer.pal(9, "Set1")
  plot(ef[[1]], xlim=xlim, col.hor=col[1], col.vert=col[1], do.points=do.points, ...)
  m <- match.call(expand.dots = FALSE) # avoid warnings for invalid arguments
  m$... <- m$...[!names(m$...) %in% c("main", "xlab", "ylab", "ylim")]  

  for(j in seq(alog=ef)[-1]) {
    mycol = col[1+((j-1)%%length(col))]
    args <- c(list(x=ef[[j]], col.hor=mycol, col.vert=mycol, do.points=do.points), m$...)
    do.call("lines", args)

  }
}

multidensity = function(x, ...)
  UseMethod("multidensity")

multidensity.formula = function(formula, data = NULL, main, xlab, ..., na.action = NULL)
{
  if(missing(formula) || (length(formula) != 3))
    stop("'formula' missing or incorrect")

  if(missing(main))
    main = sprintf("multidensity(%s)", deparse(substitute(formula)))
  if(missing(xlab))
    xlab = deparse(substitute(formula))
  
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- m$main <- m$xlab <- NULL
  m$na.action <- na.action # force use of default for this method
  m[[1]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())
  response <- attr(attr(mf, "terms"), "response")
  multidensity(split(mf[[response]], mf[-response]), main=main, xlab=xlab, ...)
}

multidensity.matrix <- function(x, ...)
  multidensity(x~col(x), ...)

multidensity.default = function(x, bw="nrd0", xlim, ylim, col, main, xlab,  ...) {
  stopifnot(length(x)>=1)

  ef = vector(mode="list", length=length(x))
  ef[[1]] = density(x[[1]], na.rm=TRUE, bw=bw)
  bw = ef[[1]]$bw
  for(j in seq(along=x)[-1])
    ef[[j]] = density(x[[j]], na.rm=TRUE, bw=bw)
  
  if(missing(xlim))
    xlim = range(unlist(x), na.rm=TRUE)
  if(missing(ylim))
    ylim = range(unlist(lapply(ef, "[[", "y")), na.rm=TRUE)
  if(missing(col))
    col = brewer.pal(9, "Set1")
  if(missing(main))
    main = sprintf("multidensity.default(%s)", deparse(substitute(x)))
  if(missing(xlab))
    xlab = deparse(substitute(x))
  
  plot(ef[[1]], xlim=xlim, ylim=ylim, xlab=xlab, main=main, col=col[1],  ...)
  m <- match.call(expand.dots = FALSE) ## avoid warnings for invalid arguments
  m$... <- m$...[!names(m$...) %in% c("main", "xlab", "ylab", "ylim")]  
  for(j in seq(along=ef)[-1]) {
    args <- c(list(x=ef[[j]], col=col[1+((j-1)%%length(col))]), m$...)
    do.call("lines", args)
  }
  invisible(ef)
}

