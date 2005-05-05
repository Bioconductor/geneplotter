.pradaMakeX    <- function(x, y) {
  ## check correct input
  if (!is.matrix(x) || ncol(x)<2)
    x <- cbind(x, y)
  if (!is.matrix(x) && ncol(x)==2)
    stop("Wrong size of arguments 'x' (and/or 'y')")
  return(x)
}

.pradaMakeDens <- function(x, nbin, bandwidth) {
  ##for(p in c("KernSmooth", "RColorBrewer"))
  ##  if(!require(p, character.only=TRUE))
  ##    stop(paste("Required package '", p, "' could not be found", sep=""))
  
  if (length(nbin) == 1)
    nbin <- c(nbin, nbin)
  if (!is.numeric(nbin) || (length(nbin)!=2))
    stop("'nbin' must be numeric of length 1 or 2")

  if (missing(bandwidth)) {
    bandwidth <- diff(apply(x, 2, quantile, probs=c(0.05, 0.95))) / 25
  } else {
    if(!is.numeric(bandwidth))
      stop("'bandwidth' must be numeric")
  }
  ## create density map
  rv <- bkde2D(x, gridsize=nbin, bandwidth=bandwidth)
  rv$bandwidth <- bandwidth
  return(rv)
}

smoothScatter <- function(x, y, 
                          nbin=128,
                          bandwidth,
                          colramp=colorRampPalette(brewer.pal(9, "YlGnBu")),
                          nrpoints=100,
                          transformation=function(x) x^.25,
                          xlab, ylab, ...) {
  
  if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) > 1) )
    stop("'nrpoints' should be integer or inf")
  if(missing(xlab)) {
    xlab <- if(missing(y))
      colnames(x)[1]
    else
      deparse(substitute(x))
  }
  if(missing(ylab)) {
    ylab <- if(missing(y))
      colnames(x)[2]
    else
      deparse(substitute(y))
  }

  ## create density map
  x    <- .pradaMakeX(x, y)
  map  <- .pradaMakeDens(x, nbin, bandwidth)
  xm   <- map$x1
  ym   <- map$x2
  dens <- map$fhat
  dens <- array(transformation(dens), dim=dim(dens))
  
  ## plot color image
  image(xm, ym, z=dens, col=colramp(256), xlab=xlab, ylab=ylab, ...)
  box()
  
  ## plot selection of dots
  if (nrpoints!=0){
    ## we assume that map$x1 and map$x2 go linearly from
    ## their first to their last value in nbin steps
    stopifnot(length(xm)==nrow(dens), length(ym)==ncol(dens))
    ixm <- round((x[,1]-xm[1])/(xm[length(xm)]-xm[1])*(length(xm)-1))
    iym <- round((x[,2]-ym[1])/(ym[length(ym)]-ym[1])*(length(ym)-1))
    idens <- dens[1 + iym*length(xm) + ixm]
    nrpoints <- min(nrow(x), ceiling(nrpoints))
    sel <- order(idens, decreasing=FALSE)[1:nrpoints]
    points(x[sel,1:2], pch=".", col="black")
  }
}

densCols <- function(x, y,
                     nbin=128,
                     bandwidth,
                     colramp=colorRampPalette(brewer.pal(9, "YlGnBu")[-(1:3)])) {

  ## create density map 
  x    <- .pradaMakeX(x, y)
  map  <- .pradaMakeDens(x, nbin, bandwidth)

  ## bin x-values
  dx   <- diff(range(x[,1])) / (length(map$x1)-1)
  xbin <- floor(1 + (x[,1] - min(x[,1])) / dx)

  ## bin y-values
  dy   <- diff(range(x[,2])) / (length(map$x2)-1)
  ybin <- floor(1 + (x[,2] - min(x[,2])) / dy)
    
  dens <- map$fhat[xbin + nrow(map$fhat) * (ybin-1)]
  dens[is.na(dens)]<- 0

  ## transform densities to colors
  colpal <- as.integer(1 + (length(dens)-1) * dens/max(dens))
  cols   <- colramp(length(dens))[colpal]
    
  return(cols)
}
