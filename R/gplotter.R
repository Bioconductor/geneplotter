#copyright R. Gentleman, 2001, all rights reserved
#functions/ methods to plot microarray data

#Cheng Li's Red/Blue color scheme

dChip.colors <- function(n) GetColor(seq(-3,3,6/n))

GetColor <- function(value, GreenRed=FALSE, DisplayRange=3) {
    RGB <- function(x,y,z) rgb(x/255,y/255,z/255)
    missing <- is.na(value)
    good <- value[!missing]
    ans <- value
    if ( GreenRed )
        ans[missing] <- RGB(0, 0, 0)
    else
        ans[missing] <- RGB(255, 255, 255)
    tone <- abs(good) / DisplayRange * 255 + .5;
    tone[tone > 255] <- 255

    #classical: red and green on black background
    if (GreenRed)
        tone <- ifelse(good > 0, RGB(tone, 0, 0), RGB(0, 255-tone, 0))
    else
        tone <- ifelse(good > 0, RGB(255, 255 - tone, 255 - tone),
                       RGB(255 - tone, 255 - tone, 255) )
    ans[!missing]<-tone
    ans
}

#this is the wrong place but we need to put the methods stuff
#somewhere
.initClasses <- function(env) {
    setClass("uarray", representation(uexpr="character",
                                      samplenames="character",
                                      genenames="character",
                                      samplecluster="dendrogram",
                                      genecluster="dendrogram",
                                      phenodata="character") )

     #define a generic for obtaining the data
    setGeneric("uexpr", function(object) standardGeneric("uexpr"))

    setMethod("uexpr", "uarray", function(object) get(object@uexpr))

    #define a generic for obtaining the phenotypic data
    setGeneric("phenodata", function(object)
               standardGeneric("phenodata"))

    setMethod("phenodata", "uarray", function(object) {
        if( length(object@phenodata) == 1 )
            get(object@phenodata)
        else
            NULL
    })

 #deal with the names
    setGeneric("samplenames", function(object)
               standardGeneric("samplenames"))
    setMethod("samplenames", "uarray", function(object)
              object@samplenames)

    setGeneric("genenames", function(object)
               standardGeneric("genenames"))
    setMethod("genenames", "uarray", function(object) object@genenames )

                                        # deal with the clusters
    setGeneric("samplecluster", function(object)
               standardGeneric("samplecluster"))
    setMethod("samplecluster", "uarray", function(object) object@samplecluster )

    setGeneric("genecluster", function(object)
               standardGeneric("genecluster"))
    setMethod("genecluster", "uarray", function(object) object@genecluster )

                                        # plotting
    if( !isGeneric("plot") )
        setGeneric("plot")

    setMethod("plot", "uarray", function(x, ...) {
     expr <- as.matrix(uexpr(x))
     #scale
     expr <- sweep(expr, 1, apply(expr, 1, mean, na.rm = TRUE))
     f <- function(v) {
         v <- v[!is.na(v)]
         sqrt(sum(v^2)/max(1, length(v) - 1))
     }
     expr <- sweep(expr, 1, apply(expr, 1, f), "/")
     breaks <- seq(-3,3,by=.2)
     colors<- GetColor(breaks)
     breaks <- c(-100,breaks,100)
     colors <- c(colors[1], colors)
     opar<-par(mar=c(1,1,4,10))
     on.exit(par(mar=opar))
     image(1:ncol(expr), 1:nrow(expr), z = t(expr), axes = F,
           col=colors, breaks=breaks, xlab="", ylab="")
     axis(3, at=1:ncol(expr), labels=samplenames(x),tick=FALSE)
     axis(4, at=1:nrow(expr), labels=genenames(x), tick=FALSE, las=1)
 })

}
#given a set of labels, a set of weights and a set of clusters
#order the labels, by weights within clusters

order.restricted <- function(labels, weights, clusters) {
    n <- length(labels)
    if (length(weights) != n || length(clusters) != n)
        stop("all arguments must be the same length")
    which <- split(labels, clusters)
    wts <- split(weights, clusters)
    cwts <- sapply(wts, mean)
    cord <- order(cwts)
    rval <- NULL
    for(j in cord )
        rval <- c(rval, which[[j]][order(wts[[j]])])
    rval
}



##mstree, probably should get moved elsewhere
#ignore plane for now

mstree <- function(x, plane=T)
{
    require(mva)
    dmat <- dist(x)
    m <- length(dmat)
    n <- nrow(x)
    if( n < 2 )
        stop("you need at least two observations")
    Dlarge <- max(dmat) + 10
    z <- .Fortran("prtree",as.integer(n), as.integer(m),
                  A = integer(n), Dlarge, as.double(dmat),
                  B=integer(n), C=double(n), ifault=integer(1))

    return(list(mst=z$B, dist=z$C))
}

