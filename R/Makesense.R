setGeneric("Makesense", function(expr, lib, ...) standardGeneric("Makesense"))

setMethod("Makesense", signature(expr="ExpressionSet", lib="character"),
          function(expr, lib, f=1/10) {
              Makesense(exprs(expr), lib, f)
})

setMethod("Makesense", signature(expr="ExpressionSet", lib="missing"),
          function(expr, f=1/10) {
              Makesense(expr, annotation(expr), f)
          })

setMethod("Makesense", signature(expr="matrix", lib="character"),
          function(expr, lib, f=1/10) {
    if (length(lib) != 1 || nchar(lib) < 1)
      stop("'lib' argument must be length one")

    genes <- rownames(expr)
    libCHR <- getAnnMap("CHR", lib)
    libCHRLOC <- getAnnMap("CHRLOC", lib)
    ## Select genes that are annotated at exactly _one_ chromosome.
    chr <- mget(genes, envir=libCHR, ifnotfound=NA)
    oneC <- sapply(chr, function(x)
                   if (length(x) == 1 && !is.na(x)) TRUE else FALSE)
    ## Select genes that are annotated at exactly _one_ chrom location
    ##
    ## FIXME: There are many genes with multiple CHRLOC entries, is
    ##        there anything we can do to keep more of them?
    chrL <- mget(genes, envir=libCHRLOC, ifnotfound=NA)
    oneL <- sapply(chrL, function(x)
                   if (length(x) == 1 && !is.na(x)) TRUE else FALSE)
    wanted <- oneC & oneL
    chrName <- unlist(chr[wanted])
    chrPos <- unlist(chrL[wanted])

    cP <- split(chrPos, chrName)

    gE <- expr[wanted, ]
    ans2 <- vector("list", length=ncol(gE))

    for( j in 1:ncol(gE) ) {
        s1 <- split(gE[,j], chrName)
        ans <- NULL
        for (i in names(cP)) {
            d1 <- s1[[i]]
            cL <- cP[[i]]
            dp <- d1[cL>0]
            lp <- cL[cL>0]
            dn <- d1[cL<0]
            ln <- -cL[cL<0]
            if (length(lp)) {
                lw1 <- lowess(lp, dp, f=f)
                names(lw1$x) <- names(dp)[order(lp)]
            } else {
                lw2 <- list(x=numeric(0), y=numeric(0))
            }
            if (length(ln)) {
                lw2 <- lowess(ln, dn, f=f)
                names(lw2$x) <- names(dn)[order(ln)]
                lw2$y <- -lw2$y
            } else {
                lw2 <- list(x=numeric(0), y=numeric(0))
            }
            ans[[i]] <- list(posS = lw1, negS =lw2)
        }
        ans2[[j]] <- ans
    }
    list(ans2=ans2, lib=lib)
})
