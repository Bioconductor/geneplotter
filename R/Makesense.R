if (is.null(getGeneric("Makesense")))
    setGeneric("Makesense", function(object, lib, f)
               standardGeneric("Makesense"))

setMethod("Makesense", "exprSet", function(object, lib, f) {
    if (missing(lib))
        lib <- annotation(object)

    if (length(lib) == 0)
        stop("No annotation library supplied")

    if (missing(f))
        f <- 1/10
    Makesense(exprs(object), lib, f)
})

setMethod("Makesense", "matrix", function(object, lib, f) {
    if(!require(lib, character.only=TRUE))
        stop("data library ",lib," not avaiable")

    if (missing(f))
        f <- 1/10

    libCHR <- paste(lib, "CHR", sep="")
    libCHRLOC <- paste(lib, "CHRLOC", sep="")
    ##find the positions
    chr <- contents(get(libCHR))
    ##for which do we have one and only chromosome no.
    chrSet <- c(1:22, "X", "Y")
    oneC <- sapply(chr, function(x) if(x %in% chrSet &&
                                       length(x)==1) TRUE else FALSE)
    ##have loc
    chrL <- contents(get(libCHRLOC))
    oneL <- sapply(chrL, function(x) if( !is.na(x) &&
                                        length(x)==1) TRUE else FALSE)
    chrName <- unlist(chr[oneC&oneL])
    chrPos <- unlist(chrL[oneC&oneL])

    cP <- split(chrPos, chrName)

    gE <- object[names(oneC&oneL)[oneC&oneL],]
    ans2 <- vector("list", length=ncol(gE))
    libCHRLENGTHS <- get(paste(lib,"CHRLENGTHS",sep=""))

    for( j in 1:ncol(gE) ) {
        s1 <- split(gE[,j], chrName)
        ans<-NULL
        for(i in names(cP) ) {
            d1 <- s1[[i]]
            cL <- cP[[i]]
            dp <- d1[cL>0]
            lp <- cL[cL>0]
            dn <- d1[cL<0]
            ln <- -cL[cL<0]
            lw1 <- lowess(lp, dp, f=f)
            names(lw1$x) <- names(dp)[order(lp)]
            lw2 <- lowess(ln, dn, f=f)
            names(lw2$x) <- names(dn)[order(ln)]
            lw2$y <- -lw2$y
            ans[[i]] <- list(posS = lw1, negS =lw2)
        }
        ans2[[j]] <- ans
    }
    list(ans2=ans2, lib=lib)
})
