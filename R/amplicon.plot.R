##a function to get the chromosome order

make.chromOrd <- function(genome, gnames) {
    if(!is.character(genome) && length(genome != 1 ) )
        stop("need a character vector indicating the genome")
    require("annotate") || stop("need the annotate package")

    clname <- paste(genome, "chroloc", sep="")
    do.call("data", list(clname))
    allGcrloc <- mget(gnames, env=get(clname), ifnotfound=NA)
    myfun <- function(x) min(as.numeric(x))
    allGcloc <- sapply(allGcrloc, myfun)

    dname <- paste(genome, "chrom", sep="")
    if( !exists(dname, mode="environment") )
        do.call("data", list(dname))
    whichChrom <- unlist(mget(gnames, env=get(dname), ifnotfound=NA))
    byChr.cloc <- split(allGcloc, whichChrom)
    nchrom <- length(byChr.cloc)
    byChr.ord <- vector("list", length=nchrom)
    for(i in 1:nchrom ) byChr.ord[[i]] <- order(byChr.cloc[[i]])
    names(byChr.ord) <- names(byChr.cloc)
    byChr.ord$"NA" <- NULL
    byChr.ord
}

##actually do the amplicon plotting

amplicon.plot <- function(ESET, FUN, genome="hgu95A" ) {
    print("this will take a few seconds")
    tests <- esApply(ESET, 1, FUN)
    tests.pvals <- sapply(tests, function(x) x$p.value)
    tests.stats <- sapply(tests, function(x) x$statistic)

    dname <- paste(genome, "chrom", sep="")
    if( !exists(dname, mode="environment") )
        do.call("data", list(dname))

    whichChrom <- unlist(mget(geneNames(ESET), env=get(dname),
                                  ifnotfound=NA))
    ##split these by chromosome
    byChr.pv <- split(tests.pvals, whichChrom)
    byChr.stat <- split(tests.stats, whichChrom)

    byChr.pv$"NA" <- NULL
    byChr.stat$"NA" <- NULL

    chromOrd <- make.chromOrd(genome, geneNames(ESET))
    nchrom <- length(chromOrd)

    #get the names of the chromosome and their order
    #for plotting
    chromNames <- paste(genome, "chromNames", sep="")
    if( !exists(chromNames, mode="environment") )
        do.call("data", list(chromNames))
    geneOrd <- get(chromNames)

    chromOrd <- chromOrd[geneOrd]
    byChr.pv <- byChr.pv[geneOrd]
    byChr.stat <- byChr.stat[geneOrd]

    print("patience.....")
    chrlens <- sapply(chromOrd, length)

    collist <-  vector("list", length=nchrom)
    for(i in 1:nchrom) {
        smp <- ifelse(byChr.pv[[i]] < 0.05, 1, 0)
        dir <- byChr.stat[[i]]*smp
        cols <- ifelse(dir == 0 , 2, 3)
        cols <- ifelse(dir < 0, 1, cols)
        collist[[i]] <- cols[chromOrd[[i]]]
    }

    ncols <- vector("list", length=nchrom)
    maxlen <- max(chrlens)
    for(i in 1:nchrom) {
        extras<- maxlen - chrlens[i]
        ncols[[i]]<- c(collist[[i]], rep(2, extras))
    }
    z<- data.frame(ncols)
    z<- as.matrix(z)
    image(1:maxlen, 1:nchrom, z, col=c("blue","white", "red"),
    xlab="Gene location", ylab="Chromosome", axes=FALSE )
    axis(2, at = 1:nchrom, labels=names(byChr.pv))
}

