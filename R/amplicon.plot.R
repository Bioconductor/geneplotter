##a function to get the chromosome order

make.chromOrd <- function(genome) {
    if(!is.character(genome) && length(genome != 1 ) )
        stop("need a character vector indicating the genome")
    require("annotate") || stop("need the annotate package")

    hgu95cloc <- read.annotation(paste(genome,"chroloc", sep=""))
    allGcrloc <- multiget(gnames, env=hgu95cloc)
    myfun <- function(x) min(as.numeric(x))
    allGcloc <- sapply(allGcrloc, myfun)

    if( !exists("hgu95Chrom", mode="environment") )
        hgu95Chrom <- read.annotation(paste(genome,"chrom", sep=""))
    whichChrom <- unlist(multiget(gnames, env=hgu95Chrom))
    byChr.cloc <- split(allGcloc, whichChrom)
    byChr.ord <- vector("list", length=25)
    for(i in 1:25 ) byChr.ord[[i]] <- order(byChr.cloc[[i]])
    names(byChr.ord) <- names(byChr.cloc)
    byChr.ord$"NA" <- NULL
    byChr.ord
}

##actually do the amplicon plotting

amplicon.plot <- function(ESET, FUN, genome ) {
    print("this will take a few seconds")
    tests <- esApply(ESET, 1, FUN)
    tests.pvals <- sapply(tests, function(x) x$p.value)
    tests.stats <- sapply(tests, function(x) x$statistic)

    if( !exists("hgu95Chrom", mode="environment") )
        hgu95Chrom <- read.annotation(paste(genome,"chrom", sep=""))

    whichChrom <- unlist(multiget(gnames, env=hgu95Chrom))
    ##split these by chromosome
    byChr.pv <- split(igenes.pvals, whichChrom)
    byChr.stat <- split(igenes.stats, whichChrom)

    byChr.pv$"NA" <- NULL
    byChr.stat$"NA" <- NULL

    chromOrd <- make.chromOrd(genome)
    nchrom <- length(chromOrd)
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
    image(z, col=c("blue","white", "red"), xlab="Gene location",
          ylab="Chromosome", axes=FALSE )
    axis(2, at = (0:(nchrom-1)/nchrom, labels=names(byChr.pv))
}

