.getGenes <- function(eSet, chrom, specChrom) {
    ## Will acquire the set of genes used in alongChrom

    ## Extract the gene names of the chromosome of interest
    cLocs <- chromLocs(specChrom)
    genes <- cLocs[[chrom]]

    ## Extract out of the expr set the genes that belong on this chrom
    usedGenes <- genes[names(genes) %in% geneNames(eSet)]

    ## Order the genes by location
    usedGenes <- sort(abs(usedGenes))

    return(usedGenes)
}

.getExprs <- function(eSet, usedGenes,
                      plotFormat=c("cumulative","local")[1])
{
    ## Will get the expression data for the given genes out of the
    ## expr set.  If plotFormat is set to cumulative, will generate the
    ## cumulative sum of this data across the genes.

    ## Split out only the genes on the desired chrom from the exprset
    chromExprs <- eSet@exprs[names(usedGenes),]

    if (plotFormat == "cumulative") {
       chromExprs <- t(chromExprs)
       ## Fill the matrix with the cumulative sum of the expression
       chromExprs <- apply(chromExprs, 1, cumsum)
   }

   return(chromExprs)
}

identifyLines <- function(identEnvir, ...) {
    ## Will call identify() on teh alongChrom() plot to detail which
    ## lines map tow which samples

    points <- multiget(c("X","Y"), envir=identEnvir)

    xPoints <- points$X
    yPoints <- points$Y

    x <- identify(rep(xPoints,ncol(yPoints)), yPoints,
                  labels=rep(colnames(yPoints),
                  rep(nrow(yPoints),ncol(yPoints))), ...)

    return(x)
}

.scaleData <-
    function(chromData,
    method=c("none","zscale","rangescale","rankscale", "zrobustscale")[1])
{
    ## Will scale the data set to be plotted based on a variety of
    ## methods

    if (method != "none") {
        for (i in 1:nrow(chromData)) {
            x <- chromData[i,]
            if (method == "zscale") {
                chromData[i,] <- (x - mean(x))/sd(x)
            }
            else if (method == "rangescale") {
                curRange <- range(x)
                chromData[i,] <- (x - curRange[1])/(curRange[2] - curRange[1])
            }
            else if (method == "rankscale") {
                chromData[i,] <- rank(x)
            }
            else if (method == "zrobustscale") {
                chromData[i,] <- (x - median(x))/mad(x)
            }
        }
    }

    return(chromData)
}

.cullXPoints <- function(xPoints) {
    ## Will reduce the xPoints vector to a visibly manageable size
    ## Currently if the size > 40, will leave every Nth point where
    ## xPoints/40 = N.

    maxSize <- 40

    if (length(xPoints) > maxSize) {
        ## Calculate N, and then get the maxSize elements from every
        ## Nth element.  Problem: Sometiems will generate a few extra
        ## due to integer division on N.
        N <- length(xPoints) %/% maxSize

        ## Start from 2 for now as a hack to keep from getting 0th
        ## entity, which throws off the labeling.
        keep <- seq(1,length(xPoints),N)

        xPoints <- xPoints[keep]
    }

    return(xPoints)
}

alongChrom <- function(eSet, chrom, specChrom,
                       xloc=c("equispaced", "physical")[1],
                       plotFormat=c("cumulative", "local")[1],
                       scale=c("none","zscale","rankscale","rangescale",
                               "zrobustscale")[1],
                       lTypes="1", colors="red", ...) {

    ## Will plot a set of exprset samples by genes of a chromosome
    ## according to their expression levels.
    ## Get the genes to display
    usedGenes <- .getGenes(eSet, chrom, specChrom)

    ## Get the expression data, cumulative or otherwise
    chromExprs <- .getExprs(eSet, usedGenes, plotFormat)

    ## The Y axis label varies according to if we're taking
    ## cummulative sums or not
    if (plotFormat == "cumulative") {
        ylab <- "Cumulative expression levels"
    }
    else {
        ylab <- "Expression levels"
    }

    chromExprs <- .scaleData(chromExprs,scale)

    ## Plot data
    if (xloc == "equispaced") {
        xPoints <- length(names(usedGenes)) - 1
        xPoints <- 0:xPoints

        ## Build main label
        main <- paste(ylab,"by genes in chromosome")
        main <- paste(main,chrom)
        main <- paste(main,", scaling method:",sep="")
        main <- paste(main,scale)
    }
    else {
        xPoints <- as.numeric(abs(usedGenes)) + 1
        ## Build main label
        main <- paste(ylab,"in chromosome")
        main <- paste(main,chrom)
        main <- paste(main,"by relative position, scaling method:")
        main <- paste(main,scale)
    }

    ## Make sure that xPoints isn't exceeding our visual maximum.
    ## If so, reduce the number of poitns to actually be displayed.
    dispXPoints <- .cullXPoints(xPoints)

    ## Plot the graph
    matplot(xPoints, chromExprs, type="S", lty=lTypes, col=colors,
            xlab="",ylab=ylab, xaxt="n", main=main, cex.lab=0.9, ...)
    axis(1, at=dispXPoints, labels = names(usedGenes)[dispXPoints+1], las=2,
         cex.axis=0.7,)

    ## Create an environment that contains the necessary X & Y points
    ## for use with identify()
    identEnv <- new.env()
##    multiassign(c("X","Y"),c(xPoints,chromExprs),envir=identEnv)
    assign("X", xPoints, envir=identEnv)
    assign("Y", chromExprs, envir=identEnv)
    return(identEnv)
}

