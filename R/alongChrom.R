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

.identifyLines <- function(xPoints, yPoints) {
    ## Will call identify() on teh alongChrom() plot to detail which
    ## lines map tow which samples

    identify(rep(xPoints,ncol(yPoints)), yPoints,
             labels=rep(colnames(yPoints), rep(nrow(yPoints),ncol(yPoints))))
}

alongChrom <- function(eSet, chrom, specChrom,
                       xloc=c("equispaced", "physical")[1],
                       plotFormat=c("cumulative", "local")[1],
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

    ## Plot data
    if (xloc == "equispaced") {
        xPoints <- length(names(usedGenes)) - 1
        xPoints <- 0:xPoints

        ## Build main label
        main <- paste(ylab,"by genes in chromosome")
        main <- paste(main,chrom)
    }
    else {
        xPoints <- as.numeric(abs(usedGenes)) + 1
        ## Build main label
        main <- paste(ylab,"in chromosome")
        main <- paste(main,chrom)
        main <- paste(main,"by relative position.")
    }

    ## Plot the graph
    matplot(xPoints, chromExprs, type="S", lty=lTypes, col=colors,
            xlab="",ylab=ylab, xaxt="n", main=main, cex.lab=0.9, ...)
    axis(1, at=xPoints, labels = names(usedGenes), las=2,
         cex.axis=0.7,)

    ## Call identify() on the plot
    .identifyLines(xPoints, chromExprs)
}

