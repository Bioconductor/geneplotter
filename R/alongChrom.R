.getCum <- function(eSet, chrom, specChrom) {

    # Extract the gene names of the chromosome of interest
    cLocs <- chromLocs(specChrom)
    genes <- cLocs[[chrom]]

    # Extract out of the expr set the genes that belong on this chrom
    usedGenes <- genes[names(genes) %in% geneNames(eSet)]
    chromExprs <- eSet@exprs[names(usedGenes),]

    # Order the columns by the ordering (by location) of the genes
    ord <- order(abs(usedGenes))
    chromExprs <- t(chromExprs[ord,])

    # Fill the matrix with the cumulative sum of the expression
    chromExprs <- apply(chromExprs, 1, cumsum)

    return(chromExprs)
}


alongChrom <- function(eSet, chrom, specChrom, index=TRUE, cumul=TRUE,
                     lTypes="1", colors="red", ...) {
    # Extract the gene names of the chromosome of interest
    cLocs <- chromLocs(specChrom)
    genes <- cLocs[[chrom]]

    # Extract out of the expr set the genes that belong on this chrom
    usedGenes <- genes[names(genes) %in% geneNames(eSet)]
    chromExprs <- eSet@exprs[names(usedGenes),]

    # Order the columns by the ordering (by location) of the genes
    ord <- order(abs(usedGenes))

    if (cumul == TRUE) {
        chromExprs <- t(chromExprs[ord,])
        ## Fill the matrix with the cumulative sum of the expression
        chromExprs <- apply(chromExprs, 1, cumsum)
    }
    else {
        chromExprs <- chromExprs[ord,]
    }

    ## Create the Y label based on if this is cumulative or not
    if (cumul == TRUE) {
        ylab <- "Cumulative expression levels"
    }
    else {
        ylab <- "Expression levels"
    }

    ## Plot data
    if (index == TRUE) {
        xPoints <- length(names(usedGenes)) - 1
        xlim <- xPoints

        ## Build main label
        main <- paste(ylab,"by genes in chromosome")
        main <- paste(main,chrom)

        matplot(0:xPoints, chromExprs, type="s", lty=lTypes, col=colors,
                xlab="", ylab=ylab, main=main, xaxt="n", cex.lab=0.9, ...)
        axis(1, at=c(0:xPoints), labels = names(usedGenes), las=2,
             cex.axis=0.7,)
    }
    else {
        xlim <- c(0, chromLengths(specChrom)[as.numeric(chrom)])
        xPoints <- as.numeric(abs(usedGenes)) + 1

        ## Build main label
        main <- paste(ylab,"in chromosome")
        main <- paste(main,chrom)
        main <- paste(main,"by relative position.")

        matplot(xPoints, chromExprs, type="s", lty=lTypes, col=colors,
                xlab="",ylab=ylab, xaxt="n", main=main, cex.lab=0.9, ...)
        axis(1, at=xPoints, labels = names(usedGenes), las=2,
             cex.axis=0.7,)
    }
}
