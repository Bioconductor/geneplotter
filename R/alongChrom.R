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
                     lTypes="1", colors="red") {
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


    ## Plot data
    if (index == TRUE) {
        xPoints <- length(names(usedGenes)) - 1
        xlim <- xPoints
        matplot(0:xPoints, chromExprs, type="s", lty=lTypes, col=colors,
                xlab="genes", ylab="Cummulative Expression", xaxt="n",
                cex.lab=0.9)
        axis(1, at=c(0:xPoints), labels = names(usedGenes), las=2,
             cex.axis=0.7,)
    }
    else {
        xlim <- c(0, chromLengths(specChrom)[as.numeric(chrom)])
        xPoints <- as.numeric(abs(usedGenes)) + 1
        matplot(xPoints, chromExprs, type="s", lty=lTypes, col=colors,
                xlab="genes", ylab="Cummulative Expression", xaxt="n",
                cex.lab=0.9)
        axis(1, at=xPoints, labels = names(usedGenes), las=2,
             cex.axis=0.7,)
    }


}
