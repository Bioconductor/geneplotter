alongChrom <- function(eSet, chrom, specChrom, xlim, whichGenes,
                       xloc=c("equispaced", "physical"),
                       plotFormat=c("cumulative", "local","image"),
                       scale=c("none","zscale","rankscale","rangescale",
                               "zrobustscale"),
                       geneSymbols=FALSE,
                       lty=1, type="S", colors="red", ...) {

    ## Will plot a set of exprset samples by genes of a chromosome
    ## according to their expression levels.

    ## Get the genes to display
    usedGenes <- usedChromGenes(eSet, chrom, specChrom)

    ## Filter out any NA positioned genes
    usedGenes <- usedGenes[!is.na(usedGenes)]

    ## Limit genes to requested range
    if (!missing(xlim)) {
        usedGenes <- .limitXRange(xlim, usedGenes)
    }

    geneNames <- names(usedGenes)
    if (geneSymbols == TRUE) {
        geneNames <- .getGeneSyms(geneNames)
    }


    ##make sure we get the full name for all args
    xloc <- match.arg(xloc)
    plotFormat <- match.arg(plotFormat)
    scale <- match.arg(scale)

    ## Select out requested genes
    if (!missing(whichGenes)) {
        nameLocs <- geneNames %in% whichGenes
        if (!all(nameLocs)) {
            print("Warning: Not all requested genes are displayed.")
        }
        usedGenes <- usedGenes[nameLocs]
        geneNames <- names(usedGenes)
    }

    strands <- ifelse(usedGenes>0,"+","-")

    ## Check for duplicated positions
    dup <- duplicated(abs(as.numeric(usedGenes)))
    dup <- which(dup)

    nGenes <- length(usedGenes)
    if (nGenes == 0) {
        .emptyACPlot(chrom)
        return()
    }
    else if (nGenes == 1) {
        ## !!!! TODO: Plot the single value as is instead of this
        x <- paste("Only gene to be plotted: ",
                   geneNames,":",as.numeric(abs(usedGenes)),sep="")
        stop(x)
    }

    ## Get the expression data, cumulative or otherwise
    chromExprs <- .getExprs(eSet, usedGenes, plotFormat,scale)

    ## Create labels for plotting
    if (plotFormat == "cumulative") {
        ylab <- "Cumulative expression levels"
    }
    else {
        ylab <- "Expression levels"
    }
    xlab <- "Representative Genes"
    main <- .buildMainLabel(ylab, chrom, xloc, scale)

    ## If image plot was requested, split off here
    if (plotFormat == "image") {
        return(.doImagePlot(chromExprs, chrom, geneNames, strands, scale, main,
                            xlab, 10))
    }

    ## Define the points for the x axis
    if (xloc == "equispaced") {
        xPoints <- length(geneNames) - 1
        xPoints <- 0:xPoints
    }
    else if (xloc == "physical") {
        xPoints <- abs(as.numeric(usedGenes)) + 1
        if (any(dup)) {
            xPoints[dup] <- xPoints[dup] +
                (0.5 * (xPoints[dup+1]-xPoints[dup]))
        }
    }

    ## Local plots are shifted over, so create a faxe xPoint on the end
    if (plotFormat == "local") {
        chromExprs <- rbind(chromExprs,chromExprs[length(xPoints),])
        xPoints <- c(xPoints,xPoints[length(xPoints)]+1)
        xPool <- xPoints[1:length(xPoints)-1]
    }
    else {
        xPool <- xPoints
    }


    ## Plot the graph
    opar <- par(mar=c(6,5,4,1),mgp=c(4,1,0))
    on.exit(par(opar))
    matplot(xPoints, chromExprs, type=type, lty=lty, col=colors,
            xlab=xlab,ylab=ylab, xaxt="n", main=main, cex.lab=0.9,...)
    .dispXAxis(xPoints, xPool, geneNames, plotFormat,strands)
    if (any(dup)) {
        if (xloc == "equispaced") {
            segments(dup-2,1,dup-1,1,col="cyan",lwd=2)
        }
        else {
            segments(xPoints[dup-1],1,xPoints[dup],1,col="cyan",lwd=2)
        }
    }

    ## Create an environment that contains the necessary X & Y points
    ## for use with identify()
    identEnv <- new.env()
    multiassign(c("X","Y"),list(xPoints,chromExprs),envir=identEnv)

    return(identEnv)
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

.dispXAxis <- function(xPoints, xPool, geneNames, plotFormat, strands) {

    ## Make sure that xPoints isn't exceeding our visual maximum.
    ## If so, reduce the number of poitns to actually be displayed.
    dispXPoints <- .cullXPoints(xPool)
    dispPointLocs <- match(dispXPoints,xPoints)

    if (plotFormat == "local") {
        dispXPoints <- dispXPoints+0.5
    }

    axis(1, at=dispXPoints, labels = geneNames[dispPointLocs], las=2,
         cex.axis=0.7,)

    axis(3, at=dispXPoints, labels = strands[dispPointLocs], cex.axis=0.8)
}

.doImagePlot <- function(exprs,chrom, geneNames, strands, scale, main,
                         xlab, nCols) {
    ## Passed in the expression matrix, the names of the
    ## used genes, the name of the chromosome, the scaling method & the number
    ## of colours to utilize in the plot, will generate
    ## an image plot
    ngenes <- nrow(exprs)
    nsamp <- ncol(exprs)

    ## Get the colour mapping
    d <- dChip.colors(nCols)
    w <- sort(exprs)
    b <- quantile(w,probs=seq(0,1,(1/length(d))))

    ## Build the labels
    main <- main
    xlab <- xlab
    ylab="Samples"

    ## Build the plot
    xPoints <- 1:ngenes

    image(x=xPoints,y=1:(nsamp+1),z=exprs, col=d, breaks=b,
          xlab=xlab, ylab=ylab, main=main, axes=FALSE)
    axis(2, at=1:nsamp, labels=colnames(exprs))

    .dispXAxis(xPoints, xPoints, geneNames, "image", strands)

    ## Create an environment that contains the necessary X & Y points
    ## for use with identify()

    ## !!! As is, does not return the proper data
    ##    identEnv <- new.env()
    ##    multiassign(c("X","Y"),list(xPoints,exprs),envir=identEnv)
    ##    return(identEnv)
}

.limitXRange <- function(xlim, usedGenes) {

    if (!missing(xlim)) {
        if (length(xlim) == 2) {
            if (is.character(xlim)) {
                ## If a pair of gene names are provided, get hteir
                ## locations, and then use them as the xlim values.
                xlim[1] <- as.numeric(usedGenes[xlim[1]])
                xlim[2] <- as.numeric(usedGenes[xlim[2]])
                if ((is.na(xlim[1]))|(is.na(xlim[2]))) {
                    print("Error: Bad xlim parameters provided.")
                    xlim[1] = 0
                    xlim[2] = 0
                    usedGenes <- NULL
                }
                ## Place them in proper numerical order
                xlim <- xlim[order(xlim)]
            }
            ## At this point, we're dealing with a pair of numerical
            ## values to denote the location range (in base pairs).
            ## Ensure that the max is > than the min, then pick out
            ## the remaining genes
            if (xlim[2] > xlim[1]) {
                lowLim <- match(xlim[1],usedGenes)
                if (is.na(lowLim)) {
                    lowLim <- .getClosestPos(xlim[1],usedGenes)
                }

                hiLim <- match(xlim[2], usedGenes)
                if (is.na(hiLim)) {
                    hiLim <- .getClosestPos(xlim[2],usedGenes)
                }

                subs <- seq(lowLim,hiLim)
                usedGenes <- usedGenes[subs]
            }
            else {
                print("Error: Bad xlim parameters provided.")
                usedGenes <- NULL
            }
        }
        else {
            print("Error: Bad xlim parameters provided.")
            usedGenes <- NULL
        }
    }

    return(usedGenes)
}

.getGeneSyms <- function(affys) {
    data(hgu95Asym)
    syms <- multiget(affys, env=hgu95Asym)
    syms[is.na(syms)] <- affys[is.na(syms)]
    return(as.character(syms))
}

.getClosestPos <- function(val, usedGenes) {
    ## Given a value, finds the closest value in usedGenes to the
    ## passed value and returns its location in the usedGenes vector

    dists <- abs(val-abs(as.numeric(usedGenes)))
    closest <- match(min(dists), dists)
    return(closest)
}

.scaleData <-
    function(chromData,
    method=c("none","zscale","rangescale","rankscale", "zrobustscale"))
{
    ## Will scale the data set to be plotted based on a variety of
    ## methods

    method <- match.arg(method)
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
            else {
                stmt <- paste("method:", method, ", is not implemented yet")
                stop(stmt)
            }
        }
    }

    return(chromData)
}

.cullXPoints <- function(xPoints) {
    ## Will reduce the xPoints vector to a visibly manageable size
    ## Currently if the size > 40, will leave every Nth point where
    ## xPoints/maxSize = N.  Maximum number of points is determined
    ## by determining the size of the label text and filling up 65%
    ## of the axis space with labels.

    ## First get the size of the plotting region
    preg <- par('pin')[1] * 0.65
    ## Now get the font size
    strsize <- strheight("test",units="inches")
    ## Calculate the maxSize
    maxSize <- preg %/% strsize

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

.buildMainLabel <- function(ylab, chrom, xloc, scale) {
    if (xloc == "relative") {
        main <- paste(ylab, "in chromosome", chrom,
                      "by relative position\n")
    }
    else {
        main <- paste(ylab, "by genes in chromosome", chrom, "\n")
    }

    main <- paste(main,"scaling method:",scale,"\n")

    return(main)
}

.emptyACPlot <- function(chrom) {
    plot.new()
    axis(1,labels=rep("NA",6))
    axis(2, labels=rep("NA",6))
    main <- paste("Plot empty, no genes from chromosome",chrom,
                  "in exprSet provided.\n")

    title(main = main)
}

.getExprs <- function(eSet, usedGenes,
                      plotFormat=c("cumulative","local", "image"),
                      scale=c("none","zscale","rangescale","rankscale", "zrobustscale"))
{
    ## Will get the expression data for the given genes out of the
    ## expr set.  If plotFormat is set to cumulative, will generate the
    ## cumulative sum of this data across the genes.

    ## Split out only the genes on the desired chrom from the exprset
    plotFormat <- match.arg(plotFormat)
    scale <- match.arg(scale)

    chromExprs <- eSet@exprs[names(usedGenes),]

    chromExprs <- .scaleData(chromExprs,scale)

    if (plotFormat == "cumulative") {
        chromExprs <- t(chromExprs)
        ## Fill the matrix with the cumulative sum of the expression
        chromExprs <- apply(chromExprs, 1, cumsum)
    }

    return(chromExprs)
}
