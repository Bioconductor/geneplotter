alongChrom <- function(eSet, chrom, specChrom, xlim, whichGenes,
                       plotFormat=c("cumulative", "local","image"),
                       xloc=c("equispaced", "physical"),
                       scale=c("none","zscale","rankscale","rangescale",
                               "zrobustscale"),
                       geneSymbols=FALSE, byStrand=FALSE,
                       colors="red",  lty=1, type="S", ...) {

    ## Will plot a set of exprset samples by genes of a chromosome
    ## according to their expression levels.

    ##make sure we get the full name for all args
    xloc <- match.arg(xloc)
    plotFormat <- match.arg(plotFormat)
    scale <- match.arg(scale)

    ## Get plotting labels
    labEnv <- getACPlotLabs(plotFormat, chrom, xloc, scale)

    ## Get the genes to display
    usedGenes <- usedChromGenes(eSet, chrom, specChrom)
    ## Filter out any NA positioned genes
    usedGenes <- usedGenes[!is.na(usedGenes)]
    ## Limit genes to requested range
    if (!missing(xlim)) {
        usedGenes <- limitACXRange(xlim, usedGenes)
    }
    geneNames <- names(usedGenes)
    if (geneSymbols == TRUE) {
        geneNames <- getACGeneSyms(geneNames, specChrom)
    }
    ## Select out requested genes
    if (!missing(whichGenes)) {
        nameLocs <- geneNames %in% whichGenes
        if (!all(nameLocs)) {
            print("Warning: Not all requested genes are displayed.")
        }
        usedGenes <- usedGenes[nameLocs]
        geneNames <- names(usedGenes)
    }

    ## Handle cases where we have filter out all but 0 or 1 gene.
    nGenes <- length(usedGenes)
    if (nGenes == 0) {
        emptyACPlot(chrom)
        return()
    }
    else if (nGenes == 1) {
        ## !!!! TODO: Plot the single value as is instead of this
        x <- paste("Only gene to be plotted: ",
                   geneNames,":",as.numeric(abs(usedGenes)),sep="")
        stop(x)
    }

    ## Get the expression data, cumulative or otherwise
    chromExprs <- getACExprs(eSet, usedGenes, plotFormat,scale)

    ## Figure out which strands each gene is on
    strands <- ifelse(usedGenes>0,"+","-")

    ## Check for duplicated positions
    dup <- which(duplicated(abs(as.numeric(usedGenes))))
    dup <- dup[!is.na(dup)]

    dataEnv <- getACDataEnv(chromExprs, geneNames, strands,
                            byStrand, dup)

    ## If image plot was requested, split off here
    switch(plotFormat,
           "image" = return(doACImagePlot(dataEnv, labEnv, 10)),
           "local" = return(doACLocalPlot(dataEnv, labEnv, colors)),
           "cumulative" = return(doACCumPlot(dataEnv, labEnv,
                           usedGenes, xloc, colors, lty, type, ...))
           )
}

doACImagePlot <- function(dataEnv, labEnv, nCols) {
    ## Passed in the expression matrix, the names of the
    ## used genes, the name of the chromosome, the scaling method & the number
    ## of colours to utilize in the plot, will generate
    ## an image plot
    chromExprs <- dataEnv$chromExprs
    byStrand <- dataEnv$byStrand

    ngenes <- nrow(chromExprs)
    nsamp <- ncol(chromExprs)

    ## Get the colour mapping
    d <- dChip.colors(nCols)
    w <- sort(chromExprs)
    b <- quantile(w,probs=seq(0,1,(1/length(d))))

    ## retrieve the labels
    xlab <- labEnv$xlab
    ylab <- labEnv$ylab
    main <- labEnv$main

    ## Build the plot
    xPoints <- 1:ngenes

    if (byStrand==TRUE) {
        strands <- dataEnv$strands

        mfPar <- par(mfrow = c(2,1))
        on.exit(par(mfPar))
        midVal <- b[length(b)/2]
        pos <- xPoints[which(strands == "+")]
        neg <- xPoints[which(strands == "-")]
        posExprs <- chromExprs
        posExprs[neg,] <- midVal
        negExprs <- chromExprs
        negExprs[pos,] <- midVal

        image(x=xPoints,y=1:(nsamp+1),z=posExprs, col=d, breaks=b,
              xlab=xlab, ylab=ylab, main=main, axes=FALSE)
        axis(2, at=1:nsamp, labels=colnames(posExprs))
        dispACXaxis(xPoints, dataEnv, "image")
        mtext("Plus",
              side=3,line=0.35,outer=FALSE,at=mean(par("usr")[1:2]))
        image(x=xPoints,y=1:(nsamp+1),z=negExprs, col=d, breaks=b,
              xlab=xlab, ylab=ylab, axes=FALSE)
        axis(2, at=1:nsamp, labels=colnames(chromExprs))
       dispACXaxis(xPoints, dataEnv, "image")
        mtext("Minus",
              side=3,line=0.35,outer=FALSE,at=mean(par("usr")[1:2]))
    }
    else {
        image(x=xPoints,y=1:(nsamp+1),z=chromExprs, col=d, breaks=b,
              xlab=xlab, ylab=ylab, main=main, axes=FALSE)
        axis(2, at=1:nsamp, labels=colnames(chromExprs))
        dispACXaxis(xPoints, dataEnv, "image")
    }
    invisible(chromExprs)
}

doACMatPlot <- function(xPoints, dataEnv, xlim, ylim, type, lty, col,
                       labEnv, xloc, ...) {
    xlab <- labEnv$xlab
    ylab <- labEnv$ylab
    main <- labEnv$main

    chromExprs <- dataEnv$chromExprs

    matplot(xPoints, chromExprs, xlim=xlim, ylim=ylim,type=type,
            lty=lty, col=col, xlab=xlab,ylab=ylab, main=main,
            xaxt="n", cex.lab=0.9,...)

    dispACXaxis(xPoints, dataEnv, xloc, "cumulative")
}

doACLocalPlot <- function(dataEnv, labEnv, colors) {
    ## retrieve the labels
    xlab <- labEnv$xlab
    ylab <- labEnv$ylab
    main <- labEnv$main

    envTitles <- c("chromExprs", "geneNames", "strands", "dup")
    ## Retrieve data values
    envVals <- mget(c(envTitles,"byStrand"),envir=dataEnv, ifnotfound=NA)

    xPoints <- 1:nrow(envVals$chromExprs)

    if (envVals$byStrand == TRUE) {
        mfPar <- par(mfrow = c(2,1))
        on.exit(par(mfPar),add=TRUE)
        strandVals <- getACStrandVals(envVals$chromExprs,
                                      envVals$strands, xPoints,
                                      envVals$dup, envVals$geneNames,
                                      "local")
        multiassign(envTitles,list(strandVals$posExprs,
                                      strandVals$posGen,
                                      strandVals$posStr,
                                      strandVals$posDup),envir=dataEnv)
        z <- boxplot(data.frame(t(strandVals$posExprs)), plot=FALSE)
        z$stats[,strandVals$nts] <- NA
        bxp(z,col=colors, xaxt="n", xlab=xlab, ylab=ylab, main=main,
            cex.lab=0.9)
        mtext("Plus", side=3,line=0.35,outer=FALSE,
              at=mean(par("usr")[1:2]))
        dispACXaxis(strandVals$posPoints, dataEnv)
        ## Now do negative
        multiassign(envTitles,list(strandVals$negExprs,
                                      strandVals$negGen,
                                      strandVals$negStr,
                                      strandVals$negDup),envir=dataEnv)
        z <- boxplot(data.frame(t(strandVals$negExprs)), plot=FALSE)
        z$stats[,strandVals$pts] <- NA
        bxp(z,col=colors, xaxt="n", xlab=xlab, ylab=ylab, main=main,
            cex.lab=0.9)
        mtext("Minus", side=3,line=0.35,outer=FALSE,
              at=mean(par("usr")[1:2]))
        dispACXaxis(strandVals$negPoints, dataEnv)
    }
    else {
        boxplot(data.frame(t(envVals$chromExprs)), col=colors, xlab=xlab,
                ylab=ylab, main=main, cex.lab=0.9, xaxt="n")
        dispACXaxis(xPoints, dataEnv)
    }
    invisible(envVals$chromExprs)
}

doACCumPlot <- function(dataEnv, labEnv, usedGenes, xloc, colors, lty, type,
                         ...) {
    envTitles <- c("chromExprs", "dup", "geneNames", "strands",
                   "byStrand")
    envVals <- mget(envTitles, envir=dataEnv, ifnotfound=NA)

    ## Create a fictitious start & end gene to help with plots
    start <- abs(as.numeric(usedGenes[1])) * 0.8
    end <- abs(as.numeric(usedGenes[length(usedGenes)])) * 1.2
    usedGenes <- c(start,usedGenes,end)

    geneNames <- envVals$geneNames <- c("",envVals$geneNames,"")
    strands <- envVals$strands <- c("",envVals$strands,"")
    ## Also need to give them data in the chromExprs matrix
    ## just copy data from the one next to them.
    chromExprs <- envVals$chromExprs
    chromExprs <- envVals$chromExprs <- rbind(chromExprs[1,],chromExprs,
                                              chromExprs[nrow(chromExprs),])
    dup <- envVals$dup <- envVals$dup + 1

    multiassign(envTitles, envVals, envir=dataEnv)

    ## Define the points for the X axis
    if (xloc == "equispaced")
        xPoints <- 1:length(usedGenes)
    else if (xloc == "physical") {
        xPoints <- abs(as.numeric(usedGenes)) + 1
        xPoints <- fixACPhysPoints(xPoints, dup)
    }

    ## Get x & y ranges
    xlim <- range(xPoints)
    ylim <- range(chromExprs)
    ylim[1] <- ylim[1]-0.1

    ## Plot the graph
    opar <- par(mar=c(6,5,4,1),mgp=c(4,1,0))
    on.exit(par(opar),add=TRUE)

    if (envVals$byStrand == TRUE) {
        mfPar <- par(mfrow = c(2,1))
        on.exit(par(mfPar),add=TRUE)

        strandVals <- getACStrandVals(chromExprs, strands, xPoints, dup,
                                    geneNames, "cumulative", xloc)

        strandTitles <- c("chromExprs", "geneNames","strands", "dup")
        multiassign(strandTitles,list(strandVals$posExprs,
                                      strandVals$posGen,
                                      strandVals$posStr,
                                      strandVals$posDup),envir=dataEnv)
        doACMatPlot(strandVals$posPoints, dataEnv, xlim=xlim, ylim=ylim,
                   type=type, lty=lty, col=colors,
                   labEnv=labEnv, xloc=xloc, ...)
        mtext("Plus", side=3,line=0.35,outer=FALSE,
              at=mean(par("usr")[1:2]))

        multiassign(strandTitles,list(strandVals$negExprs,
                                      strandVals$negGen,
                                      strandVals$negStr,
                                      strandVals$negDup),envir=dataEnv)
        doACMatPlot(strandVals$negPoints, dataEnv, xlim=xlim, ylim=ylim,
                   type=type, lty=lty, col=colors, labEnv=labEnv,
                   xloc=xloc, ...)
        mtext("Minus", side=3,line=0.35,outer=FALSE,
              at=mean(par("usr")[1:2]))
    }
    else {
        doACMatPlot(xPoints, dataEnv, xlim=xlim, ylim=ylim,
                   type=type, lty=lty, col=colors, labEnv=labEnv,
                  xloc=xloc,  ...)
    }
    ## Create an environment that contains the necessary X & Y points
    ## for use with identify()
    identEnv <- new.env()
    multiassign(c("X","Y"),list(xPoints,chromExprs),envir=identEnv)

    return(identEnv)
}

getACStrandVals <- function(chromExprs, strands, xPoints, dup,
                           geneNames, plotFormat, xloc="equispaced") {
    ## Determine which points are on the + and which on the -
    ## strand
    posPoints <- xPoints[strands %in% "+"]
    negPoints <- xPoints[strands %in% "-"]

    if (plotFormat == "cumulative") {
        posExprs <- chromExprs[which(strands=="+"),]
        negExprs <- chromExprs[which(strands=="-"),]
    }
    else {
        posExprs <- negExprs <- chromExprs
        posExprs[negPoints,] <- 0
        negExprs[posPoints,] <- 0
    }

    if (xloc == "physical") {
        pts <- which(xPoints %in% posPoints)
        nts <- which(xPoints %in% negPoints)
        posDup <- posPoints[pts %in% dup]
        posDup <- match(posDup,posPoints)
        negDup <- negPoints[nts %in% dup]
        negDup <- match(negDup,negPoints)
    }
    else {
        pts <- posPoints
        nts <- negPoints
        posDup <- dup[dup %in% pts]
        negDup <- dup[dup %in% nts]
    }

    posGen <- geneNames[pts]
    posStr <- strands[pts]
    negGen <- geneNames[nts]
    negStr <- strands[nts]

    strandList <- list(posExprs=posExprs, negExprs=negExprs,
                       posPoints=posPoints, negPoints=negPoints,
                       pts=pts, nts=nts, posDup=posDup, negDup=negDup,
                       posGen=posGen, posStr=posStr, negGen=negGen,
                       negStr=negStr)
    return(strandList)
}

dispACXaxis <- function(xPoints, dataEnv, xloc="equispaced",
                        plotFormat="local") {
    ## Retrieve values from dataEnv
    chromExprs <- dataEnv$chromExprs
    geneNames <- dataEnv$geneNames
    strands <- dataEnv$strands
    byStrand <- dataEnv$byStrand
    dup <- dataEnv$dup

    ## Make sure that xPoints isn't exceeding our visual maximum.
    ## If so, reduce the number of poitns to actually be displayed.
    dispXPoints <- cullACXPoints(xPoints)
    dispPointLocs <- match(dispXPoints,xPoints)

    if (any(dup))
        highlightACDups(dispXPoints, chromExprs, dup, xloc)

    if (plotFormat == "cumulative") {
        ## Need to filter out the first and last tick
        dispXPoints <- dispXPoints[2:(length(dispXPoints)-1)]
        dispPointLocs <- dispPointLocs[2:(length(dispPointLocs)-1)]
    }

    axis(1, at=dispXPoints, labels = geneNames[dispPointLocs], las=2,
         cex.axis=0.7,)
    if (byStrand == FALSE) {
        axis(3, at=dispXPoints, labels = strands[dispPointLocs],
             cex.axis=0.7, tick=FALSE, mgp=c(0,0,0))
    }
}

getACPlotLabs <- function(plotFormat, chrom, xloc, scale) {
    labEnv <- new.env()

    ylab <- switch(plotFormat,
                   "cumulative"="Cumulative expression levels",
                   "local"="Expression levels",
                   "image"="Samples"
                   )

    xlab <- "Representative Genes"
    main <- buildACMainLabel(ylab, chrom, xloc, plotFormat, scale)
    multiassign(c("xlab","ylab","main"),c(xlab,ylab,main),envir=labEnv)
    return(labEnv)
}

getACDataEnv <- function(chromExprs, geneNames, strands, byStrand,
                         dup) {
    dataEnv <- new.env()
    titles <- c("chromExprs","geneNames","strands","byStrand","dup")
    vals <- list(chromExprs, geneNames, strands, byStrand, dup)
    multiassign(titles, vals, envir=dataEnv)
    return(dataEnv)
}

highlightACDups <- function(xPoints, chromExprs, dup, xloc) {
    y <- min(chromExprs)-0.2

    for (i in seq(along=dup)) {
        ## For each dup, see if both that point and the point
        ## before it are still in the displayed set of points
        cur <- dup[i]
        prev <- dup[i] - 1
        if (xloc == "equispaced") {
            curPt <- match(cur, xPoints)
            prevPt <- match(prev, xPoints)
        }
        else {
            curPt <- cur
            prevPt <- prev
        }
        if ((!is.na(curPt))&&(!is.na(prevPt))) {
            segments(xPoints[curPt],y,xPoints[prevPt],y, col="cyan",lwd=2)
        }
    }
}

fixACPhysPoints <- function(xPoints, dup) {
    ## !!!!!
    ## !!! Currently doing this in a very inefficient manner.
    ## !!! needs to be smarter
    ## !!!!!!

    if (any(dup)) {
        dupDiff <- c(1,diff(dup),2)
        tmpDup <- NULL
        for (i in 1:(length(dup)+1)) {
            if (dupDiff[i] != 1) {
                ## At end of dup run
                dist <- xPoints[tmpDup[length(tmpDup)]+1] - xPoints[tmpDup[1]]
                spacing <- dist/(length(tmpDup)+1)
                for (j in 1:length(tmpDup)) {
                    pt <- dup[match(tmpDup[j],dup)]
                    xPoints[pt] <- xPoints[pt] + (j*spacing)
                }
                tmpDup <- NULL
            }
            tmpDup <- c(tmpDup,dup[i])
        }
    }
    return(xPoints)
}

buildACMainLabel <- function(ylab, chrom, xloc, plotFormat, scale) {
    if ((xloc == "physical")&&(plotFormat=="cumulative")) {
        main <- paste(ylab, "in chromosome", chrom,
                      "by relative position\n")
    }
    else {
        main <- paste(ylab, "by genes in chromosome", chrom, "\n")
    }

    main <- paste(main,"scaling method:",scale,"\n")

    return(main)
}

identifyLines <- function(identEnvir, ...) {
    ## Will call identify() on teh alongChrom() plot to detail which
    ## lines map tow which samples

    points <- mget(c("X","Y"), envir=identEnvir, ifnotfound=NA)

    xPoints <- points$X
    yPoints <- points$Y

    x <- identify(rep(xPoints,ncol(yPoints)), yPoints,
                  labels=rep(colnames(yPoints),
                  rep(nrow(yPoints),ncol(yPoints))), ...)

    return(x)
}

limitACXRange <- function(xlim, usedGenes) {

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
                    lowLim <- getACClosestPos(xlim[1],usedGenes)
                }

                hiLim <- match(xlim[2], usedGenes)
                if (is.na(hiLim)) {
                    hiLim <- getACClosestPos(xlim[2],usedGenes)
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

getACGeneSyms <- function(affys, chrObj) {
    syms <- mget(affys, env=geneSymbols(chrObj), ifnotfound=NA)
    syms[is.na(syms)] <- affys[is.na(syms)]
    return(as.character(syms))
}

getACClosestPos <- function(val, usedGenes) {
    ## Given a value, finds the closest value in usedGenes to the
    ## passed value and returns its location in the usedGenes vector

    dists <- abs(val-abs(as.numeric(usedGenes)))
    closest <- match(min(dists), dists)
    return(closest)
}

scaleACData <- function(chromData,
                        method=c("none","zscale","rangescale","rankscale",
                        "zrobustscale"))
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

cullACXPoints <- function(xPoints) {
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

emptyACPlot <- function(chrom) {
    plot.new()
    axis(1,labels=rep("NA",6))
    axis(2, labels=rep("NA",6))
    main <- paste("Plot empty, no genes from chromosome",chrom,
                  "in exprSet provided.\n")

    title(main = main)
}

getACExprs <- function(eSet, usedGenes,
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

    chromExprs <- scaleACData(chromExprs,scale)

    if (plotFormat == "cumulative") {
        chromExprs <- t(chromExprs)
        ## Fill the matrix with the cumulative sum of the expression
        chromExprs <- apply(chromExprs, 1, cumsum)
    }

    return(chromExprs)
}
