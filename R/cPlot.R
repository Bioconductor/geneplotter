## cPlot.R
## Will read in a tab deliminated file of genome information
## and proceed to generate a visual representation of the
## chromosomes.

.plotData <- function(chromNum, locs, xPoints, chromLens, fg,
                      scale = c("max","relative")[1],glen=0.4)
{
    ## Get the scaling factor
    scale <- cScale(xPoints, chromLens, scale)

    nlocs <- length(locs)

    ## APply the scaling factor to the x positions
    locs <- locs*scale[chromNum]

    ## Determine the direction of the Y plot (+ or -)
    ypos <- rep(chromNum, nlocs)
    ytop <- ifelse(locs>0, ypos+glen, ypos-glen)

    if (scale == "max") {
        lines(c(1,xPoints-1),c(chromNum,chromNum),col="blue")
    }
    else {
        lines(c(1,max(abs(locs))),c(chromNum,chromNum),col="blue")
    }

    ## Plot
    segments(abs(locs), ypos, abs(locs), ytop, col=fg)
}

cColor <- function(genes, color, plotChroms) {
    ## Passed a vector of gene names, a color and an instance of a
    ## chromLocation class.  Will recolor the specific genes in the
    ## cPlot created plot to match the specified color

    xPoints <- 1000

    ## Get the chromLocs listing from the chromLocation class
    locList <- chromLocs(plotChroms)

    ## Look up the locations of these genes in each chromosome,
    ## plotting any results.
    for (i in 1:nChrom(plotChroms)) {
        locs <- locList[[i]][genes]
        locs <- as.numeric(locs[!is.na(locs)])
        if (length(locs) > 0) {
            .plotData(i, locs, xPoints, chromLengths(plotChroms), color)
        }
    }
}

cPlot <- function(plotChroms, useChroms=chromNames(plotChroms),
                  scale=c("max","relative")[1], fg="white",
                  bg="lightgrey", glen=0.4) {
    ## Passed an instance of a chromLocation class, and the number of
    ## points to represent on the X axis, will utilize that data
    ## to plot a set of genes on their proper chromosome locations.

    xPoints <- 1000
    glen <- glen

    chromNames <- chromNames(plotChroms)
    labs <- rev(chromNames[chromNames %in% useChroms])
    lens <- chromLengths(plotChroms)
    lens <- rev(lens[chromNames %in% labs])

    ## Build the initial plot structure
    op <- par(bg=bg)
    plot(c(1, xPoints), c(1-glen,length(labs)+glen), type="n", xlab="",
         ylab="Chromosomes", axes=FALSE)
    par(op)

    axis(2, c(1:length(labs)), labs)

    byChroms <- chromLocs(plotChroms)

    for (i in 1:length(labs)) {
        .plotData(i,byChroms[[labs[i]]], xPoints, lens, fg, scale,glen);
    }
}



