## cPlot.R
.plotData <- function(chromNum, locs, xPoints, chromLens, fg,
                      scale = c("relative","max"),glen=0.4)
{
    ## Get the scaling factor
    scale <- match.arg(scale)

    scaledX <- cScale(xPoints, chromLens, scale, chromNum)

    nlocs <- length(locs)

    ## APply the scaling factor to the x positions
    cNum <- match(chromNum, names(chromLens))
    locs <- locs*scaledX
    if (length(locs) == 0) {
        if (scale == "relative")
            return()
    }
    else {
        ## Determine the direction of the Y plot (+ or -)
        ypos <- rep(cNum, nlocs)
        ytop <- ifelse(locs>0, ypos+glen, ypos-glen)

        ## Plot
        segments(abs(locs), ypos, abs(locs), ytop, col=fg)

        ## Drawn last to ensure that that the lines are actually displayed
    }
    if (scale == "max") {
        lines(c(1,xPoints-1),c(cNum,cNum),col="blue")
    }
    else {
        lines(c(1,max(abs(locs[!is.na(locs)]))),c(cNum,cNum),col="blue")
    }
}

cColor <- function(probes, color, plotChroms,
                   scale=c("relative","max"), glen=0.4) {
    ## Passed a vector of probe names, a color and an instance of a
    ## chromLocation class.  Will recolor the specific probes in the
    ## cPlot created plot to match the specified color.  Scale should
    ## be the same as the scale from cPlot
    scale <- match.arg(scale)
    xPoints <- 1000

    gc <- unlist(multiget(probes,env=probesToChrom(plotChroms)))
    gchr <- split(names(gc),gc)

    gchr[["NA"]] <- NULL

    ## Look up the locations of these probes in each chromosome,
    ## plotting any results.
    locList <- chromLocs(plotChroms)
    lens <- chromLengths(plotChroms)
    names(lens) <- chromNames(plotChroms)

    for (cName in names(gchr)) {
        locs <- locList[[cName]][gchr[[cName]]]
        locs <- as.numeric(locs[!is.na(locs)])
        if (length(locs) > 0) {
            .plotData(cName, locs, xPoints, lens,
                      color, scale, glen)
        }
    }
}


cPlot <- function(plotChroms, useChroms=chromNames(plotChroms),
                  scale=c("relative", "max"), fg="white",
                  bg="lightgrey", glen=0.4) {
    ## Passed an instance of a chromLocation class, and the number of
    ## points to represent on the X axis, will utilize that data
    ## to plot a set of genes on their proper chromosome locations.
    scale <- match.arg(scale)

    xPoints <- 1000

    chromNames <- chromNames(plotChroms)
    labs <- chromNames[chromNames %in% useChroms]

    lens <- chromLengths(plotChroms)
    whichLabs <- chromNames %in% labs
    lens <- lens[whichLabs]
    names(lens) <- chromNames[whichLabs]

    ## Build the initial plot structure
    op <- par(bg=bg)
    plot(c(1, xPoints), c(1-glen,length(labs)+glen), type="n", xlab="",
         ylab="Chromosomes", axes=FALSE, las=2, main=organism(plotChroms))
    par(op)

    axis(2, c(1:length(labs)), labs, las=2)

    chromLocs <- chromLocs(plotChroms)
    byChroms <- chromLocs[labs]

    for (cName in labs) {
        .plotData(cName,byChroms[[cName]], xPoints,
                  lens, fg, scale,glen);
    }
}



