## cPlot.R
## Will read in a tab deliminated file of genome information
## and proceed to generate a visual representation of the
#<<
                                        # chromosomes.

.plotData <- function(chromNum, locs, xPoints, chromLens, color,
                      scale = c("max","relative")[1])
{
  ## Get the scaling factor
    scale <- cScale(xPoints, chromLens, scale)

    nlocs <- length(locs)

  ## APply the scaling factor to the x positions
  locs <- locs*scale[chromNum]

  ## Determine the direction of the Y plot (+ or -)
  ypos <- rep(24-chromNum+1, nlocs)
  ytop <- ifelse(locs>0, ypos+0.4, ypos-0.4)

  ## Plot
  segments(abs(locs), ypos, abs(locs), ytop, col=color)
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

cPlot <- function(plotChroms, scale=c("max","relative")[1], ccol="lightgrey") {
    ## Passed an instance of a chromLocation class, and the number of
    ## points to represent on the X axis, will utilize that data
    ## to plot a set of genes on their proper chromosome locations.

    xPoints <- 1000

    ## Build the initial plot structure
    plot(c(1, xPoints), c(1,24), type="n", xlab="", ylab="Chromosomes",
    axes=FALSE,)
    labs <- rev(chromNames(plotChroms))
    axis(2, c(1:nChrom(plotChroms)), labs)

    for (i in 1:nChrom(plotChroms)) lines(c(1,xPoints-1),c(i,i),col="blue")

    byChroms <- chromLocs(plotChroms)

    for (i in 1:nChrom(plotChroms)) {
        .plotData(i,byChroms[[i]], xPoints, chromLengths(plotChroms),
        ccol, scale);
    }
}



