# gPlot.R
# Will read in a tab deliminated file of genome information
# and proceed to generate a visual representation of the
# chromosomes.

plotData <- function(chromNum, locs, scale)
{
  nlocs <- length(locs)

  locs <- locs*scale

  ypos <- rep(24-chromNum+1, nlocs)
  ytop <- ifelse(locs>0, ypos+0.4, ypos-0.4)

  segments(abs(locs), ypos, abs(locs), ytop, col="grey")
}

cPlot <- function(xPoints, plotChroms) {
    # Passed an instance of a chromLocation class, and the number of
    # points to represent on the X axis, will utilize that data
    # to plot a set of genes on their proper chromosome locations.

    # Get the scaling factor
    scale <- cScale(xPoints, chromLengths(plotChroms));

    # Build the initial plot structure
    plot(c(1, xPoints), c(1,24), type="n", xlab="", ylab="Chromosomes",
    axes=FALSE,)
    labs <- chromNames(plotChroms)
    axis(2, c(1:nChrom(plotChroms)), labs)

    for (i in 1:nChrom(plotChroms)) lines(c(1,xPoints-1),c(i,i),col="blue")

    byChroms <- chromLocs(plotChroms)

    for (i in 1:nChrom(plotChroms)) {
        plotData(i,byChroms[[i]],scale[i]);
    }
}


