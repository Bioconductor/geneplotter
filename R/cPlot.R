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


cPlot <- function(xPoints, fileName) {
    # Passed the number of points on the X axis (for resolution) and
    # the data file name, and will plot the data from the chromosome file

    # Get the scaling factor
    scale <- xScale(xPoints);

    # Build the initial plot structure
    plot(c(1, xPoints), c(1,24), type="n", xlab="", ylab="Chromosomes",
    axes=FALSE,)
    labs <- c("Y","X",22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)
    axis(2, c(1:24), labs)

    for (i in 1:24) lines(c(1,xPoints-1),c(i,i),col="blue")

    # Read in the file
    load(fileName);

    for (i in 1:length(byChroms)) {
        plotData(i,byChroms[[i]],scale[i]);
    }
}


