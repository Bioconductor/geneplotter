# cScale.R
# Used to calculate scaling on the geneplots
# Uses the vector of chromosome lengths and returns a vector
# of scales.

cScale <- function(points, cLengths) {
# Passed points - the number of points to scale the chromosomes too
# and cLengths - a vector of chromosome lengths.

    nLen <- length(cLengths);
    cScales <- vector("numeric", length=nLen);

    for (i in 1:nLen) {
        cScales[i] <- points / cLengths[i];
    }

    return(cScales);
}
