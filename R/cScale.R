# cScale.R
# Used to calculate scaling on the geneplots
# Uses the vector of chromosome lengths and returns a vector
# of scales.

cScale <- function(cLength) {
# Passed cLength - the number of points to scale the chromosomes too

    data(cLengths);

    nLen <- length(cLengths);
    cScales <- vector("numeric", length=nLen);

    for (i in 1:nLen) {
        cScales[i] <- cLength / cLengths[i];
    }

    return(cScales);
}
