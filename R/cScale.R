# cScale.R
# Used to calculate scaling on the geneplots
# Uses the vector of chromosome lengths and returns a vector
# of scales.

cScale <- function(points, cLengths, method=c("max","relative")[1]) {
# Passed points - the number of points to scale the chromosomes too
# and cLengths - a vector of chromosome lengths.

    nLen <- length(cLengths);
    cScales <- vector("numeric", length=nLen);

    if (method == "max") {
        for (i in 1:nLen) {
            cScales[i] <- points / cLengths[i];
        }
    }
    else {
        scale = points / max(cLengths)
        cScales <- rep(scale,nLen)
    }
    return(cScales);
}
