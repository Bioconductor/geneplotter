# cScale.R
# Used to calculate scaling on the geneplots
# Uses the vector of chromosome lengths and returns a vector
# of scales.

cScale <- function(points, cLengths, method=c("max","relative"),
                   chrom) {
# Passed points - the number of points to scale the chromosomes too
# and cLengths - a vector of chromosome lengths.

    method <- match.arg(method)

    if (method == "max") {
            cScales <- points / cLengths[chrom];
    }
    else {
        cScales <- points / max(cLengths)
    }

    return(cScales);
}
