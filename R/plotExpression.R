plotExpressionGraph <- function(graph, nodeLLmap, exprs, LOCUSIDenvir,
                                mapFun, log=FALSE, nodeAttrs=list(), ...) {
    require("Rgraphviz") || stop("Requires Rgraphviz to continue")

    envll <- unlist(contents(LOCUSIDenvir))
    graphLLs <- unlist(lapply(nodeLLmap, function(x){x[1]}))
    graphAffys <- names(envll)[envll %in% graphLLs]

    if (missing(mapFun))
        mapFun <- defMapFun

    cols <- getPlotExpressionColors(graphAffys, exprs, mapFun, log)

    ## Vector of colors w/ affy's as names - need SYMs
    colAffys <- names(cols)
    colLLs <- envll[colAffys]
    colSyms <- names(graphLLs[graphLLs %in% colLLs])
    names(cols) <- colSyms
    nodeAttrs$fillcolor <- cols

    plot(IMCAGraph, nodeAttrs=nodeAttrs, ...)
}


getPlotExpressionColors <- function(graphAffys, exprs, mapFun, log=FALSE) {

    if (missing(mapFun))
        mapFun <- defMapFun

    affyCols <- mapFun(exprs, log)

    affyCols[names(affyCols) %in% graphAffys]
}

defMapFun <- function(exprs, log=FALSE) {
    part1 <- 100
    part2 <- 500

    if (log) {
        part1 <- log2(part1)
        part2 <- log2(part2)
    }

    cols <- unlist(lapply(exprs, function(x) {
        if (x <= part1)
            "blue"
        else if (x <= part2)
            "green"
        else
            "red"
    }))

    cols
}

