plotChr <- function(chrN, senseObj,
	cols=rep("blue", length(senseObj[[1]])),log=FALSE,
	xloc = c("equispaced", "physical"), geneSymbols=FALSE,
	ngenes=20, lines.at=NULL, lines.col="red") {
    # lines of +/- stands of a chromosome for a given sample:
    linesStrand <- function(smooths, col="blue", log=FALSE,
        smX=NULL) {
    	if (is.null(smX)) {
        	sm.px <- smooths$pos$x
        	sm.nx <- smooths$neg$x
    	} else {
        sm.px <- smX[1:length(smooths$pos$x)]
        sm.nx <- smX[-(1:length(smooths$pos$x))]
    	}
    	if( log ) {
        	lines(sm.px, log(smooths$pos$y), col=col)
        	lines(sm.nx, -log(-smooths$neg$y), col=col)
    	}
    	else {
        	lines(sm.px, smooths$pos$y, col=col)
        	lines(sm.nx, smooths$neg$y, col=col)
    	}
    }
    xloc <- match.arg(xloc)
    ans2 <- senseObj$ans2
    libCHRLENGTHS <- get(paste(senseObj$lib,"CHRLENGTHS",sep=""))
    chrRange <- function(chrN) {
      mn <- min(sapply(ans2, function(x) min(x[[chrN]]$neg$y)))
      mx <- max(sapply(ans2, function(x) max(x[[chrN]]$pos$y)))
      c(mn, mx)
    }
    xlims <- c(0, libCHRLENGTHS[chrN])
    ylims <- chrRange(chrN)
    if( log ) {
        ylims[1] <- -log(-ylims[1])
        ylims[2] <- log(ylims[2])
    }
    at.px <- ans2[[1]][[chrN]]$pos$x
    at.nx <- ans2[[1]][[chrN]]$neg$x
    X <- c(at.px, at.nx)
    uX <- !duplicated(X)
    every <- sum(uX) %/% ngenes
    if (every < 1) every <- 1
    ind.seq <- seq(1,sum(uX),by=every)
    repGx <- sort(X[uX])[ind.seq]
    probes <- names(X[uX])[order(X[uX])][ind.seq]
    if (xloc=="equispaced") {
		repGx <- 1:length(repGx)
		names(repGx) <- probes
		xlims <- c(1, length(repGx))
    }
    pn <- rep(c("+","-"), c(length(at.px),
	length(at.nx)))[uX][order(X[uX])][ind.seq]
    # probe density:
    if (xloc=="equispaced") {
	at.px <- approx(sort(X[uX])[ind.seq],repGx, xout=at.px)$y
        at.nx <- approx(sort(X[uX])[ind.seq],repGx, xout=at.nx)$y
	smX <- c(at.px, at.nx)
    } else smX <- NULL
    opar <- par(mar = c(6, 5, 4, 1), mgp = c(4, 1, 0))
    on.exit(par(opar), add = TRUE)
    if (log == TRUE)
        yLab <- "Smoothed Expression (log)"
    else
        yLab <- "Smoothed Expression"
    plot(1,1, type="n", xlim=xlims, ylim=ylims, cex.lab=0.9,
	 xlab="Representative Genes",
         ylab=yLab, main=paste("Chromosome",chrN),
        xaxt="n", yaxt="n")
    yticks <- pretty(c(0,ylims), 5)
    axis(2, at=yticks, labels=abs(yticks))
    abline(h=0, col="gray")
    axis(1, at=at.nx, pos=0, tck=-0.01,col="gray", labels=FALSE)
    axis(1, at=at.px, pos=0, tck=0.01,col="gray", las=3, labels=FALSE)
    # label representative genes:
    labs <- probes
    if(geneSymbols) labs <- unlist(multiget(labs,
	env=get(paste(senseObj$lib,"SYMBOL",sep=""))))
    axis(1, at=repGx, labels=labs,las=3, cex.axis=.7)
    axis(3, at=repGx, labels=pn,tick=FALSE, line=-1,cex.axis=.7)
    for(i in 1:length(ans2))
        linesStrand(ans2[[i]][[chrN]], cols[i], log, smX=smX)
    if (!is.null(lines.at)) {
	lineXs <- unlist(multiget(lines.at,
        	env=get(paste(senseObj$lib,"CHRLOC",sep=""))))
	lineXs <- abs(lineXs)
	if(any(is.na(lineXs))) stop("wrong probe names")
        if (xloc=="equispaced") lineXs <-
		approx(sort(X[uX])[ind.seq],repGx, xout=lineXs)$y
    	abline(v=lineXs, col=lines.col)
    }
}

