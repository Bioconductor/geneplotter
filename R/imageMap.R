if (is.null(getGeneric("imageMap")))
    setGeneric("imageMap", function(object, ...)
               standardGeneric("imageMap"))

setMethod("imageMap",
  signature=c("matrix", "connection", "list", "character"),
  definition=function(object, con, tags, imgname) {
    
  checkTags = function(x, tagname) {
    if(is.null(names(x))) {
      if(length(x)==length(nn)) {
        names(x)=nn
      } else {
        stop(paste("'tags$", tagname, "' must have names if it is not the ",
                   "same length as the number of nodes in 'object'.", sep=""))
      }
    } else {
      if(!all(names(x) %in% nn))
        stop(paste("'names(tags$", tagname, ")' must match the names of ",
                   "the nodes in 'object'", sep=""))
    }
    return(x)
  }
  for(i in seq(along=tags))
    tags[[i]] = checkTags(tags[[i]], names(tags)[i])

  if(ncol(object)!=4)
    stop("'object' must be a matrix with 4 columns.")

  mapname <- paste("map", gsub(" |/|#", "_", imgname), sep="_")
  base::writeLines(paste("<IMG SRC=\"", imgname, "\" USEMAP=\#", mapname, " BORDER=0>", 
                   "<MAP NAME=\"", mapname, "\">", sep=""), con)
  for(i in 1:nrow(object)) {
    out = paste("<AREA SHAPE=\"rect\" COORDS=\"", coords[i,1], ",", coord[i,2], ",",
                coords[i,3], ",", coords[i,4], "\"", sep="")
    for(t in seq(along=tags))
      out = paste(out, " ", names(tags)[i], "=\"", tags[[i]][t], "\"", sep="")
    out = paste(out, "\">", sep="")
    base::writeLines(out, con)
  }

} ## end of defintion
) ## end of setMethod
