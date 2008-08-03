groupedHeatmap = function(
  z, frow, fcol,
  fillcolours = c("#2166ac","#4393c3","#92c5de","#d1e5f0","#fefefe","#fddbc7","#f4a582","#d6604d","#b2182b"),
  bordercolour = "#e0e0e0",
  zlim = range(z, na.rm=TRUE)) {

  ## Define set of vertical and horizontal
  ## lines along which the plot is organised
  ## s: a character vector with the strings for the labels of the
  ##    *other* coordinate axis
  ## g: a factor with groups
  makecoords = function(s, g) {
    stopifnot(is.factor(g))
    x0 = if(is.null(s)) unit(0, "npc") else max(convertUnit(stringWidth(s), "mm"))
    wx = unit(1, "npc") - x0
    gapsize = 0.5
    dx = wx* ( 1 / (nlevels(g)*gapsize + length(g) -0.5) )
    return(list(
      pos = x0 + ((0L:(length(g)-1L))+(as.integer(g)-1L)*gapsize+0.5) * dx,
      delta = dx) )           
  }

  ## map data values into fillcolours
  colourMap = function(z, numColours = 201, na.colour="#ffffff"){
    colores = colorRampPalette(fillcolours)(numColours)
    i = as.integer(round( (z-zlim[1]) / diff(zlim) * numColours) )
    i[i<1L] = 1L
    i[i>numColours] = numColours
    list(fill = ifelse(is.na(z), na.colour, colores[i]),
         col  = ifelse(is.na(z), na.colour,
           if (is.null(bordercolour)||is.na(bordercolour)) colores[i] else bordercolour))
  }

  if(missing(frow)) {
    frow = factor(rep(1L, nrow(z)))
  } else {
    o = order(frow)
    z = z[o, ]
    frow = frow[o]
  }

  if(missing(fcol)) {
    fcol = factor(rep(1L, ncol(z)))
  } else {
    o = order(fcol)
    z = z[, o]
    fcol = fcol[o]
  }

  textx = if(is.null(colnames(z))) NULL else paste(colnames(z), "", sep=" ")
  texty = if(is.null(rownames(z))) NULL else paste(rownames(z), "", sep=" ")
  
  cx = makecoords(s=texty, g=fcol)  
  cy = makecoords(s=textx, g=frow)
  x  = cx$pos
  y  = cy$pos

  grid.rect(x = x[rep(seq(along=x), each  =length(y))], width  = cx$delta, 
            y = y[rep(seq(along=y), times =length(x))], height = cy$delta,
            just = c(0.5,0.5), gp = do.call(gpar, colourMap(z)))
  
  if(!is.null(textx)) grid.text(textx, x=x, y=y[1]-0.5*cy$delta, just=c("right", "center"), rot=90)
  if(!is.null(texty)) grid.text(texty, x=x[1]-0.5*cx$delta, y=y, just=c("right", "center"))
}
