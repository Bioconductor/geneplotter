#copyright R. Gentleman, 2001, all rights reserved
#functions/ methods to plot microarray data

#Cheng Li's Red/Blue color scheme

dChip.colors <- function(n) GetColor(seq(-3,3,6/n))

greenred.colors <- function(n) GetColor(seq(-3,3,6/n), GreenRed=TRUE)

GetColor <- function(value, GreenRed=FALSE, DisplayRange=3) {
    RGB <- function(x,y,z) rgb(x/255,y/255,z/255)
    missing <- is.na(value)
    good <- value[!missing]
    ans <- value
    if ( GreenRed )
        ans[missing] <- RGB(0, 0, 0)
    else
        ans[missing] <- RGB(255, 255, 255)
    tone <- abs(good) / DisplayRange * 255 + .5;
    tone[tone > 255] <- 255

    #classical: red and green on black background
    if (GreenRed)
        tone <- ifelse(good > 0, RGB(tone, 0, 0), RGB(0, tone, 0))
    else
        tone <- ifelse(good > 0, RGB(255, 255 - tone, 255 - tone),
                       RGB(255 - tone, 255 - tone, 255) )
    ans[!missing]<-tone
    ans
}
