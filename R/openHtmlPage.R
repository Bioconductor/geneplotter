openHtmlPage = function(name, title="") {
  con = file(paste(name, ".html", sep=""), open="wt")
  writeLines(paste("<html><head><title>", title, "</title></head><body style=\"font-family: helvetica,arial,sans-serif;\">", sep=""), con)
  return(con)
}

closeHtmlPage = function(con) {
  writeLines("</body></html>", con)
  close(con)
}

