\name{openHtmlPage}
\alias{openHtmlPage}
\alias{closeHtmlPage}
\title{Open and close an HTML file for writing.}
\description{Open and close an HTML file for writing..}
\usage{
  openHtmlPage(name, title="")
  closeHtmlPage(con)
}
\arguments{
  \item{name}{Character. File name (\emph{without} the extension '.html').}
  \item{title}{Character. Value of the \code{title} tag in the HTML header.}
  \item{con}{Connection.}
}

\details{See example.}

\value{For \code{openHtmlPage}, a \code{\link[base::connection]{connection}}. }

\seealso{}

\author{Wolfgang Huber \url{http://www.dkfz.de/abt0840/whuber}}
\keyword{IO}
\examples{
  fn <- tempfile()
  con <- openHtmlPage(fn, "My page")
  writeLines("Hello world", con)
  closeHtmlPage(con)
  readLines(paste(fn, ".html", sep=""))
}


