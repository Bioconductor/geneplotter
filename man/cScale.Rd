\name{cScale}
\alias{cScale}
\title{A function for mapping chromosome length to a number of points. }
\description {
  Given a number of points (generally representing the number of points
  on a plot's axis), will take a vector of chromosome lengths and create
  a vector of scale factors.  There is a one to one mapping between the
  scaling factors and the chromosomes.
}
\usage{
cScale(cLength)
}
\arguments{
  \item{cLength}{The number of points to scale the chromosome length to.}
}
\details{
  The scale factor is simply derived by dividing the cLength argument by
  the length of a particular chromosome.  These factors are then
  returned to the calling function.
}
\author{Jeff Gentry}
\examples {
     cScale(1000)
}
\keyword{utilities}

