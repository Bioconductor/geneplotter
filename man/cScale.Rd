\name{cScale}
\alias{cScale}
\title{A function for mapping chromosome length to a number of points. }
\description {
  Given a number of points (generally representing the number of points
  on a plot's axis), and a vector of chromosome lengths - will generate
  a vector of the same length as the one passed in containing scaling
  factors for each chromosome.
}
\usage{
cScale(points, cLengths)
}
\arguments{
  \item{points}{The number of points to scale the chromosome length to.}
  \item{cLengths}{A vector of chromosome lengths. }
}
\details{
  The scale factor is simply derived by dividing the cLength argument by
  the length of a particular chromosome.  These factors are then
  returned to the calling function.
}
\author{Jeff Gentry}
\examples {
     data(cLengths)
     scales <- cScale(1000, cLengths)
}
\keyword{utilities}

