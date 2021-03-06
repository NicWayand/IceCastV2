% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selfIntersections.R
\name{checkIntersect}
\alias{checkIntersect}
\title{Check if line segments intersect}
\usage{
checkIntersect(a, b, c, d, seq = FALSE)
}
\arguments{
\item{a}{first coordinate of first line segment}

\item{b}{second coordinate of first line segment}

\item{c}{first coordinate of second line segment}

\item{d}{second coordinate of second line segment}

\item{seq}{indicator for whether the two line segments are intersecting}
}
\value{
boolean indicating if there is an intersection
}
\description{
Find if two line segments intersect
}
\examples{
 checkIntersect(c(0, 0), c(1, 1), c(2, 2), c(3, 3))
 checkIntersect(c(0, 0), c(1, 1), c(0.5, 0.5), c(2, 2))
}
