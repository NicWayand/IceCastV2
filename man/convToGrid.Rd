% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importExport.R
\name{convToGrid}
\alias{convToGrid}
\title{Convert \code{SpatialPolygons} object to a grid}
\usage{
convToGrid(x, myLandMat = landMat)
}
\arguments{
\item{x}{\code{SpatialPolygon} object}

\item{myLandMat}{binary matrix specifying land locations}
}
\description{
Convert \code{SpatialPolygons} object to binary grid. Using the \code{rasterize} function,
grid boxes whose centers are part of the \code{SpatialPolygons} are given value 1 and all other
grid boxes are given value 0. Land values are set to NA.
}
\examples{
grid <- convToGrid(bgWater)
image(grid)
}
