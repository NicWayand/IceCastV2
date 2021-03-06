% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/storedResults.R
\docType{data}
\name{landMat}
\alias{landMat}
\title{Binary matrix indicating where there is land}
\format{304 x 448 matix}
\usage{
landMat
}
\description{
Binary matrix of dimension 304 x 448 with value for 1 for land grid boxes and 0 otherwise.
Data are on a north Polar Stereographic grid with the land mask simplified
to match model output from the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model
produced by the National Oceanic and Atmospheric Administration’s Geophysical Fluid Dynamics
Laboratory converted to a Polar Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014).
Weights for converting to a polar stereograhic grid were obtained
from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
}
\examples{
data(landMat)
image(landMat, xaxt = "n", yaxt = "n")
}
\references{
CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:

Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate remapping and interpolation package."
Los Alamos National Laboratory, Los Alamos, NM (1997).

Msadek, R., et al.
\href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
 Geophysical Research Letters 41.14 (2014): 5208-5215.

National Center for Atmospheric Research, 2017: Earth system grid at NCAR. \url{https://www.
earthsystemgrid.org/home.html}.

Vecchi, Gabriel A., et al.
\href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
}
\keyword{datasets}
