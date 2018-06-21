#' Discrepancy maps for February 1981-2011
#'
#' The object \code{discrep} is obtained from running the \code{createMapping} function
#' for the month of February from 1981-2011. The predictions are from the CM2.5
#' Forecast-oriented Low-Ocean Resolution (FLOR) model
#' produced by the National Oceanic and Atmospheric Administration’s Geophysical Fluid Dynamics
#' Laboratory converted to a Polar Stereographic grid at a 3.5-month lead time  (Vecchi et al. 2014; Msadek et al. 2014).
#'  Weights for converting to a polar stereograhic grid were obtained
#'  from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#' The observations are from the monthly sea ice concentration
#' obtained from the National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. The results
#' are distributed by the National Snow and Ice Data Center (Comiso 2000, updated 2015).
#' @docType data
#' @format Object obtained from the \code{createMapping} function (see details)
#'
#' @details The object \code{discrep} is obtained from running the \code{createMapping} function. It is a list of five objects where \code{month},
#' \code{startYear}, and \code{endYear} give the month, first year, and last year that were mapped. The variables \code{obsList} and \code{predList}
#' are lists of arrays with one 3-dimensional array for each region. The first dimension is for the year. The other two dimensions
#' are for the fixed points' y-coordinates, the mapped points' x-coordinates, the mapped points' y-coordinates, the length of the mapping vectors in the
#' x-direction, the length of the vectors in the y-direction, and the angles of the mapping vectors.
#'
#' @keywords datasets
#' @references
#' Bootstrap sea ice concentration:
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#'
#' CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:
#' Vecchi, Gabriel A., et al.
#' \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#' cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#' Msadek, R., et al.
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
#'  Geophysical Research Letters 41.14 (2014): 5208-5215.
#'
#' National Center for Atmospheric Research, 2017: Earth system grid at NCAR. \url{https://www.
#' earthsystemgrid.org/home.html}.
#'
#' Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate remapping and interpolation package."
#' Los Alamos National Laboratory, Los Alamos, NM (1997).
#' @examples
#' data(discrep)
#' names(discrep)
"discrep"


#' Observed sea ice February 1981-1982
#'
#' The object \code{observed} is an array obtained from the function \code{readMonthlyBS}
#' for startYear = 1981 and endYear = 1982. It gives the observed sea ice
#' concentrations arranged in an array of dimension of year x month x lon x lat. The observations are from the monthly sea ice concentration
#' obtained from the National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. The results
#' are distributed by the National Snow and Ice Data Center (Comiso 2000, updated 2015).
#' @docType data
#' @format array of dimension of 2 x 12 x 304 x 448 (year x month x longitude x latitude)y
#' @keywords datasets
#' @references
#' Bootstrap sea ice concentration:
#'
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#' @examples
#' data(obsFeb19811982)
#' dim(obsFeb19811982)
"obsFeb19811982"

#' Observed sea ice February 2012
#'
#' The object \code{observed} is an array obtained from the \code{readMonthlyBS} function
#' for startYear = 2012 and endYear = 2012. It gives the observed sea ice
#' concentrations arranged in an array of dimension of year x month x lon x lat
#' @docType data
#' @format array of dimension of 2 years x 12 months  x 304 longitudes x 448 latitudes
#' @keywords datasets
#' @references
#' Bootstrap sea ice concentration:
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#' @examples
#' data(obsFeb2012)
#' dim(obsFeb2012)
"obsFeb2012"

#' Ensemble mean sea ice concentration for February 1981-1982 (initialized in November)
#'
#' The object \code{emFeb19811982} is an array of the ensemble mean of the predicted sea ice concentration
#' from the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model
#' produced by the National Oceanic and Atmospheric Administration’s Geophysical Fluid Dynamics
#' Laboratory (Vecchi et al. 2014; Msadek et al. 2014).
#' The data have been converted to a Polar stereographic grid.
#' Weights for converting to a polar stereograhic grid were obtained
#' from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#' @docType data
#' @format array of dimension of 2 x 304 x 448 (corresponding to year x longitude x 448 latitude)
#' @keywords datasets
#' @references
#' Bootstrap sea ice concentration:
#'
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#'
#' CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:
#'
#' Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate remapping and interpolation package."
#' Los Alamos National Laboratory, Los Alamos, NM (1997).
#'
#' Msadek, R., et al.
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
#'  Geophysical Research Letters 41.14 (2014): 5208-5215.
#'
#' National Center for Atmospheric Research, 2017: Earth system grid at NCAR. \url{https://www.
#' earthsystemgrid.org/home.html}.
#'
#' Vecchi, Gabriel A., et al.
#' \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#' cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#' @examples
#' data(emFeb19811982)
#' dim(emFeb19811982)
"emFeb19811982"



#' Ensemble mean sea ice concentration for February 2012 (initialized in November)
#'
#' The object \code{emFeb2012} is an array of the ensemble mean of the predicted sea ice concentration
#' from the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model
#' produced by the National Oceanic and Atmospheric Administration’s Geophysical Fluid Dynamics
#' Laboratory (Vecchi et al. 2014; Msadek et al. 2014).
#' The data have been converted to the Polar stereographic grid.
#' Weights for converting to a polar stereograhic grid were obtained
#' from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#' @docType data
#' @format matrix of dimension 304 x 448 (longitude x latitude)
#' @keywords datasets
#' @references
#' Bootstrap sea ice concentration:
#'
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#'
#' CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:
#'
#' Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate remapping and interpolation package."
#' Los Alamos National Laboratory, Los Alamos, NM (1997).
#'
#' Msadek, R., et al.
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
#'  Geophysical Research Letters 41.14 (2014): 5208-5215.
#'
#' National Center for Atmospheric Research, 2017: Earth system grid at NCAR. \url{https://www.
#' earthsystemgrid.org/home.html}.
#'
#' Vecchi, Gabriel A., et al.
#' \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#' cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#' Vecchi, Gabriel A., et al.
#' \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#' cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#' @examples
#' data(emFeb2012)
#' dim(emFeb2012)
"emFeb2012"


#' Example of mapped points
#'
#' Example of a set of mapped points organized as an n x 2 matrix of coordinates.
#' This is used to demonstrate the \code{makePolygons} function.
#' @docType data
#' @format matrix of 1027 x 2
#' @keywords datasets
#' @examples
#' data(mappedPoints)
#' head(mappedPoints)
#' plot(mappedPoints, type = "l")
"mappedPoints"

#' Example of a line that contains self-intersections
#'
#' Example of a line that contains self-intersections. We will use it to
#' demonstrate the functions that address these intersections.
#' @docType data
#' @format n x 2 matrix of coordinates
#' @keywords datasets
#' @examples
#' data(interEx)
#' plot(interEx)
"interEx"

#' Coordinates of an observed line segment
#'
#' Example of the coordinates for an observed line segment. We will use it to
#' demonstrate the \code{intLine} function.
#' @docType data
#' @format n x 2 matrix of coordinates
#' @keywords datasets
#' @examples
#' data(obsLEx)
#' head(obsLEx)
"obsLEx"

#' Coordinates of a predicted line segment
#'
#' Example of the coordinates for a predicted line segment. We will use it to
#' demonstrate the \code{intLine} function.
#' @docType data
#' @format n x 2 matrix of coordinates
#' @keywords datasets
#' @examples
#' data(predLEx)
#' head(predLEx)
"predLEx"

#' Coordinates of a line segment with self-intersections
#'
#' Example of a line segment with self-intersections. We will use it
#' to demonstrate the \code{untwistSec} function.
#' @docType data
#' @format n x 2 matrix of coordiantes
#' @keywords datasets
#' @examples
#' data(currSecEx)
#' head(currSecEx)
"currSecEx"

#' Binary matrix indicating where there is land
#'
#' Binary matrix of dimension 304 x 448 with value for 1 for land grid boxes and 0 otherwise.
#' Data are on a north Polar Stereographic grid with the land mask simplified
#' to match model output from the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model
#' produced by the National Oceanic and Atmospheric Administration’s Geophysical Fluid Dynamics
#' Laboratory converted to a Polar Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014).
#' Weights for converting to a polar stereograhic grid were obtained
#' from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#' @docType data
#' @format 304 x 448 matix
#' @keywords datasets
#' @references
#' CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model:
#'
#' Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate remapping and interpolation package."
#' Los Alamos National Laboratory, Los Alamos, NM (1997).
#'
#' Msadek, R., et al.
#' \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions of Arctic sea ice extent."}
#'  Geophysical Research Letters 41.14 (2014): 5208-5215.
#'
#' National Center for Atmospheric Research, 2017: Earth system grid at NCAR. \url{https://www.
#' earthsystemgrid.org/home.html}.
#'
#' Vecchi, Gabriel A., et al.
#' \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical}
#' cyclone activity." Journal of Climate 27.21 (2014): 7994-8016.
#'
#' @examples
#' data(landMat)
#' image(landMat, xaxt = "n", yaxt = "n")
"landMat"


