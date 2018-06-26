#' Convert \code{SpatialPolygons} object to binary grid. Using the \code{rasterize} function,
#' grid boxes whose centers are part of the \code{SpatialPolygons} are given value 1 and all other
#' grid boxes are given value 0. Land values are set to NA.
#' @title Convert \code{SpatialPolygons} object to a grid
#' @param x \code{SpatialPolygon} object
#' @param myLandMat binary matrix specifying land locations
#' @export
#' @importFrom sp SpatialPolygons Polygons
#' @importFrom raster rasterize as.matrix
#' @examples
#' grid <- convToGrid(bgWater)
#' image(grid)
convToGrid <- function(x, myLandMat = landMat) {
  x <- aggregate(x)
  nPoly <- length(x@polygons[[1]]@Polygons)

  poly <- NULL
  for (i in 1:nPoly) {
    temp <- SpatialPolygons(list(Polygons(list(Polygon(x@polygons[[1]]@Polygons[[i]]@coords)), ID =  sprintf("temp%i", i))))
    temp <- untwist(temp, polyName = sprintf("temp%i", i))
    if (!is.null(poly) & !is.null(temp)) {
      poly <- spRbind(poly, temp)
    } else if (is.null(poly) & !is.null(temp)){
      poly <- temp
    }
    poly <- aggregate(poly) #makes into one polygon, so binary
  }
  rast <- raster(nrows = 448, ncols = 304, xmn = -3850, xmx = 3750,
                 ymn = -5350, ymx = 5850) #raster has columns and rows flipped from typical orientation
  rast <- rasterize(poly, rast, fun = max, background = 0)
  rast <- as.matrix(rast)
  rast <- t(rast)[,448:1]#fix weird orientation of raster
  rast[which(myLandMat == 1, arr.ind = T)] <- NA
  return(rast)
}

#' Reads in netCDF files of observations and predictions, performs bias correction,
#' and exports a new netCDF file with bias-corrected predictions
#' @title Simple evaluation of contour-shifting
#' @param obsNCDF filepath for observed data array (see details for info about array structure)
#' @param predNCDF filepath for predicted data array (see details for info about array structure)
#' @param predYears vectors of years for which to make prediction
#' @param startYear first year to use when learning model
#' @param month month of prediction
#' @param outputFile filepath for where bias-corrected netCDF file should be stored
#' @param datTypeObs string of either "bootstrap" or "simple" indicating the file type of the observation (see details for info about array structure)
#' @param level concentration level for which to build contour
#' @importFrom ncdf4 nc_open ncvar_get nc_close nc_create ncdim_def ncvar_def ncvar_put
#' @details The predicted data array, \code{predNCDF}, should be a netCDF file with a single array of dimension: years x longitude (304) x latitude (448).
#' The variable should be named \code{iceInd}. The values in the array should indicate whether each grid box is
#' categorized to contain ice (1: ice-covered, 0: no ice, NA: land).
#' The observed data array, \code{obsNCDF}, should be a netCDF file with a single array of dimension: years x longitude (304) x latitude (448).
#'  The observed data array, \code{obsNCDF}, can be formatted  the same as \code{predNCDF}
#' if \code{datTypeObs = "simple"}. Alternatively, if \code{datTypeObs = "bootstrap"} the array values can be ice concentration values
#' obtained from the National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. Data should be retained in the same format as given by bootstrap
#' (including indicators for missing data, land etc.). The variable should be named "conc".
#'
#' @references
#' Bootstrap sea ice concentration:
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#'
#'
#' @return netCDF file of dimension years by longitude (304) by latitude (448) with indicators for where ice is predicted after bias correction.
#'  (1: ice-covered, 0: not ice, NA: land). Grid boxes have been categorized as ice if their centers are ice covered (within R the bias-corrected contours are
#'  not restricted to align to a grid).
#' @examples
#' \dontrun{
#' quickRun(obsNCDF = "/obs.nc", predNCDF = "/pred.nc", predYears = c(2001:2013),
#'          startYear = 1980, month = 2, outputFile = "/outputFile.nc", level = 15,
#'          datTypeObs = "simple")
#' }
#' @export
quickRun <- function(obsNCDF, predNCDF, predYears, startYear, endYearOffset=1, month, outputFile,
                  level, datTypeObs = "bootstrap", plotting = FALSE) {

  #extract input dimensions
  obs <- nc_open(obsNCDF)
  if (datTypeObs == "bootstrap") {
    obsMat <- ncvar_get(obs, "conc")
  } else {
    obsMat <- ncvar_get(obs, "iceInd")
  }

  obsStartYear <- obs$dim$year$vals[1]
  pred <- nc_open(predNCDF)
  predMat <- ncvar_get(pred, "iceInd")
  predStartYear <- pred$dim$year$vals[1]

  #prepare for output
  yearDim <- ncdim_def("year", "years", predYears) # (name, units, vals...)
  monDim <- ncdim_def("month", "months", month)
  lonDim <- ncdim_def("lon", "longitude", 1:304)
  latDim <- ncdim_def("lat", "latitude", 1:448)
  nPredYear <- length(predYears)
  nmonths <- length(month)
  output <- array(dim = c(nPredYear, 304, 448, nmonths))

  #run mappings for all years
  print("Starting mapping...")
  discrep <- createMapping(startYear = startYear, endYear = max(predYears) - endYearOffset,
                           obsStartYear = obsStartYear, predStartYear = predStartYear,
                           observed = obsMat[,month,,], predicted = predMat[,month,,],
                           regions = regionInfo, month = month, level = level,
                           datTypeObs = datTypeObs, datTypePred = "simple", plotting = plotting)
  print("Done mapping.")
  print("Have mappings for years...")
  print(discrep$startYear)
  print(discrep$endYear)
  
  print("Starting bias correction...")
  
  #Bias correct predictions
  bgWater <- convToGrid(bgWater)
  for (k in 1:length(predYears)) {
    start.time <- Sys.time()
    adj <- contourShift(maps = discrep, predicted = predMat[length(predStartYear:predYears[k]), month,,],
                        bcYear = predYears[k], predStartYear = predStartYear, endYearOffset = endYearOffset, 
                        regions = regionInfo,
                        level = level, datTypePred = "simple")
    adj <- convToGrid(adj)
    adj[bgWater == 1] <- 2
    output[k, , , ] <- adj
    print(sprintf("Correction of year %i month %i completed", predYears[k], month))
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
  }
  print("Finished bias correction.")
  #define variables
  iceIndDef <- ncvar_def("iceInd", "indicator", list(yearDim, lonDim, latDim, monDim),
                       longname = "Indicator of if pixel is ice covered (0: not ice, 1: ice, NA: land, 2: Outside region")

  #create netCDF file and arrays
  ncOut <- nc_create(outputFile, iceIndDef)

  #put variables
  ncvar_put(ncOut, iceIndDef, output)

  #close file (writing file to disk)
  nc_close(ncOut)
}


