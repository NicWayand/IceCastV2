###Misc. functions###

#Load data files that are needed globally (avoids CRAN error about global binding)
globalVariables(c("allRegions", "bgWater", "land", "landMat", "regionInfo"))

#' Finds the euclidean distance between two points (ignoring projection)
#' @title Find euclidean distance
#' @param p1 first point,  x and y coordinate pair
#' @param p2 second point,  x and y coordinate pair
#' @examples getDist(c(1, 2), c(3, 4))
#' @export
#' @return distance value
getDist <- function(p1, p2) {
  sqrt(((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2))
}

#' Remove holes from a \code{SpatialPolygons} object.  Note that this function differs from the function \code{findHoles} in that
#' it only removes holes contained within the polygon itself, not gaps between the polygon and region boundaries
#' @title Remove holes in a polygon
#' @param myPoly \code{SpatialPolygons} object
#' @param polyName character string to name polygon (defaults to "notSpecified")
#' @return \code{SpatialPolygons} object with holes removed
#' @export
#' @importFrom sp Polygon Polygons SpatialPolygons
#' @importFrom raster aggregate
#' @importFrom maptools spRbind
#' @importFrom methods is
#' @examples
#' withHoles <- bgWater[2]
#' plot(withHoles, col = "blue", main = "Polygon with Holes")
#' noHoles <- rmHoles(withHoles)
#' plot(noHoles, col = "blue", main = "Holes removed")
rmHoles <- function(myPoly, polyName = "notSpecified") {
  ##check myPoly input
  if (!(is(myPoly)[1] == "SpatialPolygons")) {
    stop("myPoly is not a polygon object")
  }

  ##remove holes
  myPoly <- suppressWarnings(aggregate(myPoly)) #give whole polygon one ID
  holes <- sapply(myPoly@polygons[[1]]@Polygons, function(x){x@hole}) #identify holes
  for (i in which(!holes)) {
    #current part of the polygon
    coords <- myPoly@polygons[[1]]@Polygons[i][[1]]@coords
    temp <- SpatialPolygons(list(Polygons(list(Polygon(coords, hole = FALSE)), sprintf("temp%i", i))))
    if (i == 1) { #make final version of the polygon if first iteration
      newPoly <-  temp
    } else { #add to existing final version of the polygon
      newPoly <- spRbind(newPoly, temp)
    }
  }
  newPoly <- aggregate(newPoly)
  newPoly@polygons[[1]]@ID <- polyName

  return(newPoly)
}

#' Keep only \code{SpatialPolygons} from a spatial object.
#' @title Keep only spatial polygons
#' @param myPoly \code{SpatialCollections}, \code{SpatialPolygons}, \code{SpatialPoints}, or \code{SpatialLines} object
#' @return \code{SpatialPolygons} object
#' @importFrom methods is
#' @export
#' @examples
#' par(mfrow = c(1, 2))
#' plot(spatialCollEx, col = "blue", main = "Spatial Collections Object")
#' polyOnly <- keepPoly(spatialCollEx)
#' plot(polyOnly, col = "blue", main = "Spatial Polygon Only")
keepPoly <- function(myPoly) {
  if (is(myPoly)[1] == "SpatialCollections") {
    myPoly <- myPoly@polyobj
  } else if (is(myPoly)[1] == "SpatialPoints") {
    return(NULL)
  } else if (is(myPoly)[1] == "SpatialLines") {
    return(NULL)
  } else if (is.null(myPoly)) {
    return(NULL)
  } else if (is(myPoly)[1] == "SpatialPolygons") {
    return(myPoly)
  }else {
    stop("Incorrect object type")
  }
  return(myPoly)
}

#' Keep only \code{SpatialLines} from a spatial object.
#' @title Keep only spatial lines
#' @param myPoly \code{SpatialCollections}, \code{SpatialPolygons}, \code{SpatialPoints}, or \code{SpatialLines} object
#' @return \code{SpatialPolygons} object
#' @importFrom methods is
#' @export
#' @examples
#' par(mfrow = c(1, 2))
#' plot(spatialCollEx, col = "blue", main = "Spatial Collections Object")
#' lineOnly <- keepLine(spatialCollEx)
#' plot(lineOnly, col = "blue", main = "Spatial Line Only")
keepLine <- function(myPoly) {
  if (is(myPoly)[1] == "SpatialCollections") {
    myPoly <- myPoly@lineobj
  } else if (is(myPoly)[1] == "SpatialPoints") {
    return(NULL)
  } else if (is(myPoly)[1] == "SpatialLines") {
    return(myPoly)
  } else if (is.null(myPoly)) {
    return(NULL)
  } else if (is(myPoly)[1] == "SpatialPolygons") {
    return(NULL)
  } else {
    stop("Incorrect object type")
  }
  return(myPoly)
}

#' Get coordinates from a spatial object of lines and points. There is no ordering of points returned.
#'  Note: This differs from \code{extractCoords} in that the ordering of the points is NOT considered.
#' @title Extract coordinates from a spatial object of lines and points
#' @param myPoints spatial object of type \code{SpatialCollections}, \code{SpatialPoints}, or \code{SpatialLines}
#' @return n x 2 matrix of coordinates
#' @importFrom rgeos gLineMerge
#' @importFrom sp disaggregate
#' @importFrom methods as is
#' @export
#' @examples
#' #Load sample line
#' exampleLine <- as(rmHoles(bgWater[2]), "SpatialLines")
#' getCoords(exampleLine)
getCoords <- function(myPoints) {
  ##split apart lines and points
  line <- points <- NULL
  first <- TRUE
  if (is(myPoints)[1] == "SpatialCollections") {
    line <- myPoints@lineobj
    if (!is.null(line)) {
      line <- disaggregate(gLineMerge(myPoints@lineobj))
    }
    points <- myPoints@pointobj
  } else if (is(myPoints)[1] == "SpatialPoints") {
    points <- myPoints
  } else if (is(myPoints)[1] == "SpatialLines") {
    line <- disaggregate(myPoints)
  }

  ##add coordinates from line objects
  if (!is.null(line)) {
    n1 <- length(line@lines)
    for (i in 1:n1) {
      temp <- line@lines[[i]]
      n2 <- length(temp)
      for (j in 1:n2) {
        if (first) {
          coords <- temp@Lines[[j]]@coords
          first <- FALSE
        } else {
          coords <- rbind(coords, temp@Lines[[j]]@coords)
        }
      }
    }
  }

  ##add individual points
  if (!is.null(points)) {
    if (first) {
      coords <- points@coords
      first <- FALSE
    } else {
      coords <- rbind(coords, points@coords)
    }
  }

  return(coords)
}

#' Function to find to which matrix indices coordinates correspond (on a 304 x 448 grid)
#' @title Find indices in matrix
#' @param coords coordinates of interest
#' @param xmn min x (defaults to value for Northern Polar stereographic grid: -3850)
#' @param ymn min y (defaults to value for Northern Polar stereographic grid: -5350)
#' @return n x 2 matrix of coordinates on a 304 x 448 grid
#' @export
#' @examples
#' dat <- matrix(nrow = 2, ncol = 2, data = c(-2000, 0, 300, 1000))
#' getInd(dat)
getInd <- function(coords, xmn = -3850, ymn = -5350) {
  cbind(round((coords[,1] - xmn)/25 + 1, 0),  round((coords[,2] - ymn)/25 + 1, 0))
}

#' Function to extract coordinates from a \code{SpatialLines} object.
#' If there are breaks in the line, this function connects the closest points to create one line.
#' Note: This differs from the function \code{getCoords} in that the ordering of the points is considered.
#' @title Function to extract coordinates.
#' @param x \code{SpatialLines} or \code{SpatialPolygons} object
#' @return n x 2 matrix of coordinates
#' @importFrom methods as is slot
#' @export
#' @examples
#' coords <- extractCoords(regionInfo$regions[[3]])
#' par(mfrow = c(1, 2))
#' plot(regionInfo$regions[[3]], main = "Polygon Object")
#' plot(coords, type = "p", main = "Coordinates", pch = 20)
extractCoords <- function(x){
  ##convert spatial polygon to spatial line object if needed
  if (is(x)[1] == "SpatialPolygons") {
    x <- as(x, "SpatialLines")
  }

  ##extract coordinates into a list of lines where each line is composed of a list of
  #segments
  res <- lapply(slot(x, "lines"), function(x) lapply(slot(x, "Lines"), function(y) slot(y, "coords")))

  ##go through each line segment and extract the coordinates in order
  nLines <- length(res)
  lines <- list()
  for (i in 1:nLines) {
    nSegs <- length(res[[i]])
    used <- rep(FALSE, nSegs)
    first <- lapply(res[[i]], function(x){x[1, ]}) #list of first coordinates in each segment
    last <- lapply(res[[i]], function(x){x[nrow(x),]}) #list of last coordinates in each segment
    coords <- res[[i]][[1]] #start coordinate matrix with first line
    used[1] <- TRUE
    while (any(!used)) {
      #find segments that start or end as close as possible to the beginning or
      #end of the existing coordinate matrix
      matchFirstDist <- unlist(lapply(last, function(x){getDist(x, coords[1,])}))
      matchLastDist <- unlist(lapply(first, function(x){getDist(x, coords[nrow(coords),])}))
      matchFirstDist[used] <- Inf; matchLastDist[used] <- Inf
      test <- which.min(c(min(matchFirstDist), min(matchLastDist)))
      if (test == 1) { #add segment to the start of the coordinate matrix
        addTop <- which.min(matchFirstDist)
        temp <- res[[i]][[addTop]]
        coords <- rbind(temp[1:(nrow(temp) - 1), ], coords)
        used[addTop] <- TRUE
      } else { #assign coordinate to the end of the coordinate matrix
        addBottom <- which.min(matchLastDist)
        temp <- res[[i]][[addBottom]]
        coords <- rbind(coords, temp[2:nrow(temp), ])
        used[addBottom] <- TRUE
      }
    }
    lines[[i]] <- coords
  }

  #if there is more than one line, put lines together into one line
  if (nLines > 1) {
    used <- rep(FALSE, nLines)
    first <- lapply(lines, function(x){x[1, ]}) #list of first coordinates in each segment
    last <- lapply(lines, function(x){x[nrow(x),]}) #list of last coordinates in each segment
    coords <- lines[[1]] #start coordinate matrix with first line
    used[1] <- TRUE
    while (any(!used)) {
      #find lines that start or end as close as possible to the beginning or
      #end of the existing coordinate matrix
      matchFirstDist <- unlist(lapply(last, function(x){getDist(x, coords[1,])}))
      matchLastDist <- unlist(lapply(first, function(x){getDist(x, coords[nrow(coords),])}))
      matchFirstDist[used] <- Inf; matchLastDist[used] <- Inf
      test <- which.min(c(min(matchFirstDist), min(matchLastDist)))
      if (test == 1) { #add line to the start of the coordinate matrix
        addTop <- which.min(matchFirstDist)
        temp <- lines[[addTop]]
        coords <- rbind(temp[1:(nrow(temp) - 1), ], coords)
        used[addTop] <- TRUE
      } else { #add line to the bottom of the coordinate matrix
        addBottom <- which.min(matchLastDist)
        temp <- lines[[addBottom]]
        coords <- rbind(coords, temp[2:nrow(temp), ])
        used[addBottom] <- TRUE
      }
    }
  }

  return(coords)
}

#' Takes in a matrix and returns polygons representing
#' contiguous regions in the matrix. Typically these regions are either where the ice concentration
#' is above a certain level or where there is land.
#' @title Get polygons corresponding to regions
#' @param dat matrix of one of the allowed data types ("gfdl", "bootstrap", or "simple)  (see details)
#' @param datType string indicating the format of the data: either "gfdl", "bootstrap", or "simple" (see details)
#' @param level concentration level of interest
#' @param myLandMat binary matrix specifying land locations
#' @param myAllRegions \code{SpatialPolygons} object specifying region that will be considered
#' @param useAll boolean, if true indicates to use the full area (overrides \code{myLandMat})
#' @param landInd boolean, if true indicates that the region of interest is the land
#' @param xmn min x (defaults to value for polar stereographic grid, -3850)
#' @param xmx max x (defaults to value for polar stereographic grid, 3750)
#' @param ymn min y (defaults to value for polar stereographic grid, -5350)
#' @param ymx max y (defaults to value for polar stereographic grix, 5850)
#' @details For \code{datType = "simple"}  the values in the \code{dat} matrix are
#' indicators of whether the grid box contains ice (1: ice-covered, 0: no ice, NA: land).
#' If \code{datType = "gfdl"} or \code{datType  = "bootstrap"}, the values in the matrix correspond
#' to the raw ice concentrations values observed or predicted (including indicators for missing data, land etc.).
#' If \code{datType = "gfdl"}, the predictions are
#' formatted as in the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model produced by the National Oceanic and Atmospheric Administration’s
#' Geophysical Fluid Dynamics Laboratory converted to a Polar Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014).
#' If \code{datType = "bootstrap"} the array values are formatted the same as the ice concentration values obtaned from the
#'  National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm.
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
#' @return region of interest as a \code{SpatialPolygons} object
#' @importFrom raster raster rasterToPolygons
#' @importFrom rgeos gIntersection
#' @importFrom sp disaggregate
#' @importFrom maptools unionSpatialPolygons
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' obsExample <-  getRegion(dat = obsFeb2012, datType = "bootstrap", level = 15)
#' plot(land, col = 'grey', border = FALSE)
#' plot(obsExample, col = "lightblue", add = TRUE)
#' }
getRegion <- function(dat,  datType, level, myLandMat = landMat, myAllRegions = allRegions,
                      useAll = FALSE, landInd = FALSE, xmn = -3850, xmx = 3750,
                      ymn = -5350, ymx = 5850) {
  #check input
  if (!(datType == "bootstrap" || datType == "gfdl" || datType == "simple")) {
    stop("datType must be 'bootstrap', 'gfdl', or 'simple'")
  }

  if (datType != "simple") {
    ##temporary matrix to store region before converting to a polygon
    temp <- matrix(nrow = nrow(dat), ncol = ncol(dat), data = NA)
  } else {
    temp <- dat
  }

  ##find region (if not land)
  if (datType == "gfdl" & landInd == FALSE) { #ice region, GFDL data
    temp[dat >= level/100] <- 1 #above ice concentration threshold
    temp[dat < level/100] <- 0 #below ice concentration threshold
    temp[is.na(dat)] <- 0 #land is not part of contiguous region
  } else if (datType == "bootstrap" & landInd == FALSE) { #ice region, bootstrap data
    dat[myLandMat == 1] <- 120 #convert specified land in landmat to land indicator for gfdl (120)
    temp[dat >= level & dat <= 100] <- 1 #above ice concentration threshold
    temp[dat == 110] <- 1 #assume satelite hole is ice
    temp[dat == 120] <- 0 #land is not part of the main contiguous area
    temp[dat < level] <- 0 #below ice concentration threshold
  ##Find region (if land)
  } else if (datType == "gfdl" & landInd == TRUE) {
    temp[dat < 1] <- 0 #set ice to zero
    temp[is.na(dat)] <- 1 #set land to zero
  } else if (datType == "bootstrap" & landInd == TRUE) {
    dat[myLandMat == 1] <- 120 #convert to same land boundaries as gfdl model
    temp[dat <= 100] <- 0 #ice and ocean are not land
    temp[dat == 110] <- 0 #assume satelite hole is not on land
    temp[dat == 120] <- 1 #land
  }

  ##convert to polygon object
  temp2 <- t(temp[,ncol(temp):1]) #orientation used by raster
  poly <- raster(temp2, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx) #make raster
  poly <- rasterToPolygons(poly, fun = function(x){x > 0}) #make set of grid boxes corresponding to polygon
  if (landInd == TRUE) { #polygon ID
    ID <- rep("land", length(poly))
  } else {
    ID <- rep("ice", length(poly))
  }
  polyUnion <- unionSpatialPolygons(poly, ID) #aggregate into one polygon with boundaries
  indiv <- disaggregate(polyUnion, hole = FALSE) #seperate into each polygon with boundaries

  ##finalize object
  if (!landInd & !useAll) {
    indiv <- gIntersection(indiv, myAllRegions)
  }
  if (is(indiv)[1] == "SpatialCollections") {
    indiv <- indiv@polyobj
  }

  return(indiv)
}

#' Determine initialization month based on month being forecast and lag. Considers lags up to 11 months in advance.
#' @title Get initialization month
#' @param month month being forecast (integer from 1 to 12)
#' @param lag months in advance prediction is being made (integer from 1 to 11).
#' @return integer corresponding to the initialization month
#' @details Note that this calculation assumes that the prediction for a month is its first day. This differs from the convention used in our
#' paper which rounds up to the nearest full month. In practice, this may not be the case.
#' @export
#' @examples
#' initMonth <- getInitMonth(month = 10, lag = 4)
#' initMonth
getInitMonth <- function(month, lag) {
  ##check "month" input
  if (!(month%%1== 0) || month < 1 || month > 12) {
    stop("month is not an integer from 1 to 12")
  }
  ##check "lag" input
  if (!(lag%%1  == 0) || lag < 1 || lag > 11) {
    stop("lag is not an integer from 1 to 12")
  }

  ##calculate initialization month
  initMonth <- month - lag
  if (initMonth < 1) {
    initMonth <- 12 - (-initMonth)
  }

  return(initMonth)
}


