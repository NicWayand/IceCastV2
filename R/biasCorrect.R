 ###Functions to apply bias correction after mappings have been obtained###

#' Apply contour-shifting to bias correct a predicted contour using existing mappings.
#' @title Apply contour-shifting to bias correct
#' @param maps object obtained from the \code{createMappings} function (see details)
#' @param predicted array of predicted values of dimension year x month x longitude x latitude
#' @param bcYear year to be bias-corrected
#' @param predStartYear year prediction array starts in
#' @param regions region information list (see details)
#' @param level concentration level for which to build contour
#' @param datTypePred string indicating the format of the prediction: either "gfdl" or "simple" (see details)
#' @param myLandMat binary matrix specifying land locations
#' @param myAllRegions a single \code{SpatialPolygons} object given the region under consideration
#' @param myLand \code{SpatialPolygons} corresponding to the land
#' @return \code{SpatialPolygons} object of the adjusted region

#' @export
#' @details The object \code{maps} is obtained from running the \code{createMapping} function. It is a list of five objects where \code{month},
#' \code{startYear}, and \code{endYear} give the month, first year, and last year that were mapped. The variables \code{obsList} and \code{predList}
#' are lists of arrays with one 3-dimensional array for each region. The first dimension is for the year. The other two dimensions
#' are for the fixed points' y-coordinates, the mapped points' x-coordinates, the mapped points' y-coordinates, the length of the mapping vectors in the
#' x-direction, the length of the vectors in the y-direction, and the angles of the mapping vectors.
#'
#' A region information list is a list six objects: \code{regions}, \code{lines}, \code{out}, \code{centRegion}, \code{centLines}, and
#' \code{centFrom}. The first three objects are ordered lists giving information about each of the regions that will be mapped outside the
#' Central Arctic region. The variable \code{regions} gives \code{SpatialPolygons} objects for the corresponding regions and the variable \code{lines}
#'  gives the \code{SpatialLines} objects for the corresponding fixed lines. The variable \code{out} gives \code{SpatialPolygons} objects that are outside the corresponding regions,
#'  but that border the fixed lines. These are used when building new polygons to determine if points are being mapped outside the region of interest.
#'  The variable \code{centRegions} is the \code{SpatialPolygons} object corresponding to the central Arctic region, \code{centLines} is the \code{SpatialLines} object for the fixed line,
#'  and \code{centFrom} is an n x 2 matrix with each row repeatedly giving the coordinates of the center point from which mapping vectors will emanate.
#'  The package contains \code{regionInfo} in the \code{regionInfo.rda} file, which can be used unless you want to define your own regions.
#'
#' The predicted data array, \code{predicted}, should be single array of dimension: years x longitude (304) x latitude (448).
#' If \code{datTypePred = "simple"}, the values in the array should indicate whether each grid box is
#' categorized to contain ice (1: ice-covered, 0: no ice, NA: land). If \code{datTypePred = "gfdl"}.
#' #' If \code{datTypePred = "gfdl"} , the values in the \code{predicted} array correspond
#' to the raw ice concentrations values predicted (including indicators for missing data, land etc.)
#' formatted as in the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model produced by the National Oceanic and Atmospheric Administration’s
#' Geophysical Fluid Dynamics Laboratory converted to a Polar Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014). Weights for
#' converting to a polar stereograhic grid were obtained from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#'
#'
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
#'
#' @importFrom sp disaggregate Polygon
#' @importFrom rgeos gArea
#' @importFrom MASS rlm
#' @importFrom stats predict
#' @examples
#' \dontrun{
#' adj <- contourShift(maps = discrep, predicted = emFeb2012, bcYear = 2012, predStartYear = 1980,
#'                     regions = regionInfo, level = 15, datTypePred = "gfdl")
#' plot(land, col = "grey", border = FALSE)
#' plot(adj, add = TRUE, col = "blue")
#' }
contourShift <- function(maps, predicted, bcYear, predStartYear, regions, level, datTypePred, myLandMat = landMat,
                         myAllRegions = allRegions, myLand = land, endYearOffset) {

  ##Read-in and format prediction
  raw <- getRegion(dat = predicted, datType = datTypePred, level = level)
  discrepIndex <- bcYear - maps$startYear + 1 #need indices to match dimension of xDelta, yDelta, etc. set in getDiscrepPar.R script
  nReg <- length(regions$regions)

  ##find mappings for raw prediction outside central Arctic region
  rawMap <- getMap(ice = raw)

  ##find mapping for central Arctic region
  rawCent <- matrix(ncol = 4, nrow = nrow(regions$centFrom))
  rawCent[, 1:2] <- regions$centFrom
  rawCurr <- rmHoles(keepPoly(gIntersection(raw, regions$centRegion))) #ice in central Arctic region
  rawCurr <- disaggregate(rmHoles(rawCurr)) #holes ignored
  rawCurr <- rawCurr[which.max(gArea(rawCurr, byid = T))] #focus on largest region in central Arctic
  for (s in 1:length(regions$centLines)) {
    possRaw <- gIntersection(rawCurr, regions$centLines[[s]]) #intersection
    if (!is.null(possRaw)) { #get coordinates of intersection
      possRaw <- getCoords(possRaw)
      rawCent[s, 3:4] <- possRaw[which.max((possRaw[,1] - regions$centFrom[s, 1])^2 + (possRaw[,2] - regions$centFrom[s, 2])^2),]
    } else { #no intersection, store zero mapping
      rawCent[s, 3:4] <- regions$centFrom[s,]
    }
  }

  ##adjust typical regions
  yearInd <- maps$startYear:(bcYear - 1)
  end <- list()
  for (j in 1:nReg) {
    line <- regions$lines[[j]]@lines[[1]]@Lines[[1]]@coords
    new <- matrix(nrow = nrow(line), ncol = 2)
    colnames(new) <- c("x", "y")
    for (s in 1:nrow(line)) {
      #calculate pattern in x-values and bias correct
      xTempObs <- maps$obsList[[j]][1:(discrepIndex - endYearOffset), s, 3]
      lmXObs <- suppressWarnings(rlm(xTempObs ~ yearInd))
      xObsPred <- predict(lmXObs, newdata = data.frame(yearInd = bcYear))
      xTempRaw <- maps$predList[[j]][1:(discrepIndex - endYearOffset), s, 3]
      lmXRaw <- suppressWarnings(rlm(xTempRaw ~ yearInd))
      xRawPred <- predict(lmXRaw, newdata = data.frame(yearInd = bcYear))
      new[s, "x"] <- rawMap[[j]][s, "toX"] + (xObsPred - xRawPred)

      #calculate pattern in y-values and bias correct
      yTempObs <- maps$obsList[[j]][1:(discrepIndex - endYearOffset), s, 4]
      lmYObs <- suppressWarnings(rlm(yTempObs ~ yearInd))
      yObsPred <- predict(lmYObs, newdata = data.frame(yearInd = bcYear))
      yTempRaw <-  maps$predList[[j]][1:(discrepIndex - endYearOffset), s, 4]
      lmYRaw <- suppressWarnings(rlm(yTempRaw ~ yearInd))
      yRawPred <- predict(lmYRaw, newdata = data.frame(yearInd = bcYear))
      new[s, "y"] <- rawMap[[j]][s, "toY"]  + (yObsPred - yRawPred)
    }
    end[[j]] <- new
  }

  ##make adjusted points for each region into a polygon and combine polygons for all regions
  adj <- NULL
  for (j in 1:nReg) {
    new <- makePolygons(myEnd = end[[j]], myFixedLine = regions$lines[[j]],
                        myOut = regions$out[[j]], polyName = sprintf("new%i", j))
    if (is.null(adj) & !is.null(new)) {
      adj <- new
      adj@polygons[[1]]@ID <- "adj"
    } else if (!is.null(new)) {
      adj <- aggregate(spRbind(adj, new))
      adj@polygons[[1]]@ID <- "adj"
    }
  }

  ##add regions from the raw prediction that don't touch the land at all
  tempRaw <- disaggregate(keepPoly(gDifference(raw, regions$centRegion)))
  add <- tempRaw[!gIntersects(tempRaw, aggregate(myLand), byid = T)]
  adj <- aggregate(spRbind(adj, add))
  adj <- aggregate(gDifference(adj, myLand))
  adj@polygons[[1]]@ID <- "adj"

  ##adjust Central Arctic region
  new <- matrix(ncol = 2, nrow = length(regions$centLines))
  for (s in 1:length(regions$centLines)) {
    #calculate length of lines
    xTempRaw <- maps$predList[[nReg + 1]][1:(discrepIndex - endYearOffset), s, 3]
    yTempRaw <- maps$predList[[nReg + 1]][1:(discrepIndex - endYearOffset), s, 4]
    rawLength <- sqrt(xTempRaw^2 + yTempRaw^2)
    xTempObs <- maps$obsList[[nReg + 1]][1:(discrepIndex - endYearOffset), s, 3]
    yTempObs <- maps$obsList[[nReg + 1]][1:(discrepIndex - endYearOffset), s, 4]
    obsLength <- sqrt(xTempObs^2 + yTempObs^2)
    #use regression to get adjustment
    if (length(unique(rawLength)) > 1) {
      lmRaw <- suppressWarnings(rlm(rawLength ~ yearInd))
      rawPred <- predict(lmRaw, newdata = data.frame(yearInd = bcYear))
    } else {
      rawPred <- rawLength[1]
    }
    if (length(unique(obsLength)) > 1) {
      lmObs <- suppressWarnings(rlm(obsLength ~ yearInd))
      obsPred <- predict(lmObs, newdata = data.frame(yearInd = bcYear))
    } else {
      obsPred <- obsLength[1]
    }
    #calculate new length of mapping vectors
    xRaw <- rawCent[s, 3]
    yRaw <- rawCent[s, 4]
    rawDist <- sqrt(xRaw^2 + yRaw^2)
    adjDist <- rawDist + (obsPred - rawPred)
    R <- adjDist/rawDist #ratio for similar triangle
    new[s, 1] <- R*xRaw
    new[s, 2] <- R*yRaw
  }

  ##make polygon for central Arctic region
  centPoly <- SpatialPolygons(list(Polygons(list(Polygon(new, hole = FALSE)), "centerRegion")))
  if (!is.null(centPoly)) {
    if (!suppressWarnings(gIsValid(centPoly))) {
      centPoly <- untwist(centPoly)
      centPoly <- gDifference(centPoly, myLand)
    } else {
      centPoly <- gDifference(centPoly, myLand)
    }
    while (!suppressWarnings(gIsValid(centPoly))) {
      centPoly <- untwist(centPoly)
      centPoly <- gDifference(centPoly, myLand)
    }
    centPoly <- aggregate(centPoly)
    centPoly@polygons[[1]]@ID <- "centerRegion"
  }

  ##add any areas predicted by the physical model that do not touch the main region to prediction
  tempRaw <- disaggregate(keepPoly(gIntersection(raw, regions$centRegion)))
  add <- tempRaw[!suppressWarnings(gIntersects(tempRaw, centPoly, byid = T))]
  centPoly <- aggregate(spRbind(centPoly, add))

  #combine ice in central Arctic region with rest of regions
  adj <- aggregate(spRbind(adj, centPoly))

  return(adj)
}


#' Create a new polygon from the coordinates of mapped points
#' @title Create polygon from mapped points
#' @param myEnd n x 2  list of mapped points, i.e. the points to which the polygon should extend
#' @param myFixedLine \code{SpatialLines} object from which the mapping vectors emanate
#' @param myOut \code{SpatialPolygons} object outside the current region that borders the region's fixed line. It is
#' used to determine if points are being mapped outside the region of interest.
#' @param myLand \code{SpatialPolygons} object for the land
#' @param polyName character string to name the new polygon (defaults to "unspecified")
#' @return \code{SpatialPolygons} object created from the mapped points
#' @export
#' @importFrom rgeos gIntersects gIsValid gDifference
#' @importFrom raster aggregate
#' @importFrom sp SpatialPoints SpatialPolygons Polygon Polygons
#' @importFrom maptools spRbind
#' @examples
#' newPoly <- makePolygons(myEnd = mappedPoints, myFixedLine = regionInfo$lines[[1]],
#'                         myOut = regionInfo$out[[1]])
#' plot(newPoly)
makePolygons <- function(myEnd, myFixedLine, myOut, myLand = land, polyName = "unspecified") {
  ##find values of myEnd that go into the land indicating that a new polygon should be drawn
  touchInd <- which(gIntersects(SpatialPoints(myEnd), aggregate(myOut), byid = T))
  touchInd <- unique(c(1, touchInd)) #need to start building at the beginning of the will always need to start somewhere, so add first point
  nTouch <- length(touchInd)
  l1 <- myFixedLine@lines[[1]]@Lines[[1]]@coords
  new <- curr <- NULL
  myEnd <- round(myEnd) #round to the nearest planar coordinate for matching

  ##loop over values that touch myFixedLine or go into myOut and make polygons for each set of coordinates in between
  if (nTouch > 1) {
    for (s in 1:(nTouch - 1)) {
      #find point on myEnd
      firstEnd <- (touchInd[s] + 1)
      if (s != (nTouch - 1)) {
        lastEnd <- (touchInd[s + 1] - 1)
      } else {
        lastEnd <- nrow(myEnd)
      }
      #find corresponding points on myFixedLine
      firstFixed <- which.min((l1[, 1] - myEnd[firstEnd, 1])^2 +  (l1[, 2] - myEnd[firstEnd, 2])^2)
      lastFixed <- which.min((l1[, 1] - myEnd[lastEnd, 1])^2 +  (l1[, 2] - myEnd[lastEnd, 2])^2)
      #check if just making a line along the land, so we wouldn't need a polygon
      isPoly <- TRUE
      onEnd <- matrix(myEnd[firstEnd:lastEnd, ], ncol = 2)
      onFixed <-  matrix(l1[lastFixed:firstFixed, ], ncol = 2)
      if ((firstFixed == lastFixed) || (firstEnd == lastEnd)) {
        isPoly <- FALSE
      }
      if ((nrow(onEnd) == 2) & (nrow(onFixed) == 2)) {
        if (isTRUE(all.equal(onEnd, onFixed, 1))) {
          isPoly <- FALSE
        }
      }
      #make new polygon
      if (isPoly) {
        if (firstEnd != lastEnd + 1) {
          curr <- SpatialPolygons(list(Polygons(list(Polygon(rbind(onEnd, onFixed))), sprintf("new%i", s))))
        } else  {
          curr <- SpatialPolygons(list(Polygons(list(Polygon(onFixed)), sprintf("curr%i", s))))
        }
        #untwist polygon if needed or leave as just a horizontal or vectical line
        coords <- curr@polygons[[1]]@Polygons[[1]]@coords
        coords <- matrix(as.numeric(coords), ncol = 2)
        if ((length(unique(coords[, 1])) == 1) || (length(unique(coords[, 2])) == 1)) {
          curr <- NULL
        } else if (!suppressWarnings(gIsValid(curr))) {
          curr <- suppressWarnings(untwist(curr, polyName = sprintf("curr%i", s)))
        }
        #combine individual polygons into one polygon object
        if (is.null(new) & !is.null(curr)) {
          new <- curr
        } else if (!is.null(curr)) {
          new <- aggregate(spRbind(new, curr))
        }
      }
    }
  } else  { #myEnd touches myFixedLine only once, so make only one polygon
    stopifnot(nTouch == 1)
    revInd <- (sum((l1[1, ] - myEnd[1,])^2) < sum((l1[nrow(l1), ] - myEnd[1,])^2)) #check which way end points match up
    if (revInd) {
      new <-  SpatialPolygons(list(Polygons(list(Polygon(rbind(l1, myEnd[nrow(myEnd):1, ]))), "new")))
    } else {
      new <-  SpatialPolygons(list(Polygons(list(Polygon(rbind(l1, myEnd))), "new")))
    }
    if (!suppressWarnings(gIsValid(new))) {
      new <- suppressWarnings(untwist(new, polyName = "new1"))
    }
  }

  ##remove land & self-intersections and aggregate
  if (!is.null(new)) {
    if (!suppressWarnings(gIsValid(new))) {
      new <- untwist(new)
      new <- gDifference(new, myLand)
    } else {
      new <- gDifference(new, myLand)
    }
    while (!gIsValid(new)) {
      new <- untwist(new)
      new <- gDifference(new, myLand)
    }
    if (!is.null(new)) {
      new <- aggregate(new)
      new@polygons[[1]]@ID <- polyName
    }
  }

  return(new)
}

