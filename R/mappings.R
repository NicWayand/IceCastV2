###Functions for mapping##

#' Find the mapping vectors for one observation or prediction.
#' @title Map one observation or prediction
#' @param ice \code{SpatialPolygon} object corresponding to the region of ice
#' @param plotting boolean indicating if map should be plotted
#' @param reg  region information list (defaults to regionInfo$regions)
#' @param lines  lines (defaults to regionInfo$lines)
#' @param out back stop (defaults to regionInfo$out)
#' @param main character vector giving the title for the plotting map (defaults to no title)
#' @param myLand \code{SpatialPolygon} corresponding to the land
#' @param nSpace Spacing between points and lines that should be plotted (defaults to every seventh point and arrow)
#' @return List of the length of the number of regions. Each item in the list is a matrix. Each row of each matrix corresponds to a point in the region's fixed line. The seven columns give the fixed point's x-coordinate,
#' the fixed point's y-coordinate, the mapped point's x-coordinate, the mapped point's y-coordinate, the length of the mapping vectors in the
#' x-direction, the length of the vectors in the y-direction, and the angles of the mapping vectors.
#'
#' @details Often \code{reg}, \code{lines}, and \code{out} are taken from a region information list.
#' A region information list is a list of six objects: \code{regions}, \code{lines}, \code{out}, \code{centRegion}, \code{centLines}, and
#' \code{centFrom}. The first three objects are ordered lists giving information about each of the regions that will be mapped outside the
#' Central Arctic region. The variable \code{regions} gives \code{SpatialPolygons} objects for the corresponding regions and the variable \code{lines}
#'  gives the \code{SpatialLines} objects for the corresponding fixed lines. The variable \code{out} gives \code{SpatialPolygons} objects that are outside the corresponding regions,
#'  but that border the fixed lines. These are used when building new polygons to determine if points are being mapped outside the region of interest.
#'  The variable \code{centRegions} is the \code{SpatialPolygons} object corresponding to the central Arctic region, \code{centLines} is the \code{SpatialLines} object for the fixed line,
#'  and \code{centFrom} is an n x 2 matrix with each row repeatedly giving the coordinates of the center point from which mapping vectors will emanate.
#'  The package contains \code{regionInfo} in the \code{regionInfo.rda} file, which can be used unless you want to define your own regions.
#'
#' @importFrom sp SpatialPoints disaggregate
#' @importFrom rgeos gIntersection gIntersects gLineMerge
#' @importFrom maptools spRbind
#' @importFrom methods as
#' @importFrom graphics lines points
#' @importFrom sp plot
#' @export
#' @examples
#' \dontrun{
#' obs <- getRegion(dat = obsFeb19811982[1,,], datType = "bootstrap", level = 15)
#' obsMap <- getMap(ice = obs, plotting = TRUE,
#'                  main = "Observed Mapping \n February 1985")
#' }
getMap <- function(ice, plotting = FALSE, reg = regionInfo$regions, lines = regionInfo$lines,
                  out = regionInfo$out, main = "", myLand = land, nSpace = 7) {
  #colors from tim.colors(64) from the 'fields' package
  timColors <-  c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF",
                  "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF",
                  "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF",
                  "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF",
                  "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
                  "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00",
                  "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000",
                  "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000",
                  "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000",
                  "#800000")
  ##set up for angle legend by creating a square matrix with colors corresponding to angles
  if (plotting) {
    #angle from the center of the square
    x <- y <- seq(-1, 1, .01)
    grid <- expand.grid(x, y)
    angle <- apply(grid, 1, function(x){atan2(x[2], x[1])})
    legVals <- matrix(nrow = length(x), length(y), data = angle)
  }

  ##cycle through regions to get their mappings
  first <- TRUE
  nReg <- length(reg)
  mapList <- list()
  for (j in 1:nReg) {
    #set up
    SP <- as(lines[[j]], "SpatialPoints")
    points <- SP@coords
    map <- matrix(nrow = nrow(points), ncol = 7)
    colnames(map) <- c("fromX", "fromY", "toX", "toY", "deltaX", "deltaY", "angle")
    map[, c("fromX", "fromY")] <- points

    #get mapping for current region
    iceCurr <- keepPoly(gIntersection(ice, reg[[j]])) #ice in current region
    if (!is.null(iceCurr)) {
      inter <- suppressWarnings(gIntersects(disaggregate(iceCurr), lines[[j]], byid = T)) #ice intersects with fixed line
      if (any(inter)) {
        #clean up polygon
        iceCurr <- rmHoles(disaggregate(iceCurr)[inter])
        regCurr <- disaggregate(rmHoles(gIntersection(iceCurr, reg[[j]])))
        regCurr@polygons[[1]]@ID <- sprintf("regCurr%i",j)
        #loop over each section of ice
        for (m in 1:length(regCurr)) {
          temp <- findHoles(myPoly = regCurr[m], myReg = reg[[j]], myOut = out[[j]])
          ind <- which(gIntersects(temp$newPoly, SP, byid = T))
          lineCoord <- keepLine(gDifference(as(temp$newPoly, "SpatialLines"), lines[[j]]))
          if (!is.null(out) & length(ind) > 1) {
            #interpolate the points on the line to get the mapped points
            lineCoord <- extractCoords(gLineMerge(aggregate(lineCoord)))
            map[ind, c("toX", "toY")] <- intLine(points[ind,], lineCoord)
          }
        }

        #summary
        map[which(is.na(map[, "toX"])), "toX"] <- map[which(is.na(map[, "toX"])), "fromX"]
        map[which(is.na(map[, "toY"])), "toY"] <- map[which(is.na(map[, "toY"])), "fromY"]
        map[, "deltaX"] <- map[, "toX"] - map[, "fromX"]
        map[, "deltaY"] <- map[, "toY"] - map[, "fromY"]
        map[, "angle"] <- apply(map[, c("deltaX", "deltaY")], 1, function(x){atan2(x[2], x[1])})

        ##optional plotting of mapping
        if (first & plotting) {
          plot(myLand, col = "grey", main = main)
          if (requireNamespace("fields", quietly = TRUE)) {
            fields::add.image(3000, 5100, legVals, image.width = .1, image.height = .1)
          } else {
            print("Note: Fields Package is needed for adding angle legend to graph.Angle legend omitted")
          }
          first <- FALSE
        }
        if (plotting) { #plotting every nth point on the fixed line and mapping lines
          bp <- seq(-pi, pi, length.out = 65) #divide up (-pi, pi) into even break points
          angCol <- timColors[as.numeric(cut(map[, "angle"], breaks = bp))]
          plotInd <- 1:(nrow(points))%%7 == 0
          plot(iceCurr, lwd = 2, add = T)
          points(map[,c("fromX", "fromY")][plotInd, ], col = "black", pch = 20, cex = .25)
          #Note: Appropriate-sized arrows can be used in place of lines with the 'Arrows' function in the  "shape" package
          sapply((1:nrow(map))[plotInd], function(s){
            lines(x = c(map[s,"fromX"], map[s,"toX"]), y = c(map[s,"fromY"], map[s,"toY"]),
                  col = angCol[s], lwd = .7)})
        }
      } else { #otherwise use zero-vectors
        map[, c("toX", "toY")] <- map[, c("fromX", "fromY")]
        map[, c("deltaX", "deltaY")] <- 0
      }
    } else { #otherwise use zero-vectors
      map[, c("toX", "toY")] <- map[, c("fromX", "fromY")]
      map[, c("deltaX", "deltaY")] <- 0
    }
    mapList[[j]] <- map
  }
  return(mapList)
}

#' The function evenly spaces the number of points that are on one line, \code{predL}, on a
#' different line, \code{obsL}
#' @title Space points along a line
#' @param predL predicted line (n1 x 2 matrix of coordinates)
#' @param obsL predicted line (n2 x 2 matrix of coordinates)
#' @param plotting boolean indicating whether maps should be plotted
#' @return n x 2 matrix of evenly-spaced coordinates
#' @export
#' @examples
#' lineSpace <- intLine(predLEx, obsLEx, plotting = TRUE)
intLine <- function(predL, obsL, plotting = FALSE) {
  nPred <- nrow(predL); nObs <- nrow(obsL)  #number of points in prediction and observation
  #matrix giving where new points should be mapped to
  map <- matrix(nrow = nPred, ncol = 2, data = NA)
  colnames(map) <- c("x", "y")

  ##flip indexing of one line if needed
  if (getDist(obsL[1, ], predL[1, ]) > getDist(obsL[1, ], predL[nrow(predL),])) {
    obsL <- obsL[nrow(obsL):1, ]
  }

  ##identify where points on prediction should map onto observation
  if (nPred == nObs) { #same number of observation as prediction points, just map 1-to-1
    map <- obsL
  } else  {
    if ((nObs < nPred) || (nPred != 2)) { #typical case
      map[1, ] <- predL[1,]; map[nPred,] <- predL[nPred,] #assign first and last points to stay matching
      step <- (nObs - 1)/(nPred - 1) #proportion of sections of the line used per mapping
      used <- step #keep track of what proportion of the line has been used
      for (i in 2:(nPred - 1)) {
        s <- floor(used); e  <- ceiling(used) #indices of what points the dividing line for each section is between
        #assign mappings, indexing of matrices is off by one from indices of line
        map[i, 1] <- obsL[s + 1, 1] +  (used - s)*(obsL[e + 1, 1] - obsL[s + 1, 1])
        map[i, 2] <- obsL[s + 1, 2] +  (used - s)*(obsL[e + 1, 2] - obsL[s + 1, 2])
        used <- used + step #update portion of line used
      }
    } else { #have to be careful when the two
      step <- (nObs - 1)/3 #proportion of sections of the line used per mapping
      used <- step #set initial proportion used
      for (i in 1:nPred) {
        s <- floor(used); e  <- ceiling(used) #indices of what points the dividing line for each section is between
        #assign mappings, indexing of matrices is off by one from indices of line
        map[i, 1] <- obsL[s + 1, 1] +  (used - s)*(obsL[e + 1, 1] - obsL[s + 1, 1])
        map[i, 2] <- obsL[s + 1, 2] +  (used - s)*(obsL[e + 1, 2] - obsL[s + 1, 2])
        used <- used + step #update portion of line used
      }
    }
  }

  ##optional plotting of maps
  if (plotting) {
    plot(predL, col = "blue", xlim = c(min(c(obsL[,1], predL[,1])), max(c(obsL[,1], predL[,1]))),
         ylim = c(min(c(obsL[,2], predL[,2])), max(c(obsL[,2], predL[,2]))))
    points(predL, type = "l", col = "blue")
    points(obsL, col = "black")
    points(obsL, col = "black", type = "l")
    points(map, col = "red")
  }
  return(map)
}


#' Function to find and remove holes in a \code{SpatialPolygons} object. Note that this function differs from the function \code{rmHoles} in that it is considering gaps between a polygon and region boundaries in addition to holes contained within the polygon itself.
#' @title Find holes in a polygon
#' @param myPoly \code{SpatialPolygons} object of interest
#' @param myReg \code{SpatialPolygons} object for the region which myPoly is in
#' @param myOut  \code{SpatialPolygons} object that is outside the region, but that borders the regions' fixed line. This polygon is used to determine if points are being mapped outside the region of interest.
#' @return  list where \code{new$holes} gives the holes in the polygons as a \code{SpatialPolygons} object and \code{new$newPoly}
#' gives a \code{SpatialPolygons} object with the holes removed.
#' @importFrom raster aggregate
#' @importFrom sp SpatialPolygons Polygon Polygons
#' @importFrom maptools spRbind
#' @importFrom rgeos gIntersects
#' @export
#' @examples
#' holeEx <- findHoles(regionInfo$regions[[1]], regionInfo$regions[[1]], regionInfo$out[[1]])
#' plot(holeEx$newPoly)
#' plot(holeEx$holes, col = "red", add = TRUE)
findHoles <- function(myPoly, myReg, myOut) {
  holesKeep1 <- holesKeep2 <- NULL
  myPoly@polygons[[1]]@ID <- "myPoly" #re-name polygon

  ##find holes in the gap between ice and land
  if (!is.null(myOut)) {
    comb <- aggregate(spRbind(aggregate(myOut), myPoly))
    holes <- comb@polygons[[1]]@Polygons[which(sapply(comb@polygons[[1]]@Polygons, function(x){x@hole}))]
    if (length(holes) > 0) {
      poss <- SpatialPolygons(list(Polygons(list(Polygon(holes[[1]]@coords)), "poss")))
      if (length(holes) > 1) {
        for (k in 2:length(holes)) {
          poss <- spRbind(poss, SpatialPolygons(list(Polygons(list(Polygon(holes[[k]]@coords)), sprintf("poss%i", k)))))
        }
      }
      possInd <- which(gIntersects(poss, aggregate(myPoly), byid = T) & gIntersects(poss, myReg, byid = T))
      if (length(possInd) > 0) {
        holesKeep1 <- poss[possInd]
      }
    }
  }

  ##find holes in the center of the ice
  holes <- myPoly@polygons[[1]]@Polygons[which(sapply(myPoly@polygons[[1]]@Polygons, function(x){x@hole}))]
  if (length(holes) > 0) {
    poss <- SpatialPolygons(list(Polygons(list(Polygon(holes[[1]]@coords)), "holes")))
    if (length(holes) > 1) {
      for (k in 2:length(holes)) {
        poss <- spRbind(poss, SpatialPolygons(list(Polygons(list(Polygon(holes[[k]]@coords)), sprintf("poss%i", k)))))
      }
    }
    possInd <- which(gIntersects(poss, myPoly, byid = T))
    if (length(possInd) > 0) {
      holesKeep2 <- poss[possInd]
    }
  }

  ##combine polygons from both cases
  if (!is.null(holesKeep1) & !is.null(holesKeep2)) {
    holesKeep <- spRbind(holesKeep1, holesKeep2)
  } else if (!is.null(holesKeep1)) {
    holesKeep <- holesKeep1
  } else if (!is.null(holesKeep2)) {
    holesKeep <- holesKeep2
  } else {
    holesKeep <- NULL
  }

  ##remove holes to create new polygons without
  if (!is.null(holesKeep1)) {
    newPoly <- rmHoles(aggregate(spRbind(myPoly, holesKeep1)))
  } else {
    newPoly <- rmHoles(aggregate(myPoly))
  }

  new <- list("holes" = holesKeep, "newPoly" = newPoly)
  return(new)
}


#' Finds all the mappings for a set of observations and predictions often over multiple years
#' @title Map a set of observations and predictions
#' @param startYear first year to be mapped
#' @param endYear last year to be mapped
#' @param obsStartYear year in which observation array starts
#' @param predStartYear year in which prediction array starts
#' @param observed array of observed values of dimension year x longitude x latitude
#' @param predicted array of predicted values of dimension year x longitude x latitude
#' @param regions region information list (see detail)
#' @param month month under consideration
#' @param level concentration level for which to build contour
#' @param datTypeObs string of either "bootstrap" or "simple" indicating the file type of the observation (see details)
#' @param datTypePred string of either "gfdl" or "simple" indicating the file type of the prediction (see details)
#' @param plotting boolean indicatng whether maps should be plotted (defaults to false)
#' @param obsOnly run mapping only for observations
#' @param predOnly run mapping only for predictions
#' @param nX dimension in the x (defaults to value for Northern Polar stereographic grid: 304)
#' @param nY dimension in the y (defaults to value for Northern Polar stereographic grid: 448)
#' @param xmn min x value (defaults to value for Northern Polar stereographic grid: -3850)
#' @param xmx max x value (defaults to value for Northern Polar stereographic grid: 3750)
#' @param ymn min y value (defaults to value for Northern Polar stereographic grid: -5350)
#' @param ymx max y value (defaults to value for Northern Polar stereographic grid: 5850)
#' @return \code{map} object (see details)
#' @details
#' The object \code{maps} is obtained from running the \code{createMapping} function. It is a list of five objects where
#' \code{startYear} and \code{endYear} give first year, and last year that were mapped. The variables \code{obsList} and \code{predList}
#' are lists of arrays with one 3-dimensional array for each region. The first dimension is for the year. The other two dimensions
#' are for the fixed points' y-coordinates, the mapped points' x-coordinates, the mapped points' y-coordinates, the length of the mapping vectors in the
#' x-direction, the length of the vectors in the y-direction, and the angles of the mapping vectors.
#'
#'  A region information list is a list six objects: \code{regions}, \code{lines}, \code{out}, \code{centRegion}, \code{centLines}, and
#' \code{centFrom}. The first three objects are ordered lists giving information about each of the regions that will be mapped outside the
#' Central Arctic region. The variable \code{regions} gives \code{SpatialPolygons} objects for the corresponding regions and the variable \code{lines}
#'  gives the \code{SpatialLines} objects for the corresponding fixed lines. The variable \code{out} gives \code{SpatialPolygons} objects that are outside the corresponding regions,
#'  but that border the fixed lines. These are used when building new polygons to determine if points are being mapped outside the region of interest.
#'  The variable \code{centRegions} is the \code{SpatialPolygons} object corresponding to the central Arctic region, \code{centLines} is the \code{SpatialLines} object for the fixed line,
#'  and \code{centFrom} is an n x 2 matrix with each row repeatedly giving the coordinates of the center point from which mapping vectors will emanate.
#'  The package contains \code{regionInfo} in the \code{regionInfo.rda} file, which can be used unless you want to define your own regions.
#'
#' For \code{datTypeObs = "simple"} and \code{datTypePred = "simple"} the values in the \code{observed} and \code{predicted} arrays are
#' indicators of whether the grid box contains ice (1: ice-covered, 0: no ice, NA: land).
#' If \code{datTypePred = "gfdl"} or \code{datTypeObs  = "bootstrap"}, the values in the \code{observed} and \code{predicted} arrays correspond
#' to the raw ice concentrations values observed or predicted (including indicators for missing data, land etc.).
#' If \code{datTypePred = "gfdl"}, the predictions are
#' formatted as in the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model produced by the National Oceanic and Atmospheric Administration’s
#' Geophysical Fluid Dynamics Laboratory and converted to a Polar Stereographic grid (Vecchi et al. 2014; Msadek et al. 2014).
#' If \code{datTypeObs = "bootstrap"} the array values are assumed to be from the monthly sea ice concentration
#' obtained from the National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. Weights for
#' converting to a polar stereograhic grid were obtained from the spherical coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#'
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
#'
#' @importFrom sp disaggregate
#' @importFrom rgeos gArea gIntersection
#' @importFrom methods as
#' @importFrom graphics lines
#' @export
#' @examples
#' \dontrun{
#' createMapping(startYear = 1981, endYear = 1981, obsStartYear = 1981, predStartYear = 1980,
#'              observed = obsFeb19811982, predicted = emFeb19811982,
#'              regions = regionInfo, month = 2, level = 15,
#'              datTypeObs = "bootstrap", datTypePred = "gfdl", plotting = TRUE)
#' }
createMapping <- function(startYear, endYear, obsStartYear, predStartYear,
                          observed, predicted, regions, month, level,
                          datTypeObs, datTypePred, plotting = FALSE, obsOnly = FALSE,
                          predOnly = FALSE, nX = 304, nY = 448, xmn = -3850,
                          xmx = 3750, ymn = -5350, ymx = 5850) {
  #colors from tim.colors(64) from the 'fields' package
  timColors <-  c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", "#0000DF", "#0000EF",
                  "#0000FF", "#0010FF", "#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF",
                  "#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", "#00CFFF",
                  "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF",
                  "#50FFAF", "#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
                  "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", "#FFEF00", "#FFDF00",
                  "#FFCF00", "#FFBF00", "#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000",
                  "#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", "#FF0000",
                  "#EF0000", "#DF0000", "#CF0000", "#BF0000", "#AF0000", "#9F0000", "#8F0000",
                  "#800000")
  ##set up for angle legend by creating a square matrix with colors corresponding to angles
  if (plotting) {
    #angle from the center of the square
    x <- y <- seq(-1, 1, .01)
    grid <- expand.grid(x, y)
    angle <- apply(grid, 1, function(x){atan2(x[2], x[1])})
    legVals <- matrix(nrow = length(x), length(y), data = angle)
  }
  monthLab <- c("January", "February", "March", "April", "May", "June",
                "July", "August", "September", "October", "November", "December")
  nReg <- length(regions$regions)
  nYears <- length(startYear:endYear)
  predList <- obsList <- list()
  for (j in 1:nReg) {
    SP <- as(regions$lines[[j]], "SpatialPoints")@coords #need as many rows as points on the fixed line for each region
    predList[[j]] <- obsList[[j]]  <- array(dim = c(nYears, nrow(SP), 7))
  }
  obsList[[nReg + 1]] <- predList[[nReg + 1]] <- array(dim = c(nYears, nrow(regions$centFrom), 7))

  ##store mappings for all years
  for (i in startYear:endYear) {
    #read-in and format observation and prediction
    if (!predOnly) {
      obs <- getRegion(dat = observed[i - obsStartYear + 1, ,], datType = datTypeObs, level = level)
    }
    if (!obsOnly) {
      raw <- getRegion(dat = predicted[i - predStartYear + 1, ,], datType = datTypePred, level = level)
    }

    #map typical regions
    if (!predOnly) {
      obsMap <- getMap(ice = obs, plotting,
                       main = sprintf("Observed Mapping:\n %s %i", monthLab[month], i))
    }
    if (!obsOnly) {
    rawMap <- getMap(ice = raw, plotting,
                       main = sprintf("Predicted Mapping:\n %s %i", monthLab[month], i))
    }
    index <- i - startYear + 1
    for (j in 1:nReg) {
      if (!predOnly) {
        obsList[[j]][index,,] <- obsMap[[j]]
      }
      if (!obsOnly) {
        predList[[j]][index,,] <- rawMap[[j]]
      }
    }

    ##map observation in central Arctic region
    if (!predOnly) {
      obsList[[nReg + 1]][index, , 1:2] <- regions$centFrom #From
      obsCurr <- rmHoles(keepPoly(gIntersection(obs, regions$centRegion)))
      obsCurr <- disaggregate(rmHoles(obsCurr))
      obsCurr <- obsCurr[which.max(gArea(obsCurr, byid = T))]
      for (s in 1:length(regions$centLines)) {
        possObs <- gIntersection(obsCurr, regions$centLines[[s]])
        if (!is.null(possObs)) {
          possObs <- getCoords(possObs)
          obsList[[nReg + 1]][index, s, 3:4] <-
            possObs[which.max((possObs[,1] - regions$centFrom[s, 1])^2 + (possObs[,2] - regions$centFrom[s, 2])^2),]
        } else {
          obsList[[nReg + 1]][index, s, 3:4] <- regions$centFrom[s,]
        }
      }
      obsList[[nReg + 1]][index, , 5] <- obsList[[nReg + 1]][index, , 3] - obsList[[nReg + 1]][index, , 1]
      obsList[[nReg + 1]][index, , 6] <- obsList[[nReg + 1]][index, , 4] - obsList[[nReg + 1]][index, , 2]
      obsList[[nReg + 1]][index, , 7] <- apply(obsList[[nReg + 1]][index, , 5:6], 1, function(x){atan2(x[2], x[1])})

      ##optional plotting of observed central Arctic region mapping
      if (plotting) {
        mapCurr <- obsList[[nReg + 1]][index, ,]
        plotInd <- 1:(nrow(mapCurr))%%7 == 0
        plot(land, col = "grey", main = sprintf("Observed Mapping Central Arctic:\n %s %i", monthLab[month], i))
        plot(obsCurr, lwd = 2, add = T)
        if (requireNamespace("fields", quietly = TRUE)) {
          fields::add.image(3000, 5100, legVals, image.width = .1, image.height = .1)
        } else {
          print("Note: Fields Package is needed for adding angle legend to graph. Angle legend omitted")
        }
        bp <- seq(-pi, pi, length.out = 65) #divide up (-pi, pi) into even break points
        angCol <- timColors[as.numeric(cut(mapCurr[, 7], breaks = bp))]
        sapply((1:nrow(mapCurr))[plotInd], function(s){
          lines(x = c(mapCurr[s, 1], mapCurr[s, 3]), y = c(mapCurr[s, 2], mapCurr[s,4]),
                col = angCol[s], lwd = .7)})
      }
    }

    if (!obsOnly) {
      ##map prediction in central Arctic region
      predList[[nReg + 1]][index, , 1:2] <- regions$centFrom
      rawCurr <- rmHoles(keepPoly(gIntersection(raw, regions$centRegion)))
      rawCurr <- disaggregate(rmHoles(rawCurr))
      rawCurr <- rawCurr[which.max(gArea(rawCurr, byid = T))]
      for (s in 1:length(regions$centLines)) {
        possRaw <- gIntersection(rawCurr, regions$centLines[[s]])
        if (!is.null(possRaw)) {
          possRaw <- getCoords(possRaw)
          predList[[nReg + 1]][index, s, 3:4] <-
            possRaw[which.max((possRaw[,1] - regions$centFrom[s, 1])^2 + (possRaw[,2] - regions$centFrom[s, 2])^2),]
        } else {
          predList[[nReg + 1]][index, s, 3:4] <- regions$centFrom[s,]
        }
      }
      predList[[nReg + 1]][index, , 5] <- predList[[nReg + 1]][index, , 3] - predList[[nReg + 1]][index, , 1]
      predList[[nReg + 1]][index, , 6] <- predList[[nReg + 1]][index, , 4] - predList[[nReg + 1]][index, , 2]
      predList[[nReg + 1]][index, , 7] <- apply(predList[[nReg + 1]][index, , 5:6], 1, function(x){atan2(x[2], x[1])})

      ##optional plotting of predicted central Arctic region mapping
      if (plotting) {
        mapCurr <- predList[[nReg + 1]][index, ,]
        plot(land, col = "grey", main = sprintf("Predicted Mapping Central Arctic:\n %s %i", monthLab[month], i))
        plot(rawCurr, add = T, lwd = 2)
        bp <- seq(-pi, pi, length.out = 65) #divide up (-pi, pi) into even break points
        angCol <- timColors[as.numeric(cut(mapCurr[, 7], breaks = bp))]
        sapply((1:nrow(mapCurr))[plotInd], function(s){
          lines(x = c(mapCurr[s, 1], mapCurr[s, 3]), y = c(mapCurr[s, 2], mapCurr[s,4]),
                col = angCol[s], lwd = .7)})
        if (requireNamespace("fields", quietly = TRUE)) {
          fields::add.image(3000, 5100, legVals, image.width = .1, image.height = .1)
        } else {
          print("Note: Fields Package is needed for adding angle legend to graph. Angle legend omitted")
        }
      }
    }
    print(sprintf("Mapping for year %i complete", i))
  }

  return(list("startYear" = startYear, "endYear" = endYear,
              "predList" = predList, "obsList" = obsList))
}


