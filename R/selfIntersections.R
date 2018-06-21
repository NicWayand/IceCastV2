#' Find if two line segments intersect
#' @title Check if line segments intersect
#' @param a first coordinate of first line segment
#' @param b second coordinate of first line segment
#' @param c first coordinate of second line segment
#' @param d second coordinate of second line segment
#' @param seq indicator for whether the two line segments are intersecting
#' @return boolean indicating if there is an intersection
#' @export
#' @examples
#'  checkIntersect(c(0, 0), c(1, 1), c(2, 2), c(3, 3))
#'  checkIntersect(c(0, 0), c(1, 1), c(0.5, 0.5), c(2, 2))
checkIntersect <- function(a, b, c, d, seq = FALSE) {
  ##typical lines refers to lines that are not horizontal or vertical

  ##inequality functions that correctly account for floating point issues.
  ##(Care with equality of floating points is needed througout. This is why all.equal is used)
  lessEq <- function(a, b) {
    ((a < b) || isTRUE(all.equal(as.numeric(a), as.numeric(b))))
  }
  greaterEq <- function(a, b) {
    ((a > b) || isTRUE(all.equal(as.numeric(a), as.numeric(b))))
  }

  ##return true for identical points
  if (isTRUE(all.equal(a, c)) & isTRUE(all.equal(b, d))) {
    return(TRUE)

  ##two vertical lines
  } else if (isTRUE(all.equal(a[1], b[1])) & isTRUE(all.equal(c[1], d[1]))) {
    if (!isTRUE(all.equal(a[1], c[1]))) { #check if the two lines are not on the same x-coord
      return(FALSE)
      #if on same x-coord; y-coords need to overlap
    } else if ((max(a[2], b[2]) < min(c[2], d[2])) || (max(c[2], d[2]) < min(a[2], b[2]))) {
      return(FALSE) #if sequential equality means not an overlap
    } else if ((lessEq(max(a[2], b[2]), min(c[2], d[2])) || lessEq(max(c[2], d[2]), min(a[2], b[2]))) & seq == TRUE) {
      return(FALSE)
    } else {
      return(TRUE)
    }

  ##two horizontal lines
  } else if (isTRUE(all.equal(a[2], b[2])) & isTRUE(all.equal(c[2], d[2]))) {
    if (!(isTRUE(all.equal(a[2], c[2])))) { #check if the two lines are on the same y-coord
      return(FALSE)
      #if on same y-coord; x-coords need to overlap
    } else if ((max(a[1], b[1]) < min(c[1], d[1])) || (max(c[1], d[1]) < min(a[1], b[1]))) {
      return(FALSE) #if sequential equality means not an overlap
    } else if ((lessEq(max(a[1], b[1]), min(c[1], d[1])) || lessEq(max(c[1], d[1]), min(a[1], b[1]))) & seq == TRUE) {
      return(FALSE)
    } else {
      return(TRUE)
    }

 ##first line horizontal, second line vertical
  } else if (isTRUE(all.equal(a[2], b[2])) & isTRUE(all.equal(c[1], d[1]))) {
    if (greaterEq(c[1], min(a[1], b[1])) & lessEq(c[1], max(a[1], b[1])) &
        greaterEq(a[2], min(c[2], d[2])) & lessEq(a[2], max(c[2], d[2]))) {
      inter <- c(c[1], a[2])
      if (seq & isTRUE(all.equal(inter, b))) { #okay if line intersects only at end points for sequential coords
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }

  ##first line vertical, second line horizontal
  } else if (isTRUE(all.equal(c[2], d[2])) & isTRUE(all.equal(a[1], b[1]))) {
    if (greaterEq(a[1], min(c[1], d[1])) & lessEq(a[1], max(c[1], d[1])) &
        greaterEq(c[2], min(a[2], b[2])) & lessEq(c[2], max(a[2], b[2]))) {
      inter <- c(a[1], c[2])
      if (seq & isTRUE(all.equal(inter, d))) { #okay if line intersects only at end points for sequential coords
        return(FALSE)
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }

  ##first line horizontal, second typical
  } else if (isTRUE(all.equal(a[2], b[2]))) {
    #find typical line
    m2 <- (d[2] - c[2])/(d[1] - c[1])
    b2 <- c[2] - m2*c[1]
    #find possible intersection (ya = m2*xi + b2)
    xi <- (a[2] - b2)/m2
    yi <- m2*xi + b2
    #check if intersect point is in region; if so, return intersect point
    if (!(lessEq(xi, max(a[1], b[1])) & greaterEq(xi, min(a[1], b[1])) &
          lessEq(yi, max(a[2], b[2])) & greaterEq(yi, min(a[2], b[2])) &
          lessEq(xi, max(c[1], d[1])) & greaterEq(xi, min(c[1], d[1])) &
          lessEq(yi, max(c[2], d[2])) & greaterEq(yi, min(c[2], d[2])))) {
      return(FALSE)
    } else {
      inter <- c(xi, yi)
      if (seq & isTRUE(all.equal(inter, b))) { #okay if line intersects only at end points for sequential coords
        return(FALSE)
      } else {
        return(TRUE)
      }
    }

  ##first line typical, second horizontal
  } else if (isTRUE(all.equal(c[2], d[2]))) {
    #find typical line
    m2 <- (b[2] - a[2])/(b[1] - a[1])
    b2 <- a[2] - m2*a[1]
    #find possible intersection (ya = m2*xi + b2)
    xi <- (c[2] - b2)/m2
    yi <- m2*xi + b2
    #check if intersect point is in region; if so, return intersect point
    if (!(lessEq(xi, max(c[1], d[1])) & greaterEq(xi, min(c[1], d[1])) &
          lessEq(yi, max(c[2], d[2])) & greaterEq(yi, min(c[2], d[2])) &
          lessEq(xi, max(a[1], b[1])) & greaterEq(xi, min(a[1], b[1])) &
          lessEq(yi, max(a[2], b[2])) & greaterEq(yi, min(a[2], b[2])))) {
      return(FALSE)
    } else {
      inter <- c(xi, yi)
      if (seq & isTRUE(all.equal(inter, c))) { #okay if line intersects only at end points for sequential coords
        return(FALSE)
      } else {
        return(TRUE)
      }
    }

 ##first line typical, second vertical
  } else if (isTRUE(all.equal(c[1], d[1]))) {
    #find other line
    m1 <- (b[2] - a[2])/(b[1] - a[1])
    b1 <- a[2] - m1*a[1]
    #find possible intersection (yi = m1*xc + b1)
    yi <- m1*c[1] + b1
    #check if intersect point is in region; if so, return intersect point
    if (!(lessEq(yi, max(a[2], b[2])) & greaterEq(yi, min(a[2], b[2])) &
          lessEq(yi, max(c[2], d[2])) & greaterEq(yi, min(c[2], d[2])))) {
      return(FALSE)
    } else {
      inter <- c(c[1], yi)
      if (seq & isTRUE(all.equal(inter, b))) { #okay if line intersects only at end points for sequential coords
        return(FALSE)
      } else {
        return(TRUE)
      }
    }

  ##first line vertical, second line typical
  } else if (isTRUE(all.equal(a[1], b[1]))) {
    #find other line
    m2 <- (d[2] - c[2])/(d[1] - c[1])
    b2 <- c[2] - m2*c[1]
    #find possible intersection (yi = m2*xa + b2)
    yi <- m2*a[1] + b2
    #check if intersect point is in region; if so, return intersect point
    if (!(lessEq(yi, max(a[2], b[2])) & greaterEq(yi, min(a[2], b[2])) &
          lessEq(yi, max(c[2], d[2])) & greaterEq(yi, min(c[2], d[2])))) {
      return(FALSE)
    } else {
      inter <- c(a[1], yi)
      if (seq & isTRUE(all.equal(inter, b))) { #okay if line intersects only at end points for sequential coords
        return(FALSE)
      } else {
        return(TRUE)
      }
    }

  ##two typical lines
  } else {
    #find eq's of both lines
    m1 <- (b[2] - a[2])/(b[1] - a[1])
    m2 <- (d[2] - c[2])/(d[1] - c[1])
    b1 <- a[2] - m1*a[1]
    b2 <- c[2] - m2*c[1]

    #points are on the same line
    if (isTRUE(all.equal(m1, m2)) & isTRUE(all.equal(b1, b2))) {
      #points overlap when min(a, b) > max(c, d) or vice versa
      #It's sufficient to test x or y coords only, since coords are only the same line
      if (!((greaterEq(min(a[2], b[2]), max(c[2], d[2])) || greaterEq(min(c[2], d[2]), max(a[2], b[2]))))) {
        return(TRUE)
      } else if (isTRUE(all.equal(b, c)) & seq == FALSE) {
        return(TRUE)
      } else {
        return(FALSE)
      }

    } else {
      #find possible point of intersection (Yi = m1*xi + b1 & Yi = m2*xi + b2, solve for xi)
      xi <- (b2 - b1)/(m1 - m2)
      #check if possible intersect point is in region; if so, return intersect point
      if (!((lessEq(xi, max(a[1], b[1])) & greaterEq(xi, min(a[1], b[1]))) &
            (lessEq(xi, max(c[1], d[1])) & greaterEq(xi, min(c[1], d[1]))))) {
        return(FALSE)
      } else {
        inter <- c(xi, m1*xi + b1)
        if (seq & isTRUE(all.equal(inter, b))) { #okay if line intersects only at end points for sequential coords
          return(FALSE)
        } else {
          return(TRUE)
        }
      }
    }
  }
}

#' Determine if there are any intersecting line segments in a matrix of coordinates representing a line
#' @title Check if a line has intersecting segments
#' @param line matrix of coordinates corresponding to the line of interest
#' @return list where \code{list$any} is a boolean indicating if there are any intersections
#' and \code{list$val} is an index corresponding to the first intersection found
#' @importFrom utils combn
#' @export
#' @examples
#' checkResults <- anyIntersect(currSecEx)
#' checkResults$any #true/false
#' checkResults$val #indices of first intersection found
anyIntersect <- function(line) {
  ##line segments are referred to by the index of their first point
  nPoints <- nrow(line)

  ##only two points, no way to have an intersection
  if(nPoints == 2) {
    return(list("any" = FALSE, "val" = NA))
  }

  ##check all pairs of line segments for intersections
  pairs <- combn(1:(nPoints - 1), 2) #indices of all the combinations of ways pairs could be lined up
  for (i in 1:ncol(pairs)) { #loop through line seqment pairs, stop if an intersection is found
    currPair <- pairs[,i]
    seq <- ((currPair[2] - currPair[1]) == 1) #boolean to indicate if points are sequential (meaning they should have an intersection)
    if (checkIntersect(a = line[currPair[1], ], b = line[currPair[1] + 1,],
                       c = line[currPair[2],], d = line[currPair[2] + 1,], seq)) {
      #return indices of which line segments are the source of the intersection
      return(list("any" = TRUE, "val" = pairs[, i]))
    }
  }

  #return FALSE and NA if no intersections are found
  return(list("any" = FALSE, "val" = NA))
}


#' Function to remove all self-intersections from a contour.
#' @title Remove self-intersections
#' @param myPoly \code{SpatialPolygons} object from which self-intersections need to be removed
#' @param plotting boolean indicating if results should be plotted
#' @param polyName name for \code{SpatialPolygons} object to return (defaults to "unspecified")
#' @return \code{SpatialPolygons} object with self-intersections removed
#' @importFrom rgeos gIsValid
#' @importFrom sp Polygon
#' @export
#' @examples
#' \dontrun{
#' par(mfrow = c(1, 2))
#' plot(interEx, main = "Original Contour")
#' noInter <- untwist(interEx, polyName = "interEx")
#' plot(noInter, main = "Final Contour")
#' }
untwist <- function(myPoly, plotting = FALSE, polyName = "unspecified") {

  while (suppressWarnings(!gIsValid(myPoly))) {

    ##pull out coordinates
    coords <- myPoly@polygons[[1]]@Polygons[[1]]@coords
    coords <- matrix(as.numeric(coords), ncol = 2)

    ##return NULL if you start with just a horizontal or vertinal line
    if ((length(unique(coords[, 1])) == 1) || (length(unique(coords[, 2])) == 1)) {
      return (NULL)
    }

    ##find point at or near self intersection point
    prob <- suppressWarnings(gIsValid(myPoly, reason = T))
    prob <- as.numeric(unlist(strsplit(gsub("Self-intersection\\[|\\]", "", prob), " ")))
    prob <- prob[!is.na(prob)]

    ##find minimal region around intersection
    doesIntersect <- FALSE
    nNeigh <- 5
    currSec <- c()
    while (!doesIntersect) {
      index <- order(apply(coords, 1, function(x){  sqrt((x[1] - prob[1])^2 + (x[2] - prob[2])^2)}),
                     decreasing = FALSE)[1:nNeigh]
      index <- index[!is.na(index)] #remove NA's, which are put in place if more neightbors requested than points
      if (((nrow(coords) - max(index)) + min(index) <= max(index) - min(index)) & #Check if looping around end (BUT not needing all points)
          ((max(index) - min(index) + 1) !=  nrow(coords))) {
        cat <-  ((nrow(coords) - index) <= index) #categorize each point as closer to the end of the contour or the beginning
        curr <- c(min(index[cat == TRUE]):nrow(coords), 1:max(index[cat == FALSE]))
      } else { #standard case (no looping)
        curr <- min(index):max(index)
      }
      currSec <- coords[curr,]
      temp <- anyIntersect(currSec)
      if (temp$any) {

        ##if we've found the intersection, keep the minimal number of line segments
        ep <- rep(NA, 4)
        ep[1] <- curr[temp$val[1]] #first coord of intersecting line segment 1
        ep[3] <- curr[temp$val[2]] #first coord of intersecting line segment 2
        if (ep[1] != nrow(coords)) {#use max + 1, since line segments are indexed by their first point (checking for looping over)
          ep[2] <- ep[1] + 1
        } else {
          ep[2] <- 1
        }
        if (ep[3] != nrow(coords)) {
          ep[4] <- ep[3] + 1
        } else {
          ep[4] <- 1
        }
        if ((nrow(coords) - max(ep)) + min(ep) < max(ep) - min(ep) & (nrow(coords) > length(ep) + 1)) {
          #Check if indices loop back to beginning. The second section condition ensures that there is actually
          #a line segment to add back not just a point, in which case all points should be replaced
          cat <-  ((nrow(coords) - ep) <= ep) #categorize each point as closer to the end of the contour or the beginning
          curr <- c(min(ep[cat]):nrow(coords), 1:max(ep[!cat]))
        } else { #standard case (no looping)
          curr <- min(ep):max(ep)
        }
        currSec <- coords[curr, ]
        stopifnot(anyIntersect(currSec)$any) #you should always have an intersection at this point
        doesIntersect <- TRUE
      } else { #otherwise expand the region
        nNeigh <- nNeigh + 5
        if (nNeigh > nrow(coords)) {
          nNeigh <- nrow(coords)
        }
      }
    }
    new <- untwistSec(currSec)
    if (curr[1] > curr[length(curr)]) {
      #remove section being adjusted and put new points on the end
      cat <-  ((nrow(coords) - curr) <= curr) #categorize each point as closer to the end of the contour or the beginning
      coords <- rbind(coords[(max(curr[cat == FALSE]) + 1):(min((curr[cat == TRUE])) - 1),],  new)
    } else if (curr[length(curr)] == nrow(coords)) { #last point
      coords <- rbind(coords[1:(min(curr) - 1),], new)
    } else if (curr[1] == 1) { #check if indices start at one
      coords <- rbind(new, coords[(max(curr) + 1):nrow(coords),])
    } else { #no looping
      coords <- rbind(coords[1:(min(curr) - 1),], new, coords[(max(curr) + 1):nrow(coords), ])
    }

    ##if all the algorithm leads to is a single vertical line or horiztonal line,  assume no polygon
    if (length(unique(coords[,1])) == 1 || length(unique(coords[,2])) == 1) {
      return(NULL)
    }

    #If all the algorithm, leads to is a single set of two connected points, assume no polygon
    if (nrow(unique(round(coords, 2))) <= 2) {
      return(NULL)
    }

    ##make new myPoly and start over
    myPoly <- SpatialPolygons(list(
      Polygons(list(Polygon(coords)), ID = polyName))
    )
  }

  return(myPoly)
}

#' Function to correct self-intersections in a section of a line.
#' @title Remove self-intersections from one section of a contour
#' @param line  N x 2 matrix of coordinates
#' @param tol how much of a difference between the original line and the simplified line is allowed
#' @param eps how much to increase \code{tol} by on each iteration
#' @return n x 2 matrix of the new coordinates with self-intersections removed
#' @importFrom sp Line Lines SpatialLines
#' @importFrom rgeos gSimplify
#' @export
#' @examples
#' par(mfrow = c(1, 2))
#' plot(currSecEx, type = "l", main = "Original Line Section", xlab = "", ylab = "")
#' newSec <- untwistSec(currSecEx)
#' plot(newSec, type = "l", main = "New Line Section", xlab = "", ylab =  "")
untwistSec <- function(line, tol = 0, eps = .25) {

  l <- SpatialLines(list(Lines(Line(line), ID = "temp")))
  while (anyIntersect(line)$any) {
    tol <- tol + eps
    l <- gSimplify(l, tol)
    line <- l@lines[[1]]@Lines[[1]]@coords
  }

  return(line)
}
