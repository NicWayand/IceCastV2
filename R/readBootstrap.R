#' Read in individual binary files of monthly observation data. The observations are from the monthly sea ice concentration
#' obtained from the National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. The results
#' are distributed by the National Snow and Ice Data Center (Comiso 2000, updated 2015).
#' Functions assume file name conventions are the same as used by NSIDC.
#' @title Read individual bootstrap binary file
#' @references
#' Bootstrap sea ice concentration:
#'
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#' @param fileName  File name for binary bootstrap data
#' @param nX dimension in the x (defaults to value for Northern Polar stereographic grid: 304)
#' @param nY dimension in the y (defaults to value for Northern Polar stereographic grid: 448)
#' @return numeric vector of concentrations
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' #fileName should be the binary file
#' rawData <- readBootstrap(fileName)
#'}
readBootstrap <- function(fileName, nX = 304, nY = 448) {
  to.read <- file(fileName, "rb")
  dat <- readBin(to.read, integer(), n = nX*nY, size = 2, endian = "little")/10
  close(to.read)
  return(dat)
}

#' Function to process monthly bootstrap data over multiple years. The observations are from the monthly sea ice concentration
#' obtained from the National Aeronautics and Space Administration (NASA) satellites Nimbus-7
#' SMMR and DMSP SSM/I-SSMIS and processed by the bootstrap algorithm. The results
#' are distributed by the National Snow and Ice Data Center (Comiso 2000, updated 2015).
#' Functions assume file name conventions are the same as used by NSIDC.
#' @title Read in a set of bootstrap observations over a set of year
#' @references
#' Bootstrap sea ice concentration:
#'
#' Comiso, J., 2000, updated 2015: Bootstrap sea ice concentrations from Nimbus-7 SMMR and
#' DMSP SSM/I-SSMIS. version 2. \url{http://nsidc.org/data/nsidc-0079}
#' @param startYear first year to read in
#' @param endYear lastYear to read in
#' @param fileFolder Folder in which binary files are stored
#' @param version Either 2 or 3 indicating which version of the bootstrap data you are using
#' @param nX longitude dimension
#' @param nY latitude dimension
#' @details Raw binary files for 2012-2013 are included in the package
#' @export
#' @return Bootstrap observations sorted into array of dimension: year x month x lon x lat
#' @examples
#' \dontrun{
#' #myFilePath should be a file path where the 1983 binary files are stored
#' observedDemo <- readMonthlyBS(startYear = 1983, endYear = 1983, fileFolder = myFilePath)
#' }
readMonthlyBS <- function(startYear, endYear, fileFolder, version, nX = 304, nY = 448) {
  years <- startYear:endYear; nYears <- length(years)
  obs <- array(dim = c(nYears, 12, nX, nY))
  stopifnot(version == 2 || version == 3)
  for (i in 1:nYears) {
    for (j in 1:12) {
      if (version == 2) { #no missing data in V2
        fileName <- Sys.glob(paste(fileFolder, sprintf('bt_%i%02d_*_v02_n.bin', years[i], j), sep = ""))
        obs[i, j, ,nY:1] <- readBootstrap(fileName)
      } else if (version == 3 & !(j == 12 & years[i] == 1987) & !(j == 1 & years[i] == 1988)) {
        #missing Dec 1987 and Jan 1988 in V3
        #"major data gap in the SSM/I data occurs from 03 December 1987 to 13 January 1988"
        fileName <- Sys.glob(paste(fileFolder, sprintf('bt_%i%02d_*_v03_n.bin', years[i], j), sep = ""))
        obs[i, j, ,nY:1] <- readBootstrap(fileName)
      } else {
        stopifnot(version == 3)
        obs[i, j, ,nY:1] <- NA
      }
    }
  }
  return(obs)
}
