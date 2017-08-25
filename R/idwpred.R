#' @title Generate spatial predictions using inverse distance weighting (IDW)
#'
#' @description This function is to make spatial predictions using inverse
#' distance weighting.
#'
#' @param longlat a dataframe contains longitude and latitude of samples.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in longlat.
#' @param longlat2 a dataframe contains longitude and latitude of point locations
#' (i.e., the centres of grids) to be predicted.
#' @param nmax for a local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param idp numeric; specify the inverse distance weighting power.
#' @param ... other arguments passed on to gstat.
#'
#' @return A dataframe of longitude, latitude and predictions.
#' @note This function is largely based on library(gstat).
#'
#' @references Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
#' package. Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' library(sp)
#' data(swmud)
#' data(sw)
#' idwpred1 <- idwpred(swmud[, c(1,2)], swmud[, 3], sw, nmax = 12, idp = 2)
#' names(idwpred1)
#' }
#'
#' @export
idwpred <- function (longlat, trainy, longlat2, nmax = 12, idp = 2, ...) {
  classRF <- is.factor(trainy)
  if (classRF) {
    stop ("This function is not for categorical response variable")
  }
  names(longlat) <- c("LON", "LAT")
  names(longlat2) <- c("LON", "LAT")
  data.dev <- longlat
  data.pred <- longlat2
  data.dev$var1 <- trainy
  gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev, set = list(idp = idp), nmax=nmax)
  idw.pred <- stats::predict(gstat1, data.pred)
  idw.pred
}
