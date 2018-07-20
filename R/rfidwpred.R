#' @title Generate spatial predictions using the hybrid method of random forest and
#' inverse distance weighting (RFIDW)
#'
#' @description This function is to make spatial predictions using the hybrid
#' method of random forest and inverse distance weighting (RFIDW).
#'
#' @param longlat a dataframe contains longitude and latitude of point
#' samples (i.e., trainx and trainy).
#' @param trainx a dataframe or matrix contains columns of predictive variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param longlatpredx a dataframe contains longitude and latitude of point locations
#' (i.e., the centres of grids) to be predicted.
#' @param predx a dataframe or matrix contains columns of predictive variables for
#' the grids to be predicted.
#' @param mtry a function of number of remaining predictor variables to use as
#' the mtry parameter in the randomForest call.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param idp numeric; specify the inverse distance weighting power.
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param ... other arguments passed on to randomForest or gstat.
#'
#' @return A dataframe of longitude, latitude and predictions.
#'
#' @references Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#' data(petrel.grid)
#' rfidwpred1 <- rfidwpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 3],
#' petrel.grid[, c(1,2)], petrel.grid, ntree = 500, idp = 2, nmax = 12)
#' names(rfidwpred1)
#' }
#'
#' @export
rfidwpred <- function (longlat, trainx, trainy, longlatpredx, predx, mtry =
  function(p) max(1, floor(sqrt(p))), ntree = 500, idp = 2, nmax = 12, ...) {
  names(longlat) <- c("LON", "LAT")
  names(longlatpredx) <- c("LON", "LAT")
  p <- ncol(trainx)
  rf.1 <- randomForest::randomForest(trainx, trainy, mtry = mtry(p), ntree =
    ntree)
  rf.pred <- stats::predict(rf.1, predx)
  data.dev <- longlat
  data.pred <- longlatpredx
  data.dev$var1 <- trainy - stats::predict(rf.1, trainx)
  gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev, set = list(idp = idp), nmax=nmax)
  idw.pred <- stats::predict(gstat1, data.pred)$var1.pred
  rfidw.pred1 <- rf.pred + idw.pred
  rfidw.pred <- cbind(longlatpredx, rfidw.pred1)
  names(rfidw.pred) <- c("LON", "LAT", "Predictions")
  rfidw.pred
}
