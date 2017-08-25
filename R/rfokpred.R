#' @title Generate spatial predictions using the hybrid method of random forest and ordinary kriging (RFOK)
#'
#' @description This function is to make spatial predictions using the hybrid
#' method of random forest and ordinary kriging (RFOK).
#'
#' @param longlat a dataframe contains longitude and latitude of validation
#' samples.
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
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12.
#' @param vgm.args arguments for vgm, e.g. variogram model of response
#' variable and anisotropy parameters. see notes vgm in gstat for details.
#' By default, "Sph" is used.
#' @param block block size. see krige in gstat for details.
#' @param ... other arguments passed on to randomForest or gstat.
#'
#' @return A dataframe of longitude, latitude and predictions.
#'
#' @note This function is largely based rfcv in randomForest.  When 'A zero or
#' negative range was fitted to variogram' occurs, to allow OK running, the
#' range was set to be positive by using min(vgm1$dist). In this case, caution
#' should be taken in applying this method, although sometimes it can still
#' outperform IDW and OK.
#'
#' @references Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#' data(petrel.grid)
#' rfokpred1 <- rfokpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 3],
#' petrel.grid[, c(1,2)], petrel.grid, ntree = 500, nmax = 12, vgm.args =
#' ("Sph"),)
#' names(rfokpred1)
#' }
#'
#' @export
rfokpred <- function (longlat, trainx, trainy, longlatpredx, predx, mtry =
  function(p) max(1, floor(sqrt(p))), ntree = 500, nmax = 12, vgm.args =
  ("Sph"), block = 0, ...) {
  names(longlat) <- c("LON", "LAT")
  names(longlatpredx) <- c("LON", "LAT")
  p <- ncol(trainx)
  rf.1 <- randomForest::randomForest(trainx, trainy, mtry = mtry(p), ntree =
    ntree)
  rf.pred <- stats::predict(rf.1, predx)
  data.dev <- longlat
  data.pred <- longlatpredx
  data.dev$var1 <- trainy - stats::predict(rf.1, trainx)
  sp::coordinates(data.dev) = ~ LON + LAT
  vgm1 <- gstat::variogram(var1 ~ 1, data.dev)
  model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args))
  if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to
    variogram. To allow gstat running, the range was set to be positive by
    using min(vgm1$dist). ", "\n"))
  if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
  sp::coordinates(data.pred) = ~ LON + LAT
  ok.pred <- gstat::krige(var1 ~ 1, data.dev, data.pred, model = model.1,
    nmax = nmax, block = block)$var1.pred
  # rfok predictions
  pred.rfok <- rf.pred + ok.pred
  rfok.pred <- cbind(longlatpredx, pred.rfok)
  names(rfok.pred) <- c("LON", "LAT", "Predictions")
  rfok.pred
}
