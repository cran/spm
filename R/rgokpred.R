#' @title Generate spatial predictions using the hybrid method of random forest in
#' ranger and ordinary kriging (RGOK)
#'
#' @description This function is to make spatial predictions using the hybrid
#' method of random forest in ranger and ordinary kriging (RGOK).
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
#' @param num.trees number of trees. By default, 500 is used.
#' @param min.node.size Default 1 for classification, 5 for regression.
#' @param type Type of prediction. One of 'response', 'se', 'terminalNodes' with
#' default 'response'. See ranger::predict.ranger for details.
#' @param num.threads number of threads. Default is number of CPUs available.
#' @param verbose Show computation status and estimated runtime.Default is FALSE.
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12.
#' @param vgm.args arguments for vgm, e.g. variogram model of response
#' variable and anisotropy parameters. see notes vgm in gstat for details.
#' By default, "Sph" is used.
#' @param block block size. see krige in gstat for details.
#' @param ... other arguments passed on to randomForest or gstat.
#'
#' @return A dataframe of longitude, latitude, predictions and variances. The
#' variances are produced by OK based on the residuals of rf.
#'
#' @note This function is largely based rfokpred.  When 'A zero or
#' negative range was fitted to variogram' occurs, to allow OK running, the
#' range was set to be positive by using min(vgm1$dist). In this case, caution
#' should be taken in applying this method, although sometimes it can still
#' outperform IDW and OK.
#'
#' @references Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation
#' of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17.
#' http://dx.doi.org/10.18637/jss.v077.i01.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#' data(petrel.grid)
#' rgokpred1 <- rgokpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 3],
#' petrel.grid[, c(1,2)], petrel.grid, num.trees = 500, nmax = 12, vgm.args =
#' ("Sph"))
#' names(rgokpred1)
#' }
#'
#' @export
rgokpred <- function (longlat, trainx, trainy, longlatpredx, predx, mtry = function(p)
                      max(1, floor(sqrt(p))), num.trees = 500, min.node.size = NULL,
                      type = "response", num.threads = NULL, verbose = FALSE, nmax = 12,
                      vgm.args = ("Sph"), block = 0, ...) {
  names(longlat) <- c("LON", "LAT")
  names(longlatpredx) <- c("LON", "LAT")
  p <- ncol(trainx)
  data.dev <- trainx
  data.dev$var1 <- trainy
  rf1 <- ranger::ranger(var1 ~ ., data = data.dev, mtry = mtry (p), num.trees = num.trees,
                        min.node.size = min.node.size, num.threads = num.threads,
                        verbose = verbose)
  rf.pred <- stats::predict(rf1, data = predx, num.trees = num.trees, type = type,
                            num.threads = num.threads)$predictions
  data.dev1 <- longlat
  data.pred1 <- longlatpredx
  data.dev1$var1 <- trainy - stats::predict(rf1, data.dev)$predictions

  sp::coordinates(data.dev1) = ~ LON + LAT
  vgm1 <- gstat::variogram(var1 ~ 1, data.dev1)
  model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args))
  if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to
    variogram. To allow gstat running, the range was set to be positive by
    using min(vgm1$dist). ", "\n"))
  if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
  sp::coordinates(data.pred1) = ~ LON + LAT
  ok.pred <- gstat::krige(var1 ~ 1, data.dev1, data.pred1, model = model.1,
    nmax = nmax, block = block)
  # rfok predictions
  pred.rfok <- rf.pred + ok.pred$var1.pred
  rfok.pred <- cbind(longlatpredx, pred.rfok, ok.pred$var1.var)
  names(rfok.pred) <- c("LON", "LAT", "Predictions", "Variances")
  rfok.pred
}
