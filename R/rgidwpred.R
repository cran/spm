#' @title Generate spatial predictions using the hybrid method of random forest
#' in ranger and inverse distance weighting (RGIDW)
#'
#' @description This function is to make spatial predictions using the hybrid
#' method of random forest in ranger and inverse distance weighting (RGIDW).
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
#' @param num.trees number of trees. By default, 500 is used.
#' @param min.node.size Default 1 for classification, 5 for regression.
#' @param type Type of prediction. One of 'response', 'se', 'terminalNodes' with
#' default 'response'. See ranger::predict.ranger for details.
#' @param num.threads number of threads. Default is number of CPUs available.
#' @param verbose Show computation status and estimated runtime.Default is FALSE.
#' @param idp numeric; specify the inverse distance weighting power.
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param ... other arguments passed on to randomForest or gstat.
#'
#' @return A dataframe of longitude, latitude and predictions.
#'
#' @note This function is largely based on rfidwpred.
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
#' rgidwpred1 <- rgidwpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 3],
#' petrel.grid[, c(1,2)], petrel.grid, num.trees = 500, idp = 2, nmax = 12)
#' names(rgidwpred1)
#' }
#'
#' @export
rgidwpred <- function (longlat, trainx, trainy, longlatpredx, predx, mtry =
  function(p) max(1, floor(sqrt(p))), num.trees = 500, min.node.size = NULL,
  type = "response", num.threads = NULL, verbose = FALSE, idp = 2, nmax = 12, ...) {
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
  gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev1, set = list(idp = idp), nmax=nmax)
  idw.pred <- stats::predict(gstat1, data.pred1)$var1.pred
  rfidw.pred1 <- rf.pred + idw.pred
  rfidw.pred <- cbind(longlatpredx, rfidw.pred1)
  names(rfidw.pred) <- c("LON", "LAT", "Predictions")
  rfidw.pred
}
