#' @title Generate spatial predictions using random forest in ranger (RG)
#'
#' @description This function is to make spatial predictions using random forest
#' in ranger.
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param longlatpredx a dataframe contains longitude and latitude of point
#' locations (i.e., the centres of grids) to be predicted.
#' @param predx a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param mtry Number of variables to possibly split at in each node. Default is the
#' (rounded down) square root of the number variables.
#' @param num.trees number of trees. By default, 500 is used.
#' @param min.node.size Default 1 for classification, 5 for regression.
#' @param type Type of prediction. One of 'response', 'se', 'terminalNodes' with
#' default 'response'. See ranger::predict.ranger for details.
#' @param num.threads number of threads. Default is number of CPUs available.
#' @param verbose Show computation status and estimated runtime.Default is FALSE.
#' @param ... other arguments passed on to randomForest.
#'
#' @return A dataframe of longitude, latitude and predictions.
#'
#' @note This function is largely based on rfpred.
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
#' set.seed(1234)
#' rgpred1 <- rgpred(petrel[, c(1,2, 6:9)], petrel[, 5], petrel.grid[, c(1,2)],
#' petrel.grid, num.trees = 500)
#' names(rgpred1)
#' }
#'
#' @export
rgpred <- function (trainx, trainy, longlatpredx, predx,
  mtry = if (!is.null(trainy) && !is.factor(trainy)) max(floor(ncol(trainx) / 3)
  , 1) else floor(sqrt(ncol(trainx))), num.trees = 500, min.node.size = NULL,
  type = "response", num.threads = NULL, verbose = FALSE, ...) {
  data.dev <- trainx
  data.dev$var1 <- trainy
  rf1 <- ranger::ranger(var1 ~ ., data = data.dev, mtry = mtry, num.trees = num.trees,
                        min.node.size = min.node.size, num.threads = num.threads,
                        verbose = verbose)
  pred.rf1 <- stats::predict(rf1, data = predx, num.trees = num.trees, type = type,
                       num.threads = num.threads)$predictions
  rf.pred <- cbind(longlatpredx, pred.rf1)
  names(rf.pred) <- c("LON", "LAT", "Predictions")
  rf.pred
}
