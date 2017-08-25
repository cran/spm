#' @title Generate spatial predictions using random forest (RF)
#'
#' @description This function is to make spatial predictions using random forest.
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param longlatpredx a dataframe contains longitude and latitude of point
#' locations (i.e., the centres of grids) to be predicted.
#' @param predx a dataframe or matrix contains columns of predictive variables
#' for the grids to be predicted.
#' @param mtry a function of number of remaining predictor variables to use as
#' the mtry parameter in the randomForest call.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param ... other arguments passed on to randomForest.
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
#' rfpred1 <- rfpred(petrel[, c(1,2, 6:9)], petrel[, 5], petrel.grid[, c(1,2)],
#' petrel.grid, ntree = 500)
#' names(rfpred1)
#' }
#'
#' @export
rfpred <- function (trainx, trainy, longlatpredx, predx,
  mtry = if (!is.null(trainy) && !is.factor(trainy)) max(floor(ncol(trainx) / 3)
  , 1) else floor(sqrt(ncol(trainx))), ntree = 500,  ...) {
  rf.1 <- randomForest::randomForest(trainx, trainy, mtry=mtry, ntree=ntree)
  pred.rf1 <- stats::predict(rf.1, predx)
  rf.pred <- cbind(longlatpredx, pred.rf1)
  names(rf.pred) <- c("LON", "LAT", "Predictions")
  rf.pred
}
