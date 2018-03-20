#' @title Cross validation, n-fold for random forest in ranger (RG)
#'
#' @description This function is a cross validation function for random forest in ranger.
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param mtry Number of variables to possibly split at in each node. Default is the
#' (rounded down) square root of the number variables.
#' @param num.trees number of trees. By default, 500 is used.
#' @param min.node.size Default 1 for classification, 5 for regression.
#' @param num.threads number of threads. Default is number of CPUs available.
#' @param verbose Show computation status and estimated runtime.Default is FALSE.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to randomForest.
#'
#' @return A list with the following components:
#' for numerical data: me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv.
#' for categorical data: correct classification rate (ccr), kappa (kappa), sensitivity (sens),
#' specificity (spec) and true skill statistic (tss)
#'
#' @note This function is largely based on RFcv.
#'
#' @references Li, J. 2013. Predicting the spatial distribution of seabed gravel content
#' using random forest, spatial interpolation methods and their hybrid methods.
#' Pages 394-400  The International Congress on Modelling and Simulation
#' (MODSIM) 2013, Adelaide.
#'
#' Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation
#' of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17.
#' http://dx.doi.org/10.18637/jss.v077.i01.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(hard)
#' data(petrel)
#'
#' rgcv1 <- rgcv(petrel[, c(1,2, 6:9)], petrel[, 5], predacc = "ALL")
#' rgcv1
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' rgcv1 <- rgcv(petrel[, c(1,2,6:9)], petrel[, 5], predacc = "VEcv")
#' VEcv [i] <- rgcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for RF", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' rgcv1 <- rgcv(hard[, c(4:6)], hard[, 17])
#' measures <- rbind(measures, rgcv1$ccr) # for kappa, replace ccr.cv with kappa.cv
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for RF", ylab = "Correct
#' classification rate  (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
rgcv <- function (trainx, trainy, cv.fold = 10, mtry = if (!is.null(trainy) &&
  !is.factor(trainy)) max(floor(ncol(trainx) / 3), 1) else floor(sqrt(ncol(trainx))),
  num.trees = 500, min.node.size = NULL, num.threads = NULL, verbose = FALSE,
  predacc = "ALL", ...) {
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  if (classRF) {
    f <- trainy
  }     else {
    f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4/5), Inf))
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  cv.pred <- NULL
  for (i in 1:cv.fold) {
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    data.dev$var1 <- trainy[idx != i]
    all.rf <- ranger::ranger(var1 ~ ., data = data.dev, mtry = mtry, num.trees = num.trees,
                             min.node.size = min.node.size, num.threads = num.threads,
                             verbose = verbose)
    cv.pred[idx == i] <- stats::predict(all.rf, data = data.pred)$predictions
  }
  predictive.accuracy <- NULL
  if (predacc == "VEcv") {predictive.accuracy = vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
