#' @title Cross validation, n-fold for the hybrid method of random forest and
#' inverse distance weighting (RFIDW)
#'
#' @description This function is a cross validation function for the hybrid
#' method of random forest and inverse distance weighting (RFIDW).
#'
#' @param longlat a dataframe contains longitude and latitude of point
#' samples (i.e., trainx and trainy).
#' @param trainx a dataframe or matrix contains columns of predictive variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
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
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to randomForest or gstat.
#'
#' @return A list with the following components:
#' for numerical data: me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv.
#'
#' @note This function is largely based on rf.cv (see Li et al. 2013) and
#' rfcv in randomForest.
#'
#' @references Li, J. 2013. Predicting the spatial distribution of seabed
#' gravel content using random forest, spatial interpolation methods and their
#' hybrid methods. Pages 394-400  The International Congress on Modelling and
#' Simulation (MODSIM) 2013, Adelaide.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#'
#' rfidwcv1 <- rfidwcv(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 5],
#' predacc = "ALL")
#' rfidwcv1
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' rfidwcv1 <- rfidwcv(petrel[, c(1,2)], petrel[, c(1,2,6:9)], petrel[, 5],
#' predacc = "VEcv")
#' VEcv [i] <- rfidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for RFIDW", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' rfidwcv1 <- rfidwcv(petrel[, c(1,2)], petrel[, c(1,2,6:9)], petrel[, 5],
#' predacc = "ALL")
#' measures <- rbind(measures, rfidwcv1$vecv)
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for RFIDW", ylab = "VEcv (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
rfidwcv <- function (longlat, trainx, trainy, cv.fold = 10, mtry = function(p)
  max(1, floor(sqrt(p))),  ntree = 500, idp = 2, nmax = 12, predacc = "VEcv"
  , ...) {
  names(longlat) <- c("LON", "LAT")
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (classRF) {
    stop ("This function is not for categorical response variable")
  }     else {
    f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4 / 5), Inf))
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  # cross validation
  cv.pred <- NULL
  for (i in 1:cv.fold) {
    all.rf1 <- randomForest::randomForest(trainx[idx != i, , drop = FALSE],
    trainy[idx != i], mtry = mtry(p), ntree=ntree)
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw
  # rf predictions
    pred.rf1 <- stats::predict(all.rf1, data.pred)
  # residuals of rf
    data.dev1$var1 <- trainy[idx != i] - stats::predict(all.rf1, data.dev)
  # idw of the residuals
    gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev1, set = list(idp = idp), nmax=nmax)
    pred.idw1 <- stats::predict(gstat1, data.pred1)
  # rfidw predictions
    cv.pred[idx == i] <- pred.idw1$var1.pred + pred.rf1
  }
  # predicitve accuracy assessment
  predictive.accuracy <- NULL
  if (predacc == "VEcv") {predictive.accuracy = vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
