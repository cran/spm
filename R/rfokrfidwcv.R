#' @title Cross validation, n-fold for the average of the hybrid method of random forest and ordinary kriging and the hybrid method of random forest and inverse distance weighting  (RFOKRFIDW)
#'
#' @description This function is a cross validation function for the average of
#' the hybrid method of random forest and ordinary kriging and the hybrid
#' method of random forest and inverse distance weighting  (RFOKRFIDW).
#'
#' @param longlat a dataframe contains longitude and latitude of validation
#' samples.
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
#' @param nmaxidw for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used for IDW.
#' @param nmaxok for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used for OK.
#' @param vgm.args arguments for vgm, e.g. variogram model of response
#' variable and anisotropy parameters. see notes vgm in gstat for details.
#' By default, "Sph" is used.
#' @param block block size. see krige in gstat for details.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to randomForest or gstat.
#'
#' @return A list with the following components:
#' for numerical data: me, rme, mae, rmae, mse, rmse, rrmse and vecv; or vecv.
#'
#' @note This function is largely based on rf.cv (see Li et al. 2013) and
#' rfcv in randomForest.  When 'A zero or negative range was fitted to
#' variogram' occurs, to allow gstat running, the range was set to be positive by
#' using min(vgm1$dist). In this case, caution should be taken in applying this
#' method, although sometimes it can still outperform IDW and OK.
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
#' rfokrfidwcv1 <- rfokrfidwcv(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 5],
#' predacc = "ALL")
#' rfokrfidwcv1
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' rfokrfidwcv1 <- rfokrfidwcv(petrel[, c(1,2)], petrel[, c(1,2,6:9)], petrel[, 5],
#' predacc = "VEcv")
#' VEcv [i] <- rfokrfidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for RFOKRFIDW", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' rfokrfidwcv1 <- rfokrfidwcv(petrel[, c(1,2)], petrel[, c(1,2,6:9)], petrel[, 5],
#' predacc = "ALL")
#' measures <- rbind(measures, rfokrfidwcv1$vecv)
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for RFOKRFIDW", ylab = "VEcv (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
rfokrfidwcv <- function (longlat, trainx, trainy, cv.fold = 10, mtry = function(p)
  max(1, floor(sqrt(p))), ntree = 500, idp = 2, nmaxok = 12, nmaxidw = 12,
  vgm.args = ("Sph"), block = 0, predacc = "VEcv", ...) {
  names(longlat) <- c("LON", "LAT")
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (classRF) {
    stop ("This function is not for categorical response variable")
  } else {
    f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4 / 5), Inf))
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  # cross validation
  rf.pred <- NULL
  ok.pred <- NULL
  idw.pred <- NULL
  for (i in 1:cv.fold) {
    all.rf1 <- randomForest::randomForest(trainx[idx != i, , drop = FALSE],
    trainy[idx != i], mtry = mtry(p), ntree=ntree)
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    data.dev1 <- longlat[idx != i, , drop = FALSE] # for OK
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for OK
  # rf predictions
    rf.pred[idx == i] <- stats::predict(all.rf1, data.pred)
  # the residuals of rf
    data.dev1$var1 <- trainy[idx != i] - stats::predict(all.rf1, data.dev)
  # idw of the residuals
    gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev1, set = list(idp = idp), nmax=nmaxidw)
  # idw predictions
    idw.pred[idx == i] <- stats::predict(gstat1, data.pred1)$var1.pred
  # ok of the residuals
    sp::coordinates(data.dev1) = ~ LON + LAT
    vgm1 <- gstat::variogram(var1 ~ 1, data.dev1)
    model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args))
    if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to
      variogram. To allow gstat running, the range was set to be positive by
      using min(vgm1$dist). ", "\n"))
    if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
    sp::coordinates(data.pred1) = ~ LON + LAT
    ok.pred[idx == i] <- gstat::krige(var1 ~ 1, data.dev1, data.pred1,
    model = model.1, nmax = nmaxok, block = block)$var1.pred
  }
  # rfokrfidw predictions
  cv.pred <- rf.pred + (ok.pred + idw.pred) / 2
  # predicitve accuracy assessment
  predictive.accuracy <- NULL
  if (predacc == "VEcv") {predictive.accuracy = vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
