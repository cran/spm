#' @title Cross validation, n-fold for the average of the hybrid method of random
#' forest in ranger (RG) and ordinary kriging and the hybrid method of RG
#' and inverse distance weighting  (RGOKRGIDW)
#'
#' @description This function is a cross validation function for the average of
#' the hybrid method of random forest in ranger (RG) and ordinary kriging and the hybrid
#' method of RG and inverse distance weighting  (RGOKRGIDW).
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
#' @param num.trees number of trees. By default, 500 is used.
#' @param min.node.size Default 1 for classification, 5 for regression.
#' @param num.threads number of threads. Default is number of CPUs available.
#' @param verbose Show computation status and estimated runtime.Default is FALSE.
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
#' for numerical data: me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv.
#'
#' @note This function is largely based on rfokrfidw.  When 'A zero or negative
#' range was fitted to variogram' occurs, to allow gstat running, the range was
#' set to be positive by using min(vgm1$dist). In this case, caution should be
#' taken in applying this method, although sometimes it can still outperform IDW
#' and OK.
#'
#' @references Li, J. 2013. Predicting the spatial distribution of seabed
#' gravel content using random forest, spatial interpolation methods and their
#' hybrid methods. Pages 394-400  The International Congress on Modelling and
#' Simulation (MODSIM) 2013, Adelaide.
#'
#' Wright, M. N. & Ziegler, A. (2017). ranger: A Fast Implementation
#' of Random Forests for High Dimensional Data in C++ and R. J Stat Softw 77:1-17.
#' http://dx.doi.org/10.18637/jss.v077.i01.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#'
#' rgokrgidwcv1 <- rgokrgidwcv(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 5],
#' predacc = "ALL")
#' rgokrgidwcv1
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' rgokrgidwcv1 <- rgokrgidwcv(petrel[, c(1,2)], petrel[, c(1,2,6:9)], petrel[, 5],
#' predacc = "VEcv")
#' VEcv [i] <- rgokrgidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for RFOKRFIDW", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' rgokrgidwcv1 <- rgokrgidwcv(petrel[, c(1,2)], petrel[, c(1,2,6:9)], petrel[, 5],
#' predacc = "ALL")
#' measures <- rbind(measures, rgokrgidwcv1$vecv)
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for RFOKRFIDW", ylab = "VEcv (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
rgokrgidwcv <- function (longlat, trainx, trainy, cv.fold = 10, mtry = function(p)
  max(1, floor(sqrt(p))),  num.trees = 500, min.node.size = NULL, num.threads = NULL,
  verbose = FALSE, idp = 2, nmaxok = 12, nmaxidw = 12, vgm.args = ("Sph"), block = 0,
  predacc = "VEcv" , ...) {
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
    data.dev <- trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    data.dev$var1 <- trainy[idx != i]
    all.rf1 <- ranger::ranger(var1 ~ ., data = data.dev, mtry = mtry (p), num.trees = num.trees,
                              min.node.size = min.node.size, num.threads = num.threads,
                              verbose = verbose)
    data.dev1 <- longlat[idx != i, , drop = FALSE] # for idw
    data.pred1 <- longlat[idx == i, , drop = FALSE] # for idw
    # rg predictions
    pred.rf1 <- stats::predict(all.rf1, data = data.pred)$predictions
    # residuals of rg
    data.dev1$var1 <- trainy[idx != i] - stats::predict(all.rf1, data.dev)$predictions
    # idw of the residuals
    gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
                             LAT, data = data.dev1, set = list(idp = idp), nmax=nmaxidw)
    pred.idw1 <- stats::predict(gstat1, data.pred1)
    # ok of the residuals
    sp::coordinates(data.dev1) = ~ LON + LAT
    vgm1 <- gstat::variogram(var1 ~ 1, data.dev1)
    model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args))
    if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to
                                    variogram. To allow gstat running, the range was set to be positive by
                                    using min(vgm1$dist). ", "\n"))
    if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
    sp::coordinates(data.pred1) = ~ LON + LAT
    pred.ok1 <- gstat::krige(var1 ~ 1, data.dev1, data.pred1,
                             model = model.1, nmax = nmaxok, block = block)

    # rgokrgidw predictions
    cv.pred[idx == i] <- (pred.idw1$var1.pred + pred.ok1$var1.pred) / 2 + pred.rf1
  }
  # predicitve accuracy assessment
  predictive.accuracy <- NULL
  if (predacc == "VEcv") {predictive.accuracy = vecv(trainy, cv.pred)} else (
    if (predacc == "ALL") {predictive.accuracy = pred.acc(trainy, cv.pred)} else (
      stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
