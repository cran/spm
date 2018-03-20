#' @title Cross validation, n-fold for inverse distance weighting (IDW)
#'
#' @description This function is a cross validation function for inverse
#' distance weighting.
#'
#' @param longlat a dataframe contains longitude and latitude of validation
#' samples.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in longlat.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param nmax for a local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param idp numeric; specify the inverse distance weighting power.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to gstat.
#'
#' @return A list with the following components:
#'  me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only.
#' @note This function is largely based on rfcv in randomForest and
#' some functions in library(gstat).
#'
#' @references Li, J., 2013. Predictive Modelling Using Random Forest and Its
#' Hybrid Methods with Geostatistical Techniques in Marine Environmental
#' Geosciences, In: Christen, P., Kennedy, P., Liu, L., Ong, K.-L., Stranieri,
#' A., Zhao, Y. (Eds.), The proceedings of the Eleventh Australasian Data
#' Mining Conference (AusDM 2013), Canberra, Australia, 13-15 November 2013.
#' Conferences in Research and Practice in Information Technology, Vol. 146.
#'
#' A. Liaw and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
#' package. Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' library(sp)
#' data(swmud)
#' data(petrel)
#'
#' idwcv1 <- idwcv(swmud[, c(1,2)], swmud[, 3], nmax = 12, idp = 2)
#' idwcv1
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' idwcv1 <- idwcv(petrel[, c(1,2)], petrel[, 3], nmax = 12, predacc = "VEcv")
#' VEcv [i] <- idwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for IDW", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd=2)
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' idwcv1 <- idwcv(swmud[, c(1,2)], swmud[, 3], predacc = "ALL")
#' measures <- rbind(measures, idwcv1$vecv)
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for IDW", ylab="VEcv (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
idwcv <- function (longlat, trainy, cv.fold = 10, nmax = 12, idp = 2, predacc =
  "VEcv", ...) {
  classRF <- is.factor(trainy)
  n <- nrow(longlat)
  p <- ncol(longlat)
  if (classRF) {
    stop ("This function is not for categorical response variable")
  } else {
    f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4/5), Inf))
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,
    length = nlvl[i]))
  }
  names(longlat) <- c("LON", "LAT")
  # cross validation
  cv.pred <- NULL
  for (i in 1 : cv.fold) {
    data.dev <- longlat[idx != i, , drop = FALSE]
    data.pred <- longlat[idx == i, , drop = FALSE]
    data.dev$var1 <- trainy[idx != i]
    # ok
    gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev, set = list(idp = idp), nmax=nmax)
    cv.pred[idx == i] <- stats::predict(gstat1, data.pred)$var1.pred
  }
  if (predacc == "VEcv") {predictive.accuracy = vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
