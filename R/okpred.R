#' @title Generate spatial predictions using ordinary kriging (OK)
#'
#' @description This function is to make spatial predictions using ordinary
#' kriging.
#'
#' @param longlat a dataframe contains longitude and latitude of point
#' samples.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in longlat.
#' @param longlat2 a dataframe contains longitude and latitude of point locations
#' (i.e., the centres of grids) to be predicted.
#' @param nmax for local kriging: the number of nearest observations that
#' should be used for a kriging prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param transformation transform the response variable to normalise the data;
#' can be "sqrt" for square root, "arcsine" for arcsine, "log" or "none"
#' for non transformation. By default, "none" is used.
#' @param delta numeric; to avoid log(0) in the log transformation.
#' @param vgm.args arguments for vgm, e.g. variogram model of response
#' variable and anisotropy parameters. see notes vgm in gstat for details.
#' By default, "Sph" is used.
#' @param anis anisotropy parameters: see notes vgm in gstat for details.
#' @param alpha direction in plane (x,y). see variogram in gstat for details.
#' @param block block size. see krige in gstat for details.
#' @param ... other arguments passed on to gstat.
#'
#' @return A dataframe of longitude, latitude, predictions and variances.
#'
#' @references Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
#' package. Computers & Geosciences, 30: 683-691.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' library(sp)
#' data(swmud)
#' data(sw)
#' okpred1 <- okpred(swmud[, c(1,2)], swmud[, 3], sw, nmax = 7, transformation =
#' "arcsine", vgm.args = ("Sph"))
#' names(okpred1)
#' }
#'
#' @export
okpred <- function (longlat, trainy, longlat2, nmax = 12, transformation =
  "none", delta = 1, vgm.args = ("Sph"), anis = c(0, 1), alpha = 0, block = 0,  ...) {
  if (transformation == "none") {trainy1 = trainy} else (
  if (transformation == "sqrt") {trainy1 = sqrt(trainy)} else (
  if (transformation == "arcsine") {trainy1 = asin(sqrt(trainy / 100))} else (
  if (transformation == "log") {trainy1 = log(trainy + delta)} else (
  stop ("This transfromation is not supported in this version!")))))
  names(longlat) <- c("LON", "LAT")
  names(longlat2) <- c("LON", "LAT")
  data.dev <- longlat
  data.pred <- longlat2
  data.dev$var1 <- trainy1
  sp::coordinates(data.dev) = ~ LON + LAT
  vgm1 <- gstat::variogram(var1 ~ 1, data.dev, alpha = alpha)
  model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args, anis = anis))
  sp::coordinates(data.pred) = ~LON + LAT
  ok.pred1 <- gstat::krige(var1 ~ 1, data.dev, data.pred,
    model = model.1, nmax = nmax, block = block)
  if (transformation == "none") {ok.pred2 = ok.pred1$var1.pred}
  if (transformation == "sqrt") {ok.pred2 = ok.pred1$var1.pred ^ 2}
  if (transformation == "arcsine") {ok.pred2 = (sin(ok.pred1$var1.pred)) ^ 2 * 100}
  if (transformation == "log") {ok.pred2 = exp(ok.pred1$var1.pred)-delta}
  ok.pred <- cbind(longlat2, ok.pred2, ok.pred1$var1.var)
  ok.pred
}
