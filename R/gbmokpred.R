#' @title Generate spatial predictions using the hybrid method of generalized boosted
#' regression modeling and ordinary kriging (gbmok)
#'
#' @description This function is to make spatial predictions using the hybrid
#' method of generalized boosted regression modeling and ordinary kriging.
#' @param longlat a dataframe contains longitude and latitude of point
#' samples (i.e., trainx and trainy).
#' @param trainx a dataframe or matrix contains columns of predictive variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param longlatpredx a dataframe contains longitude and latitude of point locations
#' (i.e., the centres of grids) to be predicted.
#' @param predx a dataframe or matrix contains columns of predictive variables for
#' the grids to be predicted.
#' @param var.monotone an optional vector, the same length as the number of
#' predictors, indicating which variables have a monotone increasing (+1),
#' decreasing (-1), or arbitrary (0) relationship with the outcome. By default,
#' a vector of 0 is used.
#' @param family either a character string specifying the name of the distribution to
#' use or a list with a component name specifying the distribution and any
#' additional parameters needed. See gbm for details. By default, "gaussian" is
#' used.
#' @param n.trees the total number of trees to fit. This is equivalent to the
#' number of iterations and the number of basis functions in the additive
#' expansion. By default, 3000 is used.
#' @param learning.rate a shrinkage parameter applied to each tree in the
#' expansion. Also known as step-size reduction.
#' @param interaction.depth the maximum depth of variable interactions.
#' 1 implies an additive model, 2 implies a model with up to 2-way
#' interactions, etc. By default, 2 is used.
#' @param bag.fraction the fraction of the training set observations randomly
#' selected to propose the next tree in the expansion. By default, 0.5 is used.
#' @param train.fraction The first train.fraction * nrows(data) observations
#' are used to fit the gbm and the remainder are used for computing
#' out-of-sample estimates of the loss function.
#' @param n.minobsinnode minimum number of observations in the trees terminal
#' nodes. Note that this is the actual number of observations not the total
#' weight. By default, 10 is used.
#' @param cv.fold integer; number of cross-validation folds to perform within
#' gbm.
#' @param weights an optional vector of weights to be used in the fitting
#' process. Must be positive but do not need to be normalized.
#' If keep.data = FALSE in the initial call to gbm then it is the user's
#' responsibility to resupply the weights to gbm.more. By default, a vector of
#' 1 is used.
#' @param keep.data a logical variable indicating whether to keep the data and
#' an index of the data stored with the object. Keeping the data and index
#' makes subsequent calls to gbm.more faster at the cost of storing an extra
#' copy of the dataset. By default, 'FALSE' is used.
#' @param verbose If TRUE, gbm will print out progress and performance
#' indicators. By default, 'TRUE' is used.
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param vgm.args arguments for vgm, e.g. variogram model of response
#' variable and anisotropy parameters. see notes vgm in gstat for details.
#' By default, "Sph" is used.
#' @param block block size. see krige in gstat for details.
#' @param n.cores The number of CPU cores to use. See gbm for details. By
#' default, 6 is used.
#' @param ... other arguments passed on to gbm.
#'
#' @return A dataframe of longitude, latitude, predictions and variances. The
#' variances are produced by OK based on the residuals of gbm.
#'
#' @note This function is largely based on gbm. When 'A zero or negative range
#' was fitted to variogram' occurs, to allow OK running, the range was set
#' to be positive by using min(vgm1$dist). In this case, caution should be
#' taken in applying this method, although sometimes it can still outperform
#' IDW and OK.
#'
#' @references Greg Ridgeway with contributions from others (2015). gbm:
#' Generalized Boosted Regression Models. R package version 2.1.1.
#' https://CRAN.R-project.org/package=gbm
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#' data(petrel.grid)
#' gbmokpred1 <- gbmokpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 3],
#'   petrel.grid[, c(1,2)], petrel.grid, family = "gaussian", n.cores=6,
#'   nmax = 12, vgm.args = ("Sph"))
#' names(gbmokpred1)
#' }
#'
#' @export
gbmokpred <- function (longlat, trainx, trainy, longlatpredx, predx,
  var.monotone = rep(0, ncol(trainx)),
  family = "gaussian",
  n.trees = 3000,          # default number of trees
  learning.rate = 0.001,
  interaction.depth = 2,
  bag.fraction = 0.5,
  train.fraction = 1.0,
  n.minobsinnode = 10,
  cv.fold = 10, # becuase of the way used to resample data, we can not do leave-one-out cv.
  weights = rep(1, nrow(trainx)),   # by default set equal
  keep.data = FALSE,
  verbose = TRUE,
  nmax = 12,
  vgm.args = ("Sph"),
  block = 0,
  n.cores = 6,  ...) {
  names(longlat) <- c("LON", "LAT")
  names(longlatpredx) <- c("LON", "LAT")
  gbm1 <- gbm::gbm(trainy ~ ., data=trainx,
    var.monotone = var.monotone,
    distribution = as.character(family),
    n.trees = n.trees,
    shrinkage = learning.rate,
    interaction.depth = interaction.depth,
    bag.fraction = bag.fraction,
    train.fraction = train.fraction,
    n.minobsinnode = n.minobsinnode,
    weights = weights,
    cv.folds = cv.fold,
    n.cores = n.cores)
  best.iter <- gbm::gbm.perf(gbm1, method = "cv")
  print(best.iter)
  pred.gbm1 <- gbm::predict.gbm(gbm1, predx, n.trees =
  best.iter, type = "response")
  data.dev <- longlat
  data.pred <- longlatpredx
  data.dev$var1 <- trainy - gbm::predict.gbm(gbm1, trainx, n.trees = best.iter,
   type = "response")
  sp::coordinates(data.dev) = ~ LON + LAT
  vgm1 <- gstat::variogram(var1 ~ 1, data.dev)
  model.1 <- gstat::fit.variogram(vgm1, gstat::vgm(vgm.args))
  if (model.1$range[2] <= 0) (cat("A zero or negative range was fitted to
    variogram. To allow gstat running, the range was set to be positive by
    using min(vgm1$dist). ", "\n"))
  if (model.1$range[2] <= 0) (model.1$range[2] <- min(vgm1$dist)) # set negative range to be positive
  sp::coordinates(data.pred) = ~ LON + LAT
  ok.pred <- gstat::krige(var1 ~ 1, data.dev, data.pred, model = model.1,
    nmax = nmax, block = block)
  # gbmok predictions
  pred.gbmok <- pred.gbm1 + ok.pred$var1.pred
  gbmok.pred <- cbind(longlatpredx, pred.gbmok, ok.pred$var1.var)
  names(gbmok.pred) <- c("LON", "LAT", "Predictions", "Variances")
  gbmok.pred
}
