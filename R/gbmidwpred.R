#' @title Generate spatial predictions using the hybrid method of generalized boosted regression
#' modeling and inverse distance weighting (gbmidw)
#'
#' @description This function is to make spatial predictions using the hybrid
#' method of generalized boosted regression modeling and inverse distance
#' weighting.
#'
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
#' @param cv.fold integer; number of folds in the cross-validation. it is also
#' the number of cross-validation folds to perform within gbm. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
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
#' @param idp numeric; specify the inverse distance weighting power.
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param n.cores The number of CPU cores to use. See gbm for details. By
#' default, 6 is used.
#' @param ... other arguments passed on to gbm.
#'
#' @return A list with the following components:
#' for numerical data: me, rme, mae, rmae, mse, rmse, rrmse and vecv; or vecv
#' for categorical data: correct classification rate (ccr.cv) and kappa (kappa.cv)
#'
#' @note This function is largely based on gbm.
#'
#' @references Li, J., J. Siwabessy, M. Tran, Z. Huang, and A. Heap. 2013.
#' Predicting Seabed Hardness Using Random Forest in R. Pages 299-329 in Y.
#' Zhao and Y. Cen, editors. Data Mining Applications with R. Elsevier.
#'
#' Li, J. 2013. Predicting the spatial distribution of seabed gravel content
#' using random forest, spatial interpolation methods and their hybrid methods.
#' Pages 394-400  The International Congress on Modelling and Simulation
#' (MODSIM) 2013, Adelaide.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Greg Ridgeway with contributions from others (2015). gbm: Generalized
#' Boosted Regression Models. R package version 2.1.1.
#' https://CRAN.R-project.org/package=gbm
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#' data(petrel.grid)
#' gbmidwpred1 <- gbmidwpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)], petrel[, 3],
#'   petrel.grid[, c(1,2)], petrel.grid, family = "gaussian", n.cores=6,
#'   nmax = 12)
#' names(gbmidwpred1)
#' }
#'
#' @export
gbmidwpred <- function (longlat, trainx, trainy, longlatpredx, predx,
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
  idp = 2,
  nmax = 12,
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
  gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev, set = list(idp = idp), nmax=nmax)
  idw.pred <- stats::predict(gstat1, data.pred)$var1.pred
  gbmidw.pred1 <- pred.gbm1 + idw.pred
  gbmidw.pred <- cbind(longlatpredx, gbmidw.pred1)
  names(gbmidw.pred) <- c("LON", "LAT", "Predictions")
  gbmidw.pred
}

