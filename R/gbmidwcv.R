#' @title Cross validation, n-fold for the hybrid method of generalized boosted regression modeling and inverse distance weighting (gbmidw)
#'
#' @description This function is a cross validation function for the hybrid
#' method of generalized boosted regression modeling and inverse distance
#' weighting.
#'
#' @param longlat a dataframe contains longitude and latitude of validation
#' samples.
#' @param trainx a dataframe or matrix contains columns of predictive variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
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
#' If keep.data=FALSE in the initial call to gbm then it is the user's
#' responsibility to resupply the weights to gbm.more. By default, a vector of
#' 1 is used.
#' @param keep.data a logical variable indicating whether to keep the data and
#' an index of the data stored with the object. Keeping the data and index
#' makes subsequent calls to gbm.more faster at the cost of storing an extra
#' copy of the dataset. By default, 'FALSE' is used.
#' @param verbose If TRUE, gbm will print out progress and performance
#' indicators. By default, 'TRUE' is used.
#' @param n.cores The number of CPU cores to use. See gbm for details. By
#' default, 6 is used.
#' @param idp numeric; specify the inverse distance weighting power.
#' @param nmax for local predicting: the number of nearest observations that
#' should be used for a prediction or simulation, where nearest is
#' defined in terms of the space of the spatial locations. By default, 12
#' observations are used.
#' @param predacc can be either "VEcv" for vecv or "ALL" for all measures
#' in function pred.acc.
#' @param ... other arguments passed on to gbm.
#'
#' @return A list with the following components:
#' for numerical data: me, rme, mae, rmae, mse, rmse, rrmse and vecv; or vecv
#' for categorical data: correct classification rate (ccr.cv) and kappa (kappa.cv)
#'
#' @note this function is largely based on rf.cv (see Li et al. 2013),
#' rfcv{randomForest} and gbm
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
#' data(sponge)
#'
#' gbmidwcv1 <- gbmidwcv(sponge[, c(1,2)], sponge[, -c(3)], sponge[, 3],
#' cv.fold = 10, family = "poisson", n.cores=2, predacc = "ALL")
#' gbmidwcv1
#'
#' n <- 20 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' gbmidwcv1 <- gbmidwcv(sponge[, c(1,2)], sponge[, -c(3)], sponge[, 3],
#' cv.fold = 10, family = "poisson", n.cores=2, predacc = "VEcv")
#' VEcv [i] <- gbmidwcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for gbmidw", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#' }
#'
#' @export
gbmidwcv <- function (longlat, trainx, trainy, var.monotone = rep(0, ncol(trainx)),
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
  predacc = "VEcv",
  n.cores = 6, ...) {
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  cv.pred <- NULL
  if (classRF) {
  stop ("This function is not for categorical response variable")
  }     else {
  f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4/5), Inf))
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  names(longlat) <- c("LON", "LAT")
  # cross validation
  gbm.pred <- NULL
  idw.pred <- NULL
  for (i in 1:cv.fold) {
    all.gbm1 <- gbm::gbm(trainy[idx != i] ~ ., data=trainx[idx != i, , drop = FALSE],
    var.monotone = var.monotone,
    distribution = as.character(family),
    n.trees = n.trees,
    shrinkage = learning.rate,
    interaction.depth = interaction.depth,
    bag.fraction = bag.fraction,
    train.fraction = train.fraction,
    n.minobsinnode = n.minobsinnode,
    weights = weights[idx != i],
    cv.folds = cv.fold,
    keep.data = keep.data,
    verbose = verbose,
    n.cores = n.cores)
    # gbm predictions
    data.dev<-trainx[idx != i, , drop = FALSE]
    data.pred <- trainx[idx == i, , drop = FALSE]
    data.dev1<-longlat[idx != i, , drop = FALSE] # for idw
    data.pred1<-longlat[idx == i, , drop = FALSE] # for idw
    best.iter <- gbm::gbm.perf(all.gbm1, method = "cv")
    print(best.iter)
    gbm.pred[idx == i] <- gbm::predict.gbm(all.gbm1, data.pred, n.trees = best.iter, type
    = "response")
    # residuals of gbm
    data.dev1$var1 <- trainy[idx != i] - gbm::predict.gbm(all.gbm1, data.dev, n.trees =
    best.iter, type = "response")
    # idw of the residuals
    gstat1 <- gstat::gstat(id = "var1", formula = var1 ~ 1, locations = ~ LON +
    LAT, data = data.dev1, set = list(idp = idp), nmax=nmax)
    idw.pred[idx == i] <- stats::predict(gstat1, data.pred1)$var1.pred
  }
  # gbmidw predictions
    cv.pred <- gbm.pred + idw.pred
# predicitve accuracy assessment
  predictive.accuracy <- NULL
  if (predacc == "VEcv") {predictive.accuracy = vecv(trainy, cv.pred)} else (
  if (predacc == "ALL") {predictive.accuracy = pred.acc(trainy, cv.pred)} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
