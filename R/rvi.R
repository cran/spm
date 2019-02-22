#' @title Relative variable influence based on generalized boosted regression
#' modeling (gbm)
#'
#' @description This function is to to derive a relative variable influence based
#' on generalized boosted regression modeling.
#'
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
#' @param cv.fold integer; number of cross-validation folds to perform within
#' gbm. if > 1, then apply n-fold cross validation; the default is 10, i.e.,
#' 10-fold cross validation that is recommended.
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
#' @param n.cores The number of CPU cores to use. See gbm for details. By
#' default, 6 is used.
#' @param ... other arguments passed on to gbm.
#'
#' @return A list of column number of importance variable in trainx arranged from the most
#' influential to the least influential (impvar), and a dataframe of variables (var), and
#' relative influence (rel.inf)
#'
#' @note This function is largely based on gbm.
#'
#' @references Greg Ridgeway with contributions from others (2015). gbm:
#' Generalized Boosted Regression Models. R package version 2.1.1.
#' https://CRAN.R-project.org/package=gbm
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(sponge)
#' set.seed(1234)
#' rvi1 <- rvi(sponge[, -c(3)], sponge[, 3], family = "poisson", n.cores=2)
#' names(ri1)
#' impvar <- (1:ncol(sponge[, -c(3)]))[ri1$var]
#' }
#'
#' @export
rvi <- function (trainx, trainy,
  var.monotone = rep(0, ncol(trainx)),
  family = "gaussian",
  n.trees = 3000,          # default number of trees
  learning.rate = 0.001,
  interaction.depth = 2,
  bag.fraction = 0.5,
  train.fraction = 1.0,
  n.minobsinnode = 10,
  cv.fold = 10,
  weights = rep(1, nrow(trainx)),   # by default set equal
  keep.data = FALSE,
  verbose = TRUE,
  n.cores = 6, ...) {
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
    keep.data = keep.data,
    cv.folds = cv.fold,
    n.cores = n.cores)
  best.iter <- gbm::gbm.perf(gbm1, method = "cv")
  print(best.iter)
  gbm.rvi <- summary(gbm1, n.trees=best.iter, las=2)
  impvar <- (1:ncol(trainx))[gbm.rvi$var]
  list(impvar = impvar, gbm.rvi = gbm.rvi)
}
