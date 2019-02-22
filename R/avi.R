#' @title Averaged variable importance based on random forest
#'
#' @description This function is to derive an averaged variable importance based
#' on random forest
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param mtry a function of number of remaining predictor variables to use as
#' the mtry parameter in the randomForest call.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param importance imprtance of predictive variables.
#' @param maxk maxk split value. By default, 4 is used.
#' @param nsim iteration number. By default, 100 is used.
#' @param corr.threshold correlation threshold and the defaults value is 0.5.
#' @param ... other arguments passed on to randomForest.
#'
#' @return A list with the following components: averaged variable importance
#' (avi), column number of importance variable in trainx arranged from the most
#' important to the least important (impvar), names of importance variable
#' arranged from the most important to the least important (impvar2)
#'
#' @references Smith, S.J., Ellis, N., Pitcher, C.R., 2011. Conditional variable
#' importance in R package extendedForest.
#'
#' Li, J. 2013. Predicting the spatial distribution of seabed gravel content
#' using random forest, spatial interpolation methods and their hybrid methods.
#' Pages 394-400  The International Congress on Modelling and Simulation
#' (MODSIM) 2013, Adelaide.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' \dontrun{
#' data(petrel)
#' set.seed(1234)
#' avi1 <- avi(petrel[, c(1,2, 6:9)], petrel[, 5], nsim = 10)
#' avi1
#'
#' avi1 <- avi(petrel[, c(1), drop = FALSE], petrel[, 5], nsim = 10)
#' avi1
#' }
#'
#' @export
avi <- function (trainx, trainy, mtry = if (!is.null(trainy) &&
  !is.factor(trainy)) max(floor(ncol(trainx) / 3), 1) else
  floor(sqrt(ncol(trainx))), ntree = 500, importance=TRUE,
  maxk = c(4), nsim = 100, corr.threshold=0.5, ...) {
  imp <- array(0,dim=c(dim(trainx)[2], length(maxk), nsim))
  for (sim in 1:nsim) {
    for (lev in 1:length(maxk)) {
      RF <- randomForest::randomForest(trainx, trainy, maxLevel = maxk[lev],
      importance = importance, ntree = ntree, mtry = mtry,
      corr.threshold = corr.threshold)
     imp[,lev,sim] <- RF$importance[,1]
    }
  }
  dimnames(imp) <- list(rownames(RF$importance), as.character(maxk), NULL)
  imp <- as.data.frame.table(imp)
  names(imp) <- c("var","maxk","sim","importance")
  avi <- stats::aggregate(. ~ var, data = imp[, c(1,4)], mean)
  impvar <- (1:ncol(trainx))[order(avi$importance, decreasing = TRUE)]
  impvar2 <- names(trainx)[impvar]
  list(avi = avi, impvar = impvar, impvar2 = impvar2)
 }
