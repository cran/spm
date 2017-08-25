#' @title Predictive error and accuracy measures for predictive models based on cross-validation
#'
#' @description This function is used to calculate the mean error (me), mean absolute error
#' (mae), mean squared error (mse), relative me (rme), relative mae (rmae),
#' root mse (rmse), relative rmse (rrmse), variance explained by predictive
#' models based on cross-validation (vecv). They are based on the differences
#' between the predicted values for and the observed values of validation
#' samples for cross-validation.
#'
#' @param obs a vector of observation values of validation samples.
#' @param pred a vector of prediction values of predictive models for validation samples.
#' @return A list with the following components:
#' me, rme, mae, rmae, mse, rmse, rrmse and vecv.
#' @references Li, J., 2016. Assessing spatial predictive models in the environmental sciences: accuracy
#' measures, data variation and variance explained. Environmental Modelling & Software 80 1-8.
#' @author Jin Li
#' @examples
#' set.seed(1234)
#' x <- sample(1:30, 30)
#' e <- rnorm(30, 1)
#' y <- x + e
#' pred.acc(x, y)
#'
#' y <- 0.8 * x + e
#' pred.acc(x, y)
#'
#' @export
pred.acc <- function (obs, pred) {
  mu <- mean(obs)
  me <- mean(obs - pred)
  mae <- mean(abs(obs - pred))
  mse <-  mean((obs - pred) ^ 2)
  rme <- me / mu * 100
  rmae <- mae / mu * 100
  rmse <- sqrt(mse)
  rrmse <- rmse / mu * 100
  vecv <- (1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2)) * 100
  list(me = me, rme = rme, mae = mae, rmae = rmae, mse = mse, rmse = rmse,
    rrmse = rrmse, vecv = vecv)
}
