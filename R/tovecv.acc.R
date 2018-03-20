#' @title Convert error measures to vecv
#'
#' @description tovecv can be used to convert existing predictive error measures to vecv.
#' For the definition of vecv, please see function vecv in library (spm). The error measures considered are
#' mean square error (mse), root mse (rmse), relative rmse (rrmse), standardised rmse (srmse) and
#' mean square reduced error (msre).
#'
#' @param n sample number of validation samples.
#' @param mu mean of validation samples.
#' @param s standard deviation of validation samples.
#' @param m value of an error measure.
#' @param measure a type of error measure (i.e. "mse", "rmse", "rrmse", "srmse" or "msre").
#' @return a numeric number.
#'
#' @references Li, J., 2016. Assessing spatial predictive models in the environmental sciences: accuracy.
#' measures, data variation and variance explained. Environmental Modelling & Software 80 1-8.
#'
#' Li, J., 2017. Assessing the accuracy of predictive models for numerical data: Not r nor r2, why not?
#' Then what? PLOS ONE 12 (8): e0183250.
#'
#' @author Jin Li
#' @examples
#' n <- 300
#' mu <- 15.5
#' sd <- 8.80
#' mse <- 50.43
#' rmse <- sqrt(mse)
#' rrmse <- rmse / mu * 100
#' srmse <- rmse / sd
#' msre <- mse / sd ^ 2
#' tovecv(n=n, mu=mu, s=sd, m=mse, measure="mse")
#'
#' tovecv(n=n, mu=mu, s=sd, m=rmse, measure="rmse")
#'
#' tovecv(n=n, mu=mu, s=sd, m=rrmse, measure="rrmse")
#'
#' tovecv(n=n, mu=mu, s=sd, m=srmse, measure="srmse")
#'
#' tovecv(n=n, mu=mu, s=sd, m=msre, measure="msre")
#'
#' @export
tovecv <- function(n, mu, s, m, measure = c("mse", "rmse", "rrmse", "srmse", "msre")) {
    if(measure == "mse")
      ve <- (1 - (n / ((n-1) * s ^ 2)) * m) * 100
    if(measure == "rmse")
      ve <- (1 - (n / ((n-1) * s ^ 2)) * m ^ 2) * 100
    if(measure == "rrmse")
      ve <- (1 - (n * mu ^ 2 / ((n-1) * s ^ 2 * 100 ^ 2)) * m ^ 2) * 100
    if(measure == "srmse")
      ve <- (1 - (n / (n-1)) * m ^ 2) * 100
    if(measure == "msre")
      ve <- (1 - (n / (n-1)) * m) * 100
    ve
}

