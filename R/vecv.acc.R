#' @title Variance explained by predictive models based on cross-validation
#'
#' @description vecv is used to calculate the variance explained by predictive
#' models based on cross-validation. The vecv is based on the differences between
#' the predicted values for, and the observed values of, validation samples
#' for cross-validation. It measures the proportion of variation in the
#' validation data explained by the predicted values obtained from predictive
#' models based on cross-validation.
#'
#' @param obs observation values of validation samples.
#' @param pred prediction values of predictive models for validation samples.
#' @return a numeric number.
#'
#' @references Li, J., 2016. Assessing spatial predictive models in the environmental sciences: accuracy.
#' measures, data variation and variance explained. Environmental Modelling & Software 80 1-8.
#'
#' @author Jin Li
#' @examples
#' set.seed(1234)
#' x <- sample(1:30, 30)
#' e <- rnorm(30, 1)
#' y <- x + e
#' vecv(x, y)
#'
#' y <- 0.8 * x + e
#' vecv(x, y)
#'
#' @export
vecv <- function (obs, pred) {
    (1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2)) * 100
}
