#' @title Predictive error and accuracy measures for predictive models based on cross-validation
#'
#' @description This function is used to calculate the mean error (me), mean absolute error
#' (mae), mean squared error (mse), relative me (rme), relative mae (rmae),
#' root mse (rmse), relative rmse (rrmse), variance explained by predictive
#' models based on cross-validation (vecv), and Legates and McCabe's E1 (e1) for numerical data; and
#' it also calculates correct classification rate (ccr), kappa (kappa), sensitivity (sens), specificity
#' (spec), and true skill statistic (tss) for categorical data with the observed (obs) data specified
#' as factor. They are based on the differences between the predicted values for and the observed values
#' of validation samples for cross-validation. For 0 and 1 data, the observed values need to be specified
#' as factor in order to use accuracy measures for categorical data. Moreover, sens, spec, tss and rmse are
#' for categorical data with two levels (e.g. presence and absence data).
#'
#' @param obs a vector of observation values of validation samples.
#' @param pred a vector of prediction values of predictive models for validation samples.
#' @return A list with the following components:
#' me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1 for numerical data;
#' ccr, kappa, sens, spec and tss for categorical data with two levels; and
#' ccr, kappa for categorical data with more than two levels.
#'
#' @references Li, J., 2016. Assessing spatial predictive models in the environmental sciences: accuracy
#' measures, data variation and variance explained. Environmental Modelling & Software 80 1-8.
#'
#' Li, J., 2017. Assessing the accuracy of predictive models for numerical data: Not r nor r2, why not?
#' Then what? PLOS ONE 12 (8): e0183250.
#'
#' Allouche, O., Tsoar, A., Kadmon, R., 2006. Assessing the accuracy of species distribution models:
#' prevalence, kappa and true skill statistic (TSS). Journal of Applied Ecology 43 1223-1232.
#'
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
  if (is.factor(obs)) {
    if (length(levels(obs)) == 2) {
      data1 <- as.data.frame(cbind(pred, obs))
      kappa <- psy::ckappa(data1)$kappa
      ccr <- sum(diag(table(data1))) / sum(table(data1)) * 100
      sensitivity <- table(data1)[1, 1] / sum(table(data1)[,1])
      specificity <- table(data1)[2, 2] / sum(table(data1)[,2])
      tss <- sensitivity + specificity - 1
      list (kappa = kappa, ccr = ccr, sens = sensitivity, spec = specificity, tss = tss)
      } else {
        data1 <- as.data.frame(cbind(pred, obs))
        kappa <- psy::ckappa(data1)$kappa
        ccr <- sum(diag(table(data1))) / sum(table(data1)) * 100
        list (kappa = kappa, ccr = ccr)
      }
    } else {
      mu <- mean(obs)
      me <- mean(obs - pred)
      mae <- mean(abs(obs - pred))
      mse <-  mean((obs - pred) ^ 2)
      rme <- me / mu * 100
      rmae <- mae / mu * 100
      rmse <- sqrt(mse)
      rrmse <- rmse / mu * 100
      vecv <- (1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2)) * 100
      e1 <- (1 - sum(abs(obs - pred)) / sum(abs(obs - mean(obs)))) * 100
      list(me = me, rme = rme, mae = mae, rmae = rmae, mse = mse, rmse = rmse,
           rrmse = rrmse, vecv = vecv, e1 = e1)
    }
}
