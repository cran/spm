\name{sponge.grid}
\alias{sponge.grid}
\docType{data}
\title{A dataset of predictors for generating sponge species richness in a selected region in the Timor Sea region, northern Australia marine margin}
\description{This dataset contains 95530 rows of 7 predictive variables including longitude (easting), latitude (northing), topographic position index (tpi3), variance of backscatter (var7), entropy (entro7), backscatter at incidence angle 11 degree (bs11), and backscatter at incidence angle 34 degree (bs34).}
\usage{data("sponge.grid")}
\format{
  A data frame with 95530 rows on the following 7 variables.
  \describe{
    \item{\code{easting}}{a numeric vector, m}
    \item{\code{northing}}{a numeric vector, m}
    \item{\code{tpi3}}{a numeric vector, no unit}
    \item{\code{var7}}{a numeric vector, dB^2}
    \item{\code{entro7}}{a numeric vector, no unit}
    \item{\code{bs11}}{a numeric vector, dB}
    \item{\code{bs34}}{a numeric vector, dB}
  }
}
\details{For details, please see the source. This dataset was used to produce the figure of predictions in the paper listed in the source.}
\source{Li, J., B. Alvarez, J. Siwabessy, M. Tran, Z. Huang, R. Przeslawski, L. Radke, F. Howard, and S. Nichol. 2017. Application of random forest, generalised linear model and their hybrid methods with geostatistical techniques to count data: Predicting sponge species richness. Environmental Modelling & Software,97: 112-129.}
\keyword{datasets}
