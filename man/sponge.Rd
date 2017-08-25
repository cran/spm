\name{sponge}
\alias{sponge}
\docType{data}
\title{A dataset of sponge species richness in the Timor Sea region, northern Australia marine margin}
\description{This dataset contains 77 samples of 8 variables including longitude (easting), latitude (northing), sponge, topographic position index (tpi3), variance of backscatter (var7), entropy (entro7), backscatter at incidence angle 11 degree (bs11), and backscatter at incidence angle 34 degree (bs34).}
\usage{data("sponge")}
\format{
  A data frame with 77 observations on the following 8 variables.
  \describe{
    \item{\code{easting}}{a numeric vector, m}
    \item{\code{northing}}{a numeric vector, m}
    \item{\code{sponge}}{a numeric vector, no unit}
    \item{\code{tpi3}}{a numeric vector, no unit}
    \item{\code{var7}}{a numeric vector, dB^2}
    \item{\code{entro7}}{a numeric vector, no unit}
    \item{\code{bs11}}{a numeric vector, dB}
    \item{\code{bs34}}{a numeric vector, dB}
  }
}
\details{For details, please see the source. This dataset was published as an appendix of the paper listed in the source. Where the long and lat were replaced with easting and northing for prediction purpose.}
\source{Li, J., B. Alvarez, J. Siwabessy, M. Tran, Z. Huang, R. Przeslawski, L. Radke, F. Howard, and S. Nichol. 2017. Application of random forest, generalised linear model and their hybrid methods with geostatistical techniques to count data: Predicting sponge species richness. Environmental Modelling & Software, 97: 112-129}
\keyword{datasets}
