\name{hard}
\alias{hard}
\docType{data}
\title{A dataset of seabed hardness in the eastern Joseph Bonaparte Golf, northern Australia marine margin}
\description{This dataset contains 137 samples of 17 variables including area surveyed (Area), easting, northing, prock, bathymetry (bathy), backscatter (bs), local Moran I (bathy.moran), plannar curvature (planar.curv), profile curvature (profile.curv), topographic relief (relief), slope (slope), surface area (surface), topographic position index (tpi), homogeneity of backscatter (homogeneity),  local Moran I of backscatter (bs.moran), variance of backscatter (variance) and seabed hardness (hardness).}
\usage{data("hard")}
\format{
  A data frame with 137 observations on the following 17 variables.
  \describe{
    \item{\code{Area}}{a catergorical vector, no unit}
    \item{\code{easting}}{a numeric vector, m}
    \item{\code{northing}}{a numeric vector, m}
    \item{\code{prock}}{a numeric vector, no unit}
    \item{\code{bathy}}{a numeric vector, meter}
    \item{\code{bs}}{a numeric vector, dB}
    \item{\code{bathy.moran}}{a numeric vector, no unit}
    \item{\code{planar.curv}}{a numeric vector, no unit}
    \item{\code{profile.curv}}{a numeric vector, no unit}
    \item{\code{relief}}{a numeric vector, meter}
    \item{\code{slope}}{a numeric vector, no unit}
    \item{\code{surface}}{a numeric vector, no unit}
    \item{\code{tpi}}{a numeric vector, no unit}
    \item{\code{homogeneity}}{a numeric vector, no unit}
    \item{\code{bs.moran}}{a numeric vector, no unit}
    \item{\code{variance}}{a numeric vector, dB^2}
    \item{\code{hardness}}{a catergorical vector, no unit}
  }
}
\details{For details, please see the source. This dataset was modified by removing 3 samples with missing values from Appendix AA of the book chapter listed in the source.}
\source{Li, J., J. Siwabessy, M. Tran, Z. Huang, and A. Heap. 2013. Predicting Seabed Hardness Using Random Forest in R. Pages 299-329 in Y. Zhao and Y. Cen, editors. Data Mining Applications with R. Elsevier.}
\keyword{datasets}
