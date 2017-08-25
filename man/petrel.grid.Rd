\name{petrel.grid}
\alias{petrel.grid}
\docType{data}
\title{A dataset of grids for producing spatial predictions of seabed sediment content in the Petrel sub-basin in Australia Exclusive Economic Zone}
\description{
This dataset contains 248675 rows of 6 variables including
longitude (long), latitude (lat), bathymetry (bathy), disttance to coast (dist), seabe relief (relief),
seabed slope (slope).
}
\usage{data("petrel")}
\format{
  A data frame with 248675 observations on the following 6 variables.
  \describe{
    \item{\code{long}}{a numeric vector, decimal degree}
    \item{\code{lat}}{a numeric vector, decimal degree}
    \item{\code{bathy}}{a numeric vector, meter bellow sea level}
    \item{\code{dist}}{a numeric vector, degree}
    \item{\code{relief}}{a numeric vector, meter}
    \item{\code{slope}}{a numeric vector, no unit}
  }
}
\details{For details, please check the reference.}
\source{Li, J., 2013. Predicting the spatial distribution of seabed gravel content using random forest, spatial interpolation methods and their hybrid methods, The International Congress on Modelling and Simulation (MODSIM) 2013: Adelaide, pp. 394-400.}
\keyword{datasets}
