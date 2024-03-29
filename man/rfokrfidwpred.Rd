% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rfokrfidwpred.R
\name{rfokrfidwpred}
\alias{rfokrfidwpred}
\title{Generate spatial predictions using the average of the hybrid method of
random forest and ordinary kriging and the hybrid method of random forest and
inverse distance weighting  (RFOKRFIDW)}
\usage{
rfokrfidwpred(
  longlat,
  trainx,
  trainy,
  longlatpredx,
  predx,
  mtry = function(p) max(1, floor(sqrt(p))),
  ntree = 500,
  idp = 2,
  nmaxok = 12,
  nmaxidw = 12,
  vgm.args = ("Sph"),
  block = 0,
  ...
)
}
\arguments{
\item{longlat}{a dataframe contains longitude and latitude of point
samples (i.e., trainx and trainy).}

\item{trainx}{a dataframe or matrix contains columns of predictive variables.}

\item{trainy}{a vector of response, must have length equal to the number of
rows in trainx.}

\item{longlatpredx}{a dataframe contains longitude and latitude of point
locations (i.e., the centres of grids) to be predicted.}

\item{predx}{a dataframe or matrix contains columns of predictive variables for
the grids to be predicted.}

\item{mtry}{a function of number of remaining predictor variables to use as
the mtry parameter in the randomForest call.}

\item{ntree}{number of trees to grow. This should not be set to too small a
number, to ensure that every input row gets predicted at least a few times.
By default, 500 is used.}

\item{idp}{numeric; specify the inverse distance weighting power.}

\item{nmaxok}{for local predicting: the number of nearest observations that
should be used for a prediction or simulation, where nearest is
defined in terms of the space of the spatial locations. By default, 12
observations are used for OK.}

\item{nmaxidw}{for local predicting: the number of nearest observations that
should be used for a prediction or simulation, where nearest is
defined in terms of the space of the spatial locations. By default, 12
observations are used for IDW.}

\item{vgm.args}{arguments for vgm, e.g. variogram model of response
variable and anisotropy parameters. see notes vgm in gstat for details.
By default, "Sph" is used.}

\item{block}{block size. see krige in gstat for details.}

\item{...}{other arguments passed on to randomForest or gstat.}
}
\value{
A dataframe of longitude, latitude, predictions and variances. The
variances are the same as the variances of rfokpred.
}
\description{
This function is to make spatial predictions using the average
of the hybrid method of random forest and ordinary kriging and the hybrid
method of random forest and inverse distance weighting  (RFOKRFIDW).
}
\note{
This function is largely based rfcv in randomForest.  When 'A zero or
negative range was fitted to variogram' occurs, to allow OK running, the
range was set to be positive by using min(vgm1$dist). In this case, caution
should be taken in applying this method, although sometimes it can still
outperform IDW and OK.
}
\examples{
\dontrun{
data(petrel)
data(petrel.grid)
rfokrfidwpred1 <- rfokrfidwpred(petrel[, c(1,2)], petrel[, c(1,2, 6:9)],
petrel[, 3], petrel.grid[, c(1,2)], petrel.grid, ntree = 500, idp = 2,
nmaxok = 12, nmaxidw = 12)
names(rfokrfidwpred1)
}

}
\references{
Liaw, A. and M. Wiener (2002). Classification and Regression by
randomForest. R News 2(3), 18-22.
}
\author{
Jin Li
}
