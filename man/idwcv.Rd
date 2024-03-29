% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idwcv.R
\name{idwcv}
\alias{idwcv}
\title{Cross validation, n-fold for inverse distance weighting (IDW)}
\usage{
idwcv(longlat, trainy, cv.fold = 10, nmax = 12, idp = 2, predacc = "VEcv", ...)
}
\arguments{
\item{longlat}{a dataframe contains longitude and latitude of point
samples.}

\item{trainy}{a vector of response, must have length equal to the number of
rows in longlat.}

\item{cv.fold}{integer; number of folds in the cross-validation. if > 1,
then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
validation that is recommended.}

\item{nmax}{for a local predicting: the number of nearest observations that
should be used for a prediction or simulation, where nearest is
defined in terms of the space of the spatial locations. By default, 12
observations are used.}

\item{idp}{numeric; specify the inverse distance weighting power.}

\item{predacc}{can be either "VEcv" for vecv or "ALL" for all measures
in function pred.acc.}

\item{...}{other arguments passed on to gstat.}
}
\value{
A list with the following components:
 me, rme, mae, rmae, mse, rmse, rrmse, vecv and e1; or vecv only.
}
\description{
This function is a cross validation function for inverse
distance weighting.
}
\note{
This function is largely based on rfcv in randomForest and
some functions in library(gstat).
}
\examples{
\dontrun{
library(sp)
data(swmud)
data(petrel)

idwcv1 <- idwcv(swmud[, c(1,2)], swmud[, 3], nmax = 12, idp = 2)
idwcv1

n <- 20 # number of iterations, 60 to 100 is recommended.
VEcv <- NULL
for (i in 1:n) {
idwcv1 <- idwcv(petrel[, c(1,2)], petrel[, 3], nmax = 12, predacc = "VEcv")
VEcv [i] <- idwcv1
}
plot(VEcv ~ c(1:n), xlab = "Iteration for IDW", ylab = "VEcv (\%)")
points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
abline(h = mean(VEcv), col = 'blue', lwd=2)

n <- 20 # number of iterations, 60 to 100 is recommended.
measures <- NULL
for (i in 1:n) {
idwcv1 <- idwcv(swmud[, c(1,2)], swmud[, 3], predacc = "ALL")
measures <- rbind(measures, idwcv1$vecv)
}
plot(measures ~ c(1:n), xlab = "Iteration for IDW", ylab="VEcv (\%)")
points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
abline(h = mean(measures), col = 'blue', lwd = 2)
}

}
\references{
Li, J., 2013. Predictive Modelling Using Random Forest and Its
Hybrid Methods with Geostatistical Techniques in Marine Environmental
Geosciences, In: Christen, P., Kennedy, P., Liu, L., Ong, K.-L., Stranieri,
A., Zhao, Y. (Eds.), The proceedings of the Eleventh Australasian Data
Mining Conference (AusDM 2013), Canberra, Australia, 13-15 November 2013.
Conferences in Research and Practice in Information Technology, Vol. 146.

A. Liaw and M. Wiener (2002). Classification and Regression by
randomForest. R News 2(3), 18-22.

Pebesma, E.J., 2004. Multivariable geostatistics in S: the gstat
package. Computers & Geosciences, 30: 683-691.
}
\author{
Jin Li
}
