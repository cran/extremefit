\name{dataWind}
\alias{dataWind}
\docType{data}
\title{Wind speed for Brest (France)}
\description{The data frame provides the wind speed of Brest from 1976 to 2005.}
\usage{data("dataWind")}
\format{The data is the wind speed in meters per second (m/s) every day from 1976 to 2005.
\describe{
\item{\code{Year}}{The year of the measure.}
\item{\code{Month}}{The month of the measure.}
\item{\code{Day}}{the day of the measure.}
\item{\code{Speed}}{The wind speed in meters per second}
}}

\examples{

library(extremefit)
data("dataWind")
attach(dataWind)

pred <- NULL
for(m in 1:12){
  indices <- which(Month == m)
  X <- Speed[indices]*60*60/1000
  H <- hill.adapt(X)
  pred[m] <- predict(H, newdata = 100, type = "survival")$y
}
plot(pred, ylab = "Estimated survival probability", xlab = "Month")

}
