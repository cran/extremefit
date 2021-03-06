% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PackageOxygen.R
\name{goftest}
\alias{goftest}
\title{Goodness of fit test statistics}
\usage{
goftest(object, ...)
}
\arguments{
\item{object}{model object.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
The form of the value returned by goftest depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.
}
\description{
goftest is a generic function whose application depends on the class of its argument.
}
\references{
Grama, I. and Spokoiny, V. (2008). Statistics of extremes by oracle estimation. Ann. of Statist., 36, 1619-1648.

Durrieu, G. and Grama, I. and Pham, Q. and Tricot, J.- M (2015). Nonparametric adaptive estimator of extreme conditional tail probabilities quantiles. Extremes, 18, 437-478.
}
\seealso{
\code{\link{goftest.hill.adapt}}, \code{\link{goftest.hill.ts}}
}
