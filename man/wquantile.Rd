% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PackageOxygen.R
\name{wquantile}
\alias{wquantile}
\title{Weighted quantile}
\usage{
wquantile(X, p, weights = rep(1, length(X)))
}
\arguments{
\item{X}{a vector of data.}

\item{p}{a vector of probabilities.}

\item{weights}{the weights assiociated to the vector \eqn{X}.}
}
\value{
A vector of quantile assiociated to the probabilities vector given in input.
}
\description{
Compute the weighted quantile of order p.
}
\details{
Give the weighted quantile for a given \eqn{p}
}
\examples{

X <- rpareto(10)
p <- seq(0.01, 0.99, 0.01)
plot(p, wquantile(X, p, rep(1,length(X))), type = "s")

}
