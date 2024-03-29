% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_extra.R
\name{outl_det_m}
\alias{outl_det_m}
\title{Multivariate outlier detection}
\usage{
outl_det_m(x, method = "mcd", conf.level = 0.95)
}
\arguments{
\item{x}{a data matrix.}

\item{method}{methods for resistant estimation of multivariate location
and scatter. Only \code{mve}, \code{mcd} and \code{classical} are supported.}

\item{conf.level}{a confidential level.}
}
\value{
a logical vector.
}
\description{
Perform multivariate outlier detection.
}
\examples{
 
set.seed(134)
x <- cbind(rnorm(80), rnorm(80), rnorm(80))
y <- cbind(rnorm(10, 5, 1), rnorm(10, 5, 1), rnorm(10, 5, 1))
x <- rbind(x, y)
outl <- outl_det_m(x, method = "mcd", conf.level = 0.95)
}
\seealso{
\code{\link[=cov.rob]{cov.rob()}} for "Resistant Estimation of Multivariate Location
and Scatter"

Other outlier detectors: 
\code{\link{outl_det_u}()}
}
\concept{outlier detectors}
