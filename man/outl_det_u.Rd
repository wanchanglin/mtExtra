% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_extra.R
\name{outl_det_u}
\alias{outl_det_u}
\title{Univariate outlier detection}
\usage{
outl_det_u(x, method = c("boxplot", "median", "mean"))
}
\arguments{
\item{x}{a numeric vector.}

\item{method}{method for univariate outlier detection. Only \code{boxplot},
\code{mean} and \code{median} are supported.}
}
\value{
a logical vector.
}
\description{
Perform outlier detection using univariate method.
}
\details{
\itemize{
\item \code{mean}: the absolute difference between that observation and the sample
mean is more than 2 times of SD.
\item \code{median}: the absolute difference between the observation and the sample
median is larger than 2 times of the Median Absolute Deviation divided
by 0.6745.
\item \code{boxplot}: either smaller than the 1st quartile minus 1.5 times of IQR,
or larger than the 3rd quartile plus 1.5 times of IQR.
}
}
\examples{
x <- c(2, 3, 4, 5, 6, 7, 8, 9, 50, 50)
outl_det_u(x, "boxplot")
}
\references{
Wilcox R R, Fundamentals of Modern Statistical Methods: Substantially
Improving Power and Accuracy, Springer 2010 (2nd edition), pages 31-35.
}
\seealso{
Other outlier detectors: 
\code{\link{outl_det_m}()}
}
\concept{outlier detectors}
