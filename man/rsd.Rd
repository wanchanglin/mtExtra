% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mt_extra.R
\name{rsd}
\alias{rsd}
\title{Relative standard deviation}
\usage{
rsd(x, na.rm = TRUE)
}
\arguments{
\item{x}{an vector, matrix or data frame.}

\item{na.rm}{a logical value indicating whether NA values should be stripped
before the computation proceeds.}
}
\value{
RSD value multiplied by 100.
}
\description{
Calculate Relative Standard Deviation(RSD). RSD is also known as the
coefficient of variation (CV)
}
\details{
Some notes:
\itemize{
\item The higher the CV, the greater the dispersion in the variable.
\item The CV is often presented as the given ratio multiplied by 100
\item Basically CV<10 is very good, 10-20 is good, 20-30 is acceptable,
and CV>30 is not acceptable.
}
}
\examples{
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

data(iris)
rsd(iris[, 1:4])

## group rsd
val <- iris \%>\% 
  group_by(Species) \%>\% 
	 group_modify(~ dat_summ(., method = rsd)) \%>\%
  pivot_longer(cols = !Species) \%>\% filter(!is.na(value))

ggplot(val, aes(x = value)) +
  geom_histogram(colour = "black", fill = "white") +
  facet_grid(Species ~ .)
ggplot(val, aes(x = Species, y = value)) + geom_boxplot()

## The plotting of missing value is similar like this.
}
