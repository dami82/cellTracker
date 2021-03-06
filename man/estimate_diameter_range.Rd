% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTracker_core.R
\name{estimate_diameter_range}
\alias{estimate_diameter_range}
\title{Detect Paricle Diameters in a Numeric matrix}
\usage{
estimate_diameter_range(x, px.margin = 2, quantile.val = 0.99,
  plot = TRUE)
}
\arguments{
\item{x}{numeric matrix corresponding to a digital image}

\item{px.margin}{integer, number of pixels used as margin while searching/filtering for neighboring particles}

\item{quantile.val}{numeric, must be bigger than 0 and smaller than 1. 
Quantile for discriminating signal and background; only pixels with intensity higher than the corresponding 
quantile will count as signal while estimating particle diameters}

\item{plot}{logial, shall a histogram of the distribution of diameters be shown}
}
\value{
list including summary stats and data about the particles found in the image
}
\description{
Estimates the diameters of particles in a numeric matrix
}
\examples{
a <- cbind(c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1), 
           c(1, 1, 0, 0, 0, 0, 0, 0, 1, 1), 
           c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
           c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0), 
           c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0))
graphics::image(a)
b <- estimate_diameter_range(a)
print(b$estim.cell.num)
print(b$raw)

}
\references{
\url{https://www.data-pulse.com/dev_site/celltracker/}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
