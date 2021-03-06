% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellTracker_core.R
\name{cell_tracker}
\alias{cell_tracker}
\title{Compute Cell Tracks}
\usage{
cell_tracker(tc_obj, lnoise = NULL, diameter = NULL,
  threshold = NULL, maxDisp = 25, memory_b = 0, goodenough = 0,
  threads = 1, show_plots = TRUE, verbose = TRUE)
}
\arguments{
\item{tc_obj}{a trackedCells object}

\item{lnoise}{numeric, lnoise parameter; can be NULL if optimize_params() has already been run}

\item{diameter}{numeric, diameter parameter; can be NULL if optimize_params() has already been run}

\item{threshold}{numeric, threshold parameter; can be NULL if optimize_params() has already been run}

\item{maxDisp}{numeric,  maximum displacement of a cell per time interval. 
When many cells are detected in each frame, small maxDisp values should be used.}

\item{memory_b}{numeric, memory_b parameter as used in the original track.m function.
In the current R implementation, only the value memory_b=0 is accepted}

\item{goodenough}{numeric, goodenough parameter as used in the original track.m function.
In the current R implementation, only the value goodenough=0 is accepted}

\item{threads}{integer, number of cores to use for parallelization}

\item{show_plots}{logical, shall cells detected in each frame of the image stack be visualized}

\item{verbose}{logical, shall info about the progress of the cell tracking job be printed}
}
\value{
a trackedCells object
}
\description{
Analyze Stacks, detect cells in each frame, and analyze cell tracks over time
}
\details{
The lnoise param is used to guide a lowpass blurring operation, while the lobject param is used
to guide a highpass background subtraction. The threshold param is used for a background correction following
the initial image convolution
\itemize{

  \item \strong{lnoise}: Characteristic lengthscale of noise in pixels. 
  Additive noise averaged over this length should vanish. May assume any positive floating value.
  May be also set to 0, in which case only the highpass "background subtraction" operation is performed.

  \item \strong{lobject} Integer length in pixels somewhat larger than a typical object. 
  Can also be set to 0, in which case only the lowpass "blurring" operation defined by lnoise is done 
  without the background subtraction defined by lobject
  
  \item \strong{threshold} Numeric. By default, after the convolution, any negative pixels are reset 
  to 0.  Threshold changes the threshhold for setting pixels to 0.  Positive values may be useful 
  for removing stray noise or small particles. 
  
}
}
\references{
\url{https://www.data-pulse.com/dev_site/celltracker/}
\url{https://www.mathworks.com/matlabcentral/fileexchange/60349-fasttracks}
}
\author{
Damiano Fantini, \email{damiano.fantini@gmail.com}
}
