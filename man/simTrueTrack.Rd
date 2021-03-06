% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simTrack.R
\name{simTrueTrack}
\alias{simTrueTrack}
\title{Simulate a known movement track for subsequent estimation using YAPS}
\usage{
simTrueTrack(model = "rw", n, deltaTime = 1, D = NULL,
  shape = NULL, scale = NULL, addDielPattern = TRUE, ss = "rw",
  start_pos = NULL)
}
\arguments{
\item{model}{Movement model: 'rw': Two-dimension random walk  (X,Y)}

\item{n}{Number of steps in the simulated track}

\item{deltaTime}{Number of time units (seconds) between each location}

\item{D}{Diffusivity of the animal movement - only used if model='rw'}

\item{shape}{Shape of the Weibull distribution - only used when model='crw'.}

\item{scale}{Scale of the Weibull distribution - only used when model='crw'.}

\item{addDielPattern}{Adds a realistic(?) diel pattern to movement. Periods of both low and high movement}

\item{ss}{Simulations model for Speed of Sound - defaults to 'rw' = RW-model.}

\item{start_pos}{Specify the starting position of the track with c(x0, y0)}
}
\value{
Dataframe containing a simulated track
}
\description{
Produces a simulated regular time-spaced track following the specified movement model. Linear movement between consecutive observations is assumed.
The output contains x, y, time and sound speed at each simulated position.
}
