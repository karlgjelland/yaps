% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepTmb.R
\name{getDatTmb}
\alias{getDatTmb}
\title{Get data for input to TMB}
\usage{
getDatTmb(hydros, toa, E_dist, n_ss, pingType, rbi_min, rbi_max,
  ss_data_what, ss_data, biTable, inp_params, z_vec)
}
\arguments{
\item{hydros}{Dataframe from simHydros() or Dataframe with columns hx and hy containing positions of the receivers. Translate the coordinates to get the grid centre close to (0;0).}

\item{toa}{TOA-matrix: matrix with receivers in rows and detections in columns. Make sure that the receivers are in the same order as in hydros, and that the matrix is very regular: one ping per column (inlude empty columns if a ping is not detected).}

\item{E_dist}{Which distribution to use in the model - "Gaus" = Gaussian, "Mixture" = mixture of Gaussian and t or "t" = pure t-distribution}

\item{n_ss}{Number of soundspeed estimates: one estimate per hour is usually enough}

\item{pingType}{Type of transmitter to simulate - either stable burst interval ('sbi'), random burst interval ('rbi') or random burst interval but where the random sequence is known a priori}

\item{rbi_min}{Minimum and maximum BI for random burst interval transmitters}

\item{rbi_max}{Minimum and maximum BI for random burst interval transmitters}

\item{ss_data_what}{What speed of sound (ss) data to be used. Default ss_data_what='est': ss is estimated by the model. Alternatively, if ss_data_what='data': ss_data must be provided and length(ss_data) == ncol(toa)}

\item{ss_data}{Vector of ss-data to be used if ss_data_what = 'est'. Otherwise ss_data <- 0 (default)}

\item{biTable}{Table of known burst intervals. Only used when pingType == "pbi". Default=NULL}

\item{inp_params}{Selection of parameters used to setup and run YAPS.}

\item{z_vec}{Vector of known depth values (positive real). Default=NULL is which case no 3D is assumed. Estimation of depth from detections is currently not supported.}
}
\value{
List for use in TMB.
}
\description{
Compile data for input to TMB.
}
