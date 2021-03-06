% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simTrack.R
\name{simTelemetryTrack}
\alias{simTelemetryTrack}
\title{Simulate telemetry track based on known true track obtained using simTrueTrack}
\usage{
simTelemetryTrack(trueTrack, pingType, sbi_mean = NULL, sbi_sd = NULL,
  rbi_min = NULL, rbi_max = NULL)
}
\arguments{
\item{trueTrack}{Know track obtained using simTrueTrack}

\item{pingType}{Type of transmitter to simulate - either stable burst interval ('sbi'), random burst interval ('rbi') or random burst interval but where the random sequence is known a priori}

\item{sbi_mean, sbi_sd}{Mean and SD of burst interval when pingType = 'sbi'}

\item{rbi_min}{Minimum and maximum BI for random burst interval transmitters}

\item{rbi_max}{Minimum and maximum BI for random burst interval transmitters}
}
\value{
Data frame containing time of ping and true positions
}
\description{
Based on a known true track obtained using simTrueTrack, this function will give true positions at time-of-pings, which are also in the output. TOPs are determined by user-specified transmitter type.
Number of pings are determined automatically based on track length and transmitter specifications.
}
