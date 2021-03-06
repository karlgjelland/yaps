% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/syncGetters.R
\name{getInpSync}
\alias{getInpSync}
\title{Get object inp for synchronization}
\usage{
getInpSync(sync_dat, max_epo_diff, min_hydros, time_keeper_idx,
  fixed_hydros_idx, n_offset_day, n_ss_day, keep_rate = 1,
  excl_self_detect = TRUE, lin_corr_coeffs = NA)
}
\arguments{
\item{sync_dat}{List containing data.tables with hydrophone information and detections. See e.g. `?ssu1` for example}

\item{max_epo_diff}{Sets the upper threshold for differences in TOA of sync tags. Best parameter value depends on burst rate of sync tags and how far apart the internal clocks of the hydros are prior to synchronization. A bit less than half of minimum sync tag burst rate is a good starting choice.}

\item{min_hydros}{Sets the lower threshold of how many hydrophones need to detect each sync tag ping in order to be included in the sync process. Should be as high as possible while observing that all hydrosphones are contributing. If too low, isolated hydrophones risk falling out completely. Future versions will work towards automising this.}

\item{time_keeper_idx}{Index of the hydrophone to use as time keeper. Could e.g. be the one with smallest overall clock-drift.}

\item{fixed_hydros_idx}{Vector of hydro idx's for all hydrophones where the position is assumed to be known with adequate accuracy and precission. Include as many as possible as fixed hydros to reduce overall computation time and reduce overall variability. As a bare minimum two hydros need to be fixed, but we strongly advice to use more than two.}

\item{n_offset_day}{Specifies the number of hydrophone specific quadratic polynomials to use per day. For PPM based systems, 1 or 2 is often adeqaute.}

\item{n_ss_day}{Specifies number of speed of sound to estimate per day. Future versions will enable use of logged water temperature instead. However, estimating SS gives an extra option for sanity-checking the final sync-model.}

\item{keep_rate}{Syncing large data sets can take a long time. However, there is typically an excess number of sync tag detections and a sub-sample is typically enough for good synchronization. This parameter specifies the proportion (0-1) of data to keep when sub-sampling.}

\item{excl_self_detect}{Logical whether to excluded detections of sync tags on the hydros they are co-located with. Sometimes self detections can introduce excessive residuals in the sync model in which case they should be excluded.}

\item{lin_corr_coeffs}{Matrix of coefficients used for pre-sync linear correction. dim(lin_corr_coeffs)=(#hydros, 2).}
}
\description{
Get object inp for synchronization
}
