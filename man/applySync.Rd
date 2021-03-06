% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/applySync.R
\name{applySync}
\alias{applySync}
\title{Apply sync model to toa matrix to obtain synced data}
\usage{
applySync(toa, hydros = "", sync_model)
}
\arguments{
\item{toa}{Object containing data to be synchronized. Typically a `data.table` as e.g. `ssu1$detections`, but can also be a matrix dim=(n_ping, n_hydo).}

\item{hydros}{data.table formatted as `ssu1$hydros`}

\item{sync_model}{Synchronization model obtained using `getSyncModel()`}
}
\description{
Apply sync model to toa matrix to obtain synced data
}
