% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R, R/summary.R
\name{print.summary.ddhazard}
\alias{print.summary.ddhazard}
\alias{summary.ddhazard}
\title{Summarizing Dynamic Hazard Models Fits}
\usage{
\method{print}{summary.ddhazard}(x, digits = getOption("digits"), ...)

\method{summary}{ddhazard}(object, var_indices = 1:ncol(object$state_vecs), max_print = 10, ...)
}
\arguments{
\item{x}{object returned from \code{summary.ddhazard}.}

\item{digits}{number of digits to print.}

\item{...}{not used.}

\item{object}{object returned from \code{\link{ddhazard}}.}

\item{var_indices}{variable indices to print for time-varying effects.}

\item{max_print}{maximum number of time points to print coefficients at.}
}
\description{
The \code{sd} printed for time-varying effects are point-wise standard deviations from the smoothed covariance matrices.
}
