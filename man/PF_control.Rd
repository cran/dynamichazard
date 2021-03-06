% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PF.R
\name{PF_control}
\alias{PF_control}
\title{Auxiliary for Controlling Particle Fitting}
\usage{
PF_control(
  N_fw_n_bw = NULL,
  N_smooth = NULL,
  N_first = NULL,
  eps = 0.01,
  forward_backward_ESS_threshold = NULL,
  method = "AUX_normal_approx_w_cloud_mean",
  n_max = 25,
  n_threads = getOption("ddhazard_max_threads"),
  smoother = "Fearnhead_O_N",
  Q_tilde = NULL,
  est_a_0 = TRUE,
  N_smooth_final = N_smooth,
  nu = 0L,
  covar_fac = -1,
  ftol_rel = 1e-08,
  averaging_start = -1L,
  fix_seed = TRUE
)
}
\arguments{
\item{N_fw_n_bw}{number of particles to use in forward and backward filter.}

\item{N_smooth}{number of particles to use in particle smoother.}

\item{N_first}{number of particles to use at time \eqn{0} and time \eqn{d + 1}.}

\item{eps}{convergence threshold in EM method.}

\item{forward_backward_ESS_threshold}{required effective sample size to not re-sample in the particle filters.}

\item{method}{method for forward, backward and smoothing filter.}

\item{n_max}{maximum number of iterations of the EM algorithm.}

\item{n_threads}{maximum number threads to use in the computations.}

\item{smoother}{smoother to use.}

\item{Q_tilde}{covariance matrix of additional error term to add to the
proposal distributions. \code{NULL} implies no additional error term.}

\item{est_a_0}{\code{FALSE} if the starting value of the state model should
be fixed. Does not apply for \code{type = "VAR"}.}

\item{N_smooth_final}{number of particles to sample with replacement from
the smoothed particle cloud with \code{N_smooth} particles using the
particles' weights. This causes additional sampling error but decreases the
computation time in the M-step.}

\item{nu}{integer with degrees of freedom to use in the (multivariate)
t-distribution used as the proposal distribution. A (multivariate) normal
distribution is used if it is zero.}

\item{covar_fac}{factor to scale the covariance matrix with. Ignored if
the values is less than or equal to zero.}

\item{ftol_rel}{relative convergence tolerance of the mode objective in mode
approximation.}

\item{averaging_start}{index to start averaging. Values less then or equal
to zero yields no averaging.}

\item{fix_seed}{\code{TRUE} if the same seed should be used. E.g., in
\code{\link{PF_EM}} the same seed will be used in each iteration of the
E-step of the MCEM algorithm.}
}
\value{
A list with components named as the arguments.
}
\description{
Auxiliary for additional settings with \code{\link{PF_EM}}.
}
\details{
The \code{method} argument can take the following values

\itemize{
\item \code{bootstrap_filter} for a bootstrap filter.
\item \code{PF_normal_approx_w_cloud_mean} for a particle filter where a
Gaussian approximation is used using a Taylor
approximation made at the mean for the current particle given the mean of the
parent particles  and/or mean of the child particles.
\item \code{AUX_normal_approx_w_cloud_mean} for an auxiliary particle filter
version of \code{PF_normal_approx_w_cloud_mean}.
\item \code{PF_normal_approx_w_particles} for a filter similar to
\code{PF_normal_approx_w_cloud_mean} and differs by making a Taylor
approximation at a mean given each sampled parent and/or child particle.
\item \code{AUX_normal_approx_w_particles} for an auxiliary particle filter
version of \code{PF_normal_approx_w_particles}.
}

The \code{smoother} argument can take the following values
\itemize{
\item \code{Fearnhead_O_N} for the smoother in Fearnhead, Wyncoll, and Tawn
(2010).
\item \code{Brier_O_N_square} for the smoother in Briers, Doucet, and
Maskell (2010).
}
}
\references{
Gordon, N. J., Salmond, D. J., and Smith, A. F. (1993) Novel approach
to nonlinear/non-Gaussian Bayesian state estimation.
\emph{In IEE Proceedings F (Radar and Signal Processing)},
(Vol. 140, No. 2, pp. 107-113). IET Digital Library.

Pitt, M. K., and Shephard, N. (1999) Filtering via simulation: Auxiliary
particle filters. \emph{Journal of the American statistical association},
\strong{94(446)}, 590-599.

Fearnhead, P., Wyncoll, D., and Tawn, J. (2010) A sequential smoothing
algorithm with linear computational cost. \emph{Biometrika}, \strong{97(2)},
447-464.

Briers, M., Doucet, A., and Maskell, S. (2010) Smoothing algorithms for
state-space models.
\emph{Annals of the Institute of Statistical Mathematics}, \strong{62(1)},
61.
}
\seealso{
\code{\link{PF_EM}}
}
